from time import time
import sys, traceback
import socket,select, os, gc
import cv2
import numpy as np
import cStringIO
import fcntl
import beast

P_MATCH_THRESH=0.9

def trace(frame, event, arg):
    print>>sys.stderr,"%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno)
    return trace

#sys.settrace(trace)

CONFIGFILE=sys.argv[1]
YEAR=float(sys.argv[2])

#set up server before we do anything else
server=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
server.bind(('127.0.0.1', 8010))
server.listen(5)
server.setblocking(0)

print "Loading config" 
beast.load_config(CONFIGFILE)
print "Loading hip_main.dat" 
S_DB=beast.star_db()
S_DB.load_catalog("hip_main.dat",YEAR)
print "Filtering stars" 
SQ_RESULTS=beast.star_query(S_DB)
SQ_RESULTS.kdmask_filter_catalog()
SQ_RESULTS.kdmask_uniform_density(beast.cvar.REQUIRED_STARS)
S_FILTERED=SQ_RESULTS.from_kdmask()
print "Generating DB" 
C_DB=beast.constellation_db(S_FILTERED,2+beast.cvar.DB_REDUNDANCY,0)
print "Ready"

#Note: SWIG's policy is to garbage collect objects created with
#constructors, but not objects created by returning from a function

class star_image:
	def __init__(self, imagefile,median_image):
		b_conf=[time(),beast.cvar.PIXSCALE,beast.cvar.BASE_FLUX]
		self.img_stars = beast.star_db()
		self.img_data = []
		self.img_const = None
		self.match=None
		self.db_stars=None
		self.match_from_lm=None
		self.db_stars_from_lm=None
		
		#Placeholders so that these don't get garbage collected by SWIG
		self.fov_db=None
		self.const_from_lm=None
		
		#TODO: improve memory efficiency
		if "://" in imagefile:
			import urllib
			img=cv2.imdecode(np.asarray(bytearray(urllib.urlopen(imagefile).read()), dtype="uint8"), cv2.IMREAD_COLOR)
		else:
			img=cv2.imread(imagefile)
			
		img=np.clip(img.astype(np.int16)-median_image,a_min=0,a_max=255).astype(np.uint8)
		img_grey = cv2.cvtColor(img,cv2.COLOR_RGB2GRAY)
		
		#removes areas of the image that don't meet our brightness threshold
		ret,thresh = cv2.threshold(img_grey,beast.cvar.THRESH_FACTOR*beast.cvar.IMAGE_VARIANCE,255,cv2.THRESH_BINARY)
		contours,heirachy = cv2.findContours(thresh,1,2);

		for c in contours:
			M = cv2.moments(c)
			if M['m00']>0:
				#this is how the x and y position are defined by cv2
				cx = M['m10']/M['m00']
				cy = M['m01']/M['m00']
				#see https://alyssaq.github.io/2015/computing-the-axes-or-orientation-of-a-blob/
				#for how to convert these into eigenvectors/values
				u20 = M["m20"]/M["m00"] - cx**2
				u02 = M["m02"]/M["m00"] - cy**2
				u11 = M["m11"]/M["m00"] - cx*cy
				#the center pixel is used as the approximation of the brightest pixel
				self.img_stars.add_star(cx-beast.cvar.IMG_X/2.0,(cy-beast.cvar.IMG_Y/2.0),float(cv2.getRectSubPix(img_grey,(1,1),(cx,cy))[0,0]),-1)
				self.img_data.append(b_conf+[cx,cy,u20,u02,u11]+cv2.getRectSubPix(img,(1,1),(cx,cy))[0,0].tolist())
				
		#self.img_const=beast.constellation_db(self.img_stars.copy(),beast.cvar.MAX_FALSE_STARS+2,1)
		self.img_const=beast.constellation_db(self.img_stars,beast.cvar.MAX_FALSE_STARS+2,1)
		
	
	def match_near(self,x,y,z,r):
		SQ_RESULTS.kdsearch(x,y,z,r,beast.cvar.THRESH_FACTOR*beast.cvar.IMAGE_VARIANCE)
		#estimate density for constellation generation
		C_DB.results.kdsearch(x,y,z,r,beast.cvar.THRESH_FACTOR*beast.cvar.IMAGE_VARIANCE)
		fov_stars=SQ_RESULTS.from_kdresults()#REE
		self.fov_db = beast.constellation_db(fov_stars,C_DB.results.kdresults_size,1)
		C_DB.results.clear_kdresults()
		SQ_RESULTS.clear_kdresults()
		
		near = beast.db_match(self.fov_db,self.img_const)
		if near.p_match>P_MATCH_THRESH:
			self.match = near
			self.db_stars = near.winner.from_match()
		
	def match_lis(self):
		lis=beast.db_match(C_DB,self.img_const)
		if lis.p_match>P_MATCH_THRESH:
			x=lis.winner.R11
			y=lis.winner.R12
			z=lis.winner.R13
			self.match_near(x,y,z,beast.cvar.MAXFOV/2)

	def match_rel(self,last_match):
		#make copy of stars from lastmatch
		img_stars_from_lm=last_match.img_stars.copy()
		#update positions to eci
		w=last_match.match.winner
		#convert the stars to ECI
		for i in range(img_stars_from_lm.map_size):
			s=img_stars_from_lm.get_star(i)
			x=s.x*w.R11+s.y*w.R21+s.z*w.R31
			y=s.x*w.R12+s.y*w.R22+s.z*w.R32
			z=s.x*w.R13+s.y*w.R23+s.z*w.R33
			s.x=x
			s.y=y
			s.z=z
		self.const_from_lm=beast.constellation_db(img_stars_from_lm,beast.cvar.MAX_FALSE_STARS+2,1)
		rel=beast.db_match(self.const_from_lm,self.img_const)
		if rel.p_match>P_MATCH_THRESH:
			self.match_from_lm = rel
			self.db_stars_from_lm = rel.winner.from_match()
				
	def print_match(self):
		db=self.db_stars
		im=self.img_stars
		if db is None:
			if self.db_stars_from_lm is None:
				#neither relative nor absolute matching could be used
				print ""
				return
			else:
				db=self.db_stars_from_lm
		assert(db.map_size==im.map_size)
		star_out=[]
		for i in range(db.map_size):
			s_im=im.get_star(i)
			s_db=db.get_star(i)
			if (s_db.id>=0):
				weight=1.0/(s_db.sigma_sq+s_im.sigma_sq)
				star_out.append(str(s_im.x)+','+str(s_im.y)+','+str(s_im.z)+','+str(s_db.x)+','+str(s_db.y)+','+str(s_db.z)+','+str(weight))
		print " ".join(star_out)

NONSTARS={}
NONSTAR_NEXT_ID=0
NONSTAR_DATAFILENAME="data"+str(time())+".txt"
NONSTAR_DATAFILE=open(NONSTAR_DATAFILENAME,"w")
class nonstar:
	def __init__(self,current_image,i,source):
		global NONSTARS,NONSTAR_NEXT_ID,NONSTAR_DATAFILENAME,NONSTAR_DATAFILE
		self.id=NONSTAR_NEXT_ID
		NONSTARS[self.id]=self
		current_image.img_stars.get_star(i).id=self.id
		NONSTAR_NEXT_ID+=1
		self.data=[]
		self.add_data(current_image,i,source)
		
	def add_data(self,current_image,i,source):
		s_im=current_image.img_stars.get_star(i)
		s_db_x=0.0
		s_db_y=0.0
		s_db_z=0.0
		w=None
		if (current_image.match != None and current_image.match.p_match>P_MATCH_THRESH):
			w=current_image.match.winner
		elif (current_image.match != None and current_image.match.p_match>P_MATCH_THRESH):
			w=current_image.match_from_lm.winner
		if w != None:
			#convert the stars to ECI
			s_db_x=s_im.x*w.R11+s_im.y*w.R21+s_im.z*w.R31
			s_db_y=s_im.x*w.R12+s_im.y*w.R22+s_im.z*w.R32
			s_db_z=s_im.x*w.R13+s_im.y*w.R23+s_im.z*w.R33
		self.data.append([source,s_im.x,s_im.y,s_im.z,s_db_x,s_db_y,s_db_z]+current_image.img_data[i])
	def write_data(self,fd):
		os.write(fd,str(self.id)+" " +str(len(self.data))+"\n")
		for i in self.data:
			s=[str(j) for j in i]
			os.write(fd," ".join(s)+"\n")
	def __del__(self):
		self.write_data(NONSTAR_DATAFILE.fileno())

def flush_nonstars():
	global NONSTARS,NONSTAR_NEXT_ID,NONSTAR_DATAFILENAME,NONSTAR_DATAFILE
	NONSTARS={}
	NONSTAR_NEXT_ID=0
	gc.collect()
	NONSTAR_DATAFILE.close()
	NONSTAR_DATAFILENAME="data"+str(time())+".txt"
	NONSTAR_DATAFILE=open(NONSTAR_DATAFILENAME,"w")
	
def update_nonstars(current_image,source):
	global NONSTARS,NONSTAR_NEXT_ID
	nonstars_next={}
	im=current_image.img_stars
	db=current_image.db_stars
	db_lm=current_image.db_stars_from_lm
	if (db!=None):
		assert(db.map_size==im.map_size)
	if (db_lm!=None):
		assert(db_lm.map_size==im.map_size)
		for i in range(im.map_size):
			im.get_star(i).id=db_lm.get_star(i).id
	for i in range(im.map_size):
		s_im=im.get_star(i)
		#is this a star? if so remove from nonstars
		if (db != None and db.get_star(i).id>=0):
			if (s_im.id in NONSTARS):
				del NONSTARS[s_im.id]
			s_im.id=-1
		#if it's already there, add the latest mesurement
		elif (s_im.id in NONSTARS):
			NONSTARS[s_im.id].add_data(current_image,i,source)
			nonstars_next[s_im.id]=NONSTARS[s_im.id]
		#otherwise add a new nonstar
		else:
			ns=nonstar(current_image,i,source)
			nonstars_next[ns.id]=ns
	NONSTARS=nonstars_next
	
	#wrap around to prevent integer overflow
	if (NONSTAR_NEXT_ID>2**30):
		flush_nonstars()



class star_camera:
	def __init__(self, median_file,source="RGB"):
		self.source=source
		self.current_image=None
		self.last_match=None
		self.median_image=cv2.imread(median_file)

	def solve_image(self,imagefile):
		starttime=time()
		self.current_image=star_image(imagefile,self.median_image)
		self.current_image.match_lis()
		if self.last_match is not None:
			self.current_image.match_rel(self.last_match)
		self.current_image.print_match()
		update_nonstars(self.current_image,self.source)
		if self.current_image.match is not None:
			self.last_match=self.current_image
		else:
			self.last_match=None
		print>>sys.stderr,"Time: "+str(time() - starttime)

#dummy for now
#TODO: add science data from IR cam
class science_camera:
	def __init__(self, median_file,source="IR"):
		self.source=source
		self.current_image=None
		self.last_match=None
		self.median_image=cv2.imread(median_file)
	def solve_image(self,imagefile):
		os.write(1,os.path.abspath(NONSTAR_DATAFILENAME))

rgb=star_camera(sys.argv[3])
ir=science_camera(sys.argv[3])

CONNECTIONS = {}
class connection:
	"""Tracks activity on a file descriptor and allows TCP read/writes"""
	def __init__(self, conn, epoll):
		"""
		Create connection to track file descriptor activity
		@note: Adds C{fd . self} to C{CONNECTIONS}
		@param conn: Any file object with the fileno() method
		@param epoll: File descriptor edge polling object
		"""
		self.conn=conn
		self.fd = self.conn.fileno()
		epoll.register(self.fd, select.EPOLLIN)
		self.epoll = epoll
		CONNECTIONS[self.fd] = self

	def read(self):
		"""
		Complete non-blocking read on file descriptor
		of an arbitrary amount of data
		@return: Entire read string
		@rtype: C{string}
		"""
		# need nonblocking for read
		fl = fcntl.fcntl(self.fd, fcntl.F_GETFL)
		fcntl.fcntl(self.fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)
		data = ''
		try:
			while True:
				lastlen=len(data)
				data += os.read(self.fd, 1024)
				if len(data)==lastlen:
					break
		except OSError, e:
			pass
			# error 11 means we have no more data to read
			if e.errno != 11:
				raise
		return data

	def write(self, data):
		"""
		Blocking read on file descriptor
		@param data: ASCII data to write
		@type data: C{string}
		"""
		if len(data) == 0: return
		if self.fd==0: #stdin
			os.write(1, data)
			return
		# need blocking IO for writing
		fl = fcntl.fcntl(self.fd, fcntl.F_GETFL)
		fcntl.fcntl(self.fd, fcntl.F_SETFL, fl & ~os.O_NONBLOCK)
		os.write(self.fd, data)

	def close(self):
		"""Close connection with file descriptor"""
		self.epoll.unregister(self.fd)
		self.conn.close()
		if CONNECTIONS[self.fd]:
			del CONNECTIONS[self.fd]

epoll = select.epoll()
epoll.register(server.fileno(), select.EPOLLIN)

try:
	connection(sys.stdin,epoll)
except IOError:
	pass
	
while True:
	events = epoll.poll()
	for fd, event_type in events:
		# Activity on the master socket means a new connection.
		if fd == server.fileno():
			conn, addr = server.accept()
			connection(conn, epoll)
		elif fd in CONNECTIONS:
			w = CONNECTIONS[fd]
			data = w.read()
			print data
			if len(data) > 0:
				stdout_redir = cStringIO.StringIO()
				stdout_old = sys.stdout
				sys.stdout = stdout_redir
				try:
					exec(data)
				except:
					traceback.print_exc(file=sys.stdout)
				sys.stdout = stdout_old
				data_out = stdout_redir.getvalue()
				print data_out
				w.write(data_out)
			else:
				w.close()

