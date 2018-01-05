from time import time
import sys
import socket,select, os
import cv2
import numpy as np
import cStringIO
import fcntl
import beast

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
		self.img_stars = beast.star_db()
		self.img_flux = []
		self.img_const = None
		self.match=None
		self.match_rel_result=None
		self.db_stars=None
		self.db_stars_rel=None
		
		#Placeholders so that these don't get garbage collected by SWIG
		self.fov_db=None
		self.lm_const_rel=None
		
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
		if len(contours)<2:
			return
		for c in contours:
			M = cv2.moments(c)
			if M['m00']>0:
				#this is how the x and y position are defined by cv2
				cx = M['m10']/M['m00']
				cy = M['m01']/M['m00']
				#the center pixel is used as the approximation of the brightest pixel
				self.img_stars.add_star(cx-beast.cvar.IMG_X/2.0,(cy-beast.cvar.IMG_Y/2.0),float(cv2.getRectSubPix(img_grey,(1,1),(cx,cy))[0,0]),self.img_stars.map_size)
				self.img_flux.append(cv2.getRectSubPix(img,(1,1),(cx,cy))[0,0].tolist())
				
		self.img_const=beast.constellation_db(self.img_stars.copy(),beast.cvar.MAX_FALSE_STARS+2,1)
		
	
	def match_near(self,x,y,z,r):
		SQ_RESULTS.kdsearch(x,y,z,r,beast.cvar.THRESH_FACTOR*beast.cvar.IMAGE_VARIANCE)
		#estimate density for constellation generation
		C_DB.results.kdsearch(x,y,z,r,beast.cvar.THRESH_FACTOR*beast.cvar.IMAGE_VARIANCE)
		fov_stars=SQ_RESULTS.from_kdresults()#REE
		self.fov_db = beast.constellation_db(fov_stars,C_DB.results.kdresults_size,1)
		C_DB.results.clear_kdresults()
		SQ_RESULTS.clear_kdresults()
		
		near = beast.db_match(self.fov_db,self.img_const)
		if near.p_match>0.9:
			self.match = near
			self.db_stars = near.winner.from_match()
		
	def match_lis(self):
		lis=beast.db_match(C_DB,self.img_const)
		if lis.p_match>0.9:
			x=lis.winner.R11
			y=lis.winner.R12
			z=lis.winner.R13
			self.match_near(x,y,z,beast.cvar.MAXFOV/2)

	def match_rel(self,last_match):
		#make copy of stars from lastmatch
		lm_db_stars=last_match.img_stars.copy()
		#update positions to eci
		w=last_match.match.winner
		#convert the stars to ECI
		for i in range(lm_db_stars.map_size):
			s=lm_db_stars.get_star(i)
			x=s.x*w.R11+s.y*w.R21+s.z*w.R31
			y=s.x*w.R12+s.y*w.R22+s.z*w.R32
			z=s.x*w.R13+s.y*w.R23+s.z*w.R33
			s.x=x
			s.y=y
			s.z=z
		self.lm_const_rel=beast.constellation_db(lm_db_stars,beast.cvar.MAX_FALSE_STARS+2,1)
		rel=beast.db_match(self.lm_const_rel,self.img_const)
		if rel.p_match>0.9:
			self.match_rel_result = rel
			self.db_stars_rel = rel.winner.from_match()
			if self.match is None:
				self.match = self.match_rel_result
				self.db_stars = self.db_stars_rel
				
	def print_match(self):
		if self.match is None:
			print ""
			return
		assert(self.db_stars.map_size==self.img_stars.map_size)
		star_out=[]
		for i in range(self.db_stars.map_size):
			im=self.img_stars.get_star(i)
			db=self.db_stars.get_star(i)
			if (db.id>=0):
				weight=1.0/(db.sigma_sq+im.sigma_sq)
				star_out.append(str(im.x)+','+str(im.y)+','+str(im.z)+','+str(db.x)+','+str(db.y)+','+str(db.z)+','+str(weight))
		print " ".join(star_out)

class star_camera:
	def __init__(self, median_file):
		self.last_match=None
		self.median_image=cv2.imread(median_file)
		
	def solve_image(self,imagefile):
		starttime=time()
		self.current_image=star_image(imagefile,self.median_image)
		self.current_image.match_lis()
		if self.last_match is not None:
			self.current_image.match_rel(self.last_match)
		self.current_image.print_match()
		if self.current_image.match is not None:
			self.last_match=self.current_image
		else:
			self.last_match=None
		print>>sys.stderr,"Time: "+str(time() - starttime)

rgb=star_camera(sys.argv[3])
ir=star_camera(sys.argv[3])


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
			if len(data) > 0:
				stdout_redir = cStringIO.StringIO()
				stdout_old = sys.stdout
				sys.stdout = stdout_redir
				try:
					exec(data)
				except Exception as e:
					print str(e)
				sys.stdout = stdout_old
				data_out = stdout_redir.getvalue()
				w.write(data_out)
			else:
				w.close()
