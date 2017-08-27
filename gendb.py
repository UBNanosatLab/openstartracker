import numpy as np
from scipy import spatial
import math
import sys
import heapq
import itertools, operator
import config

#this python file generates constelations based on what stars might be within a simulated field of view
#the purpose of this is to make the database robust against cases where part of a constelation is cut off, but we still have enough stars to form another constelation
#takes in 1 commandline argument. this is our minimum fov/2

execfile(config.PROJECT_ROOT+"catalog_gen/calibration/calibration.txt")

def minfovradius():
	if DEG_X<DEG_Y:
		return DEG_X/2.0
	else:
		return DEG_Y/2.0

def maxfovradius():
	return math.sqrt(DEG_X**2+DEG_Y**2)

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


#load our star catalog, converting from id,ra,dec to x,y,z,id
def getstardb(filename=config.PROJECT_ROOT+"catalog_gen/catalog.dat"):
	stardb={}
	starfile = open(filename)
	for line in starfile.readlines():
		fields=line.split()
		HIP_ID=int(fields[0]);
		MAG=float(fields[1]);
		DEC=float(fields[2]);
		RA=float(fields[3]);
		X=float(fields[4]);
		Y=float(fields[5]);
		Z=float(fields[6]);
		MAX_90P=float(fields[7]);
		MIN_90P=float(fields[8]);
		UNRELIABLE=int(fields[9]);
		stardb[HIP_ID]=[HIP_ID,MAG,DEC,RA,X,Y,Z,MAX_90P,MIN_90P,UNRELIABLE]
	return stardb

stardb=getstardb()

#generate a list of coordinates which covers the sky completely
#with the constraint that the space between each point may be no more than fov/2
#the idea is we are simulating the process of pointing our camera at every possible region of sky
def getglobe():
	globe=[]
	step=minfovradius()/math.sqrt(2)
	decstep=180/math.floor(180/step)
	for dec in list(np.arange(-90+decstep/2,90,decstep)):
		rastep=360/math.floor(360/(step/math.cos(math.radians(dec))))
		for ra in list(np.arange(0,360,rastep)):
			x=math.cos(math.radians(ra))*math.cos(math.radians(dec))
			y=math.sin(math.radians(ra))*math.cos(math.radians(dec))
			z=math.sin(math.radians(dec))
			globe.append([x,y,z])
	return np.array(globe)


def checkcoverage(starlist):
	global ARC_ERR
	global stardb
	sd=np.array(stardb.values(),dtype=object)
	err=ARC_ERR*2./3600.
	xyz=np.array(sd[starlist,4:7].tolist(),dtype=float)
	fovxyz=getglobe()
	result=searchxyz(xyz,fovxyz,minfovradius())
	numgood=0
	for i in range(0,len(fovxyz)):
		if len(result[i])>1:
			numgood+=1
	print >>sys.stderr,"Database coverage: "+str(100.0*numgood/len(fovxyz)) + "% percent of the sky"

#radius is in degrees
def searchxyz(starxyz,points,radius=None):
	if (radius==None):
		radius=minfovradius()
	xyz = spatial.cKDTree(starxyz)
	pts = spatial.cKDTree(points)
	return pts.query_ball_tree(xyz,2*abs(math.sin(math.radians(radius)/2)))


def filterbrightness(mag=None):
	global stardb
	if (mag==None):
		mag=MIN_MAG
	sd=np.array(stardb.values(),dtype = object)
	for i in sd:
		if i[1]>mag:
			del stardb[i[0]]

def filterunreliable():
	global stardb
	sd=np.array(stardb.values(),dtype = object)
	for i in sd:
		if i[9]:
			del stardb[i[0]]

def filterdoublestars(r=None):
	global stardb
	if (r==None):
		r=PIXSCALE*2*PSF_RADIUS
	sd=np.array(stardb.values(),dtype = object)
	xyz=np.array(sd[:,4:7].tolist(),dtype=float)
	for i in searchxyz(xyz,xyz,r/3600.):
		for j in range(1,len(i)):
			k=sd[i[j]][0]
			if k in stardb:
				del stardb[k]

def sort_uniq(sequence):
    return [i for i in itertools.imap(operator.itemgetter(0),itertools.groupby(sorted(sequence)))]

#return a set of brightest stars with average density nstars per field of view
def uniformstarlist(nstars=4):
	global stardb
	sd=np.array(stardb.values(),dtype = object)
	xyz=np.array(sd[:,4:7].tolist(),dtype=float)
	starlist=[]
	for i in searchxyz(xyz,xyz):
		starlist+=heapq.nsmallest(nstars,i,key=lambda j: sd[j][1])
	return sort_uniq(starlist)

#xyz_dist takes xyz of two points
#returns distance in degrees
def xyz_dist(xyz1,xyz2):
	return math.degrees(math.asin(np.linalg.norm(np.cross(xyz1,xyz2))))

#only do this part if we were run as a python script
if __name__ == '__main__':
	filterbrightness()
	filterdoublestars()
	filterunreliable()

	#print stars
	if (len(sys.argv)>1):
		f= open(sys.argv[1], 'w')
	else:
		f= open("/dev/stdout", 'w')
	for i in stardb.values():
		print >>f, i[0],i[1],i[4],i[5],i[6]
	f.close()

	starlist=uniformstarlist()
	#print info to the user about sky coverage
	checkcoverage(starlist)

	#print constellations
	if (len(sys.argv)>2):
		f= open(sys.argv[2], 'w')
	else:
		f= open("/dev/stdout", 'w')

	sd=np.array(stardb.values(),dtype = object)
	sd[:,0]=np.array(range(0,len(sd)))
	uni_sd=sd[starlist]
	xyz=np.array(sd[:,4:7].tolist(),dtype=float)
	uni_xyz=np.array(uni_sd[:,4:7].tolist(),dtype=float)
	#for each star in the uniform starlist, get all nearby stars in the database
	ns=searchxyz(xyz,uni_xyz,maxfovradius()/2)
	#get nearby stars uniform star database
	uni_ns=searchxyz(uni_xyz,uni_xyz,maxfovradius()/2)
	
	STARTABLE=0
	NUMCONST=0
	NUMSTARS=len(sd)
	for i in range(0,len(uni_ns)):
		for j in uni_ns[i]:
			if (j>i):
				temp=list(set(ns[i])&set(ns[j]))
				print >>f, 3600*xyz_dist(uni_xyz[i],uni_xyz[j]),uni_sd[i,0],uni_sd[j,0],len(temp)
				print >>f, " ".join([str(s) for s in temp])
				STARTABLE+=len(temp)
				NUMCONST+=1
	f.close()
	#print db info
	if (len(sys.argv)>3):
		f= open(sys.argv[3], 'w')
	else:
		f= open("/dev/stdout", 'w')
	print >>f, "STARTABLE="+str(STARTABLE)
	print >>f, "NUMCONST="+str(NUMCONST)
	print >>f, "NUMSTARS="+str(NUMSTARS)
	print >>f, "PARAM="+str(int(2+3600*maxfovradius()/(2*ARC_ERR)))
	
