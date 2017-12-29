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
					 
			 
def get_body2ECI(RA,DEC,ORIENTATION):
	result = np.array([[1,0,0],[0,1,0],[0,0,1]])
	result = np.dot(result, rotation_matrix([0,0,1], math.radians(RA)))
	result = np.dot(result, rotation_matrix([0,1,0], math.radians(-DEC)))
	result = np.dot(result, rotation_matrix([1,0,0], math.radians(-ORIENTATION)))
	result = np.transpose(result)
	return result

def get_angles(body2ECI):
	RA=np.degrees(np.arctan2(body2ECI[0,1],body2ECI[0,0]))
	DEC=np.degrees(np.arcsin(body2ECI[0,2]))
	ORIENTATION=np.degrees(-np.arctan2(body2ECI[1,2],body2ECI[2,2]))
	return (RA,DEC,ORIENTATION)

def xyz_points(image_stars_info):
    """
    Converts the (x,y,magnitude) values to x,y,z points
        Input:
            image_stars_info: list of tuples in the form (x,y,magnitude)
        Output:
            star_points: a list of tuples in the form (x,y,z)
    """
    star_points = []
    for px,px,mag in image_stars_info:
        #formula for converting from x,y, magnitude to x,y,z
        j=2*np.tan((IMG_X*PIXSCALE/3600)*np.pi/(180*2))*px/IMG_X #j=(x/z)
        k=2*np.tan((IMG_Y*PIXSCALE/3600)*np.pi/(180*2))*py/IMG_Y #k=y/z
        x=1./np.sqrt(j*j+k*k+1)
        y=-j*x
        z=-k*x
        star_points.append([x,y,z])
    return star_points
    
def rigid_transform_3D(A, B, weight=[]):
	"""
	Takes in two matrices of points and finds the attitude matrix needed to
	transform one onto the other

	Input:
		A: nx3 matrix - x,y,z in body frame
		B: nx3 matrix - x,y,z in eci
		Note: the "n" dimension of both matrices must match
	Output:
		attitude_matrix: returned as a numpy matrix
	"""
	assert len(A) == len(B)
	if (len(weight) == 0):
		weight=np.array([1]*len(A))
	# dot is matrix multiplication for array
	H =  np.dot(np.transpose(A)*weight,B)
	
	#calculate attitude matrix
	#from http://malcolmdshuster.com/FC_MarkleyMortari_Girdwood_1999_AAS.pdf
	U, S, Vt = np.linalg.svd(H)
	U[:,2]*=np.linalg.det(U)*np.linalg.det(Vt)
	
	body2ECI = np.dot(U,Vt)
	body2ECI_RA_DEC_ORI(body2ECI)
	return body2ECI

def visualize(image,contours):
	"""
	Purely used to visualize the axis work. Draws the axes on the stars in
	bright lime green.
		Input:
			image: the image to have the axes drawn on its objects
			contours: list of contours in the object
	"""
	image = cv2.cvtColor(image,cv2.COLOR_GRAY2RGB)
	for c in contours:
		moments = cv2.moments(c)
		if(moments["m00"]<=0):
			continue;
		#Calculate the eigenvectors and eigenvalues of the centroid from the moments
		#definition of covariance matrix from image moments
		cy = moments['m01']/moments['m00']
		cx = moments['m10']/moments['m00']
		u20 = moments["m20"]/moments["m00"] - cx**2
		u02 = moments["m02"]/moments["m00"] - cy**2
		u11 = moments["m11"]/moments["m00"] - cx*cy
		eig_val,eig_vec = np.linalg.eig(np.matrix([[u20,u11],[u11,u20]]))

		eig_vec1 = np.matrix.tolist(eig_vec[0])[0]
		eig_vec2 = np.matrix.tolist(eig_vec[1])[0]

		# 4 = constant of proportionality
		axis_length = [4*np.sqrt(x) for x in eig_val]

		cy = moments['m01']/moments['m00']
		cx = moments['m10']/moments['m00']
		try:
			x11=(int(cx-eig_vec1[0]*axis_length[1]/2),int(cy-eig_vec1[1]*axis_length[1]/2))
			x12=(int(cx+eig_vec1[0]*axis_length[1]/2),int(cy+eig_vec1[1]*axis_length[1]/2))
			x21=(int(cx-eig_vec2[0]*axis_length[0]/2),int(cy-eig_vec2[1]*axis_length[0]/2))
			x22=(int(cx+eig_vec2[0]*axis_length[0]/2),int(cy+eig_vec2[1]*axis_length[0]/2))
			
			cv2.line(image,x11,x12,(0,0,255))
			cv2.line(image,x21,x22,(0,0,255))
		except ValueError:
			print >>sys.stderr, "Single pixel star"
	return image


def extract_stars(img,configfile,outputimg=0):
	"""
	Takes in a an image opencv image array (median already subtracted)

	Args:
		input_file: the image
	Returns:
		An array containing tuples of star info. Tuple format :(x,y,brightest value in star)

	"""
	
	cfg={}
	execfile(configfile,cfg)

	img=img.clip(0)
	img=img.astype(np.uint8)
	img = cv2.cvtColor(img,cv2.COLOR_RGB2GRAY)
	img = cv2.GaussianBlur(img,(3,3),0)
	#removes areas of the image that don't meet our brightness threshold
	ret,thresh = cv2.threshold(img,THRESH_FACTOR*IMAGE_VARIANCE,IMAGE_MAX,cv2.THRESH_BINARY)
	contours,heirachy = cv2.findContours(thresh,1,2);
	stars = []
	for c in contours:
		M = cv2.moments(c)
		if M['m00']>0:

			#this is how the x and y position are defined by cv2
			cx = M['m10']/M['m00']
			cy = M['m01']/M['m00']
			#the center pixel is used as the approximation of the brightest pixel
			#in the star in order to speed up the function by a factor of 10 and still
			#have beast return matches for constellations
			stars.append((cx,cy,cv2.getRectSubPix(img,(1,1),(cx,cy))))

	#use astrometry calibration data to correct for image distortion
	#see http://docs.astropy.org/en/stable/api/astropy.wcs.WCS.html
	results=np.array(stars)
	results[:,0]=results[:,0]-cfg['IMG_X']/2
	results[:,1]=results[:,1]-cfg['IMG_Y']/2
	if (outputimg==1):
		return (results,visualize(img,contours))
	else:
		return results

#radius is in degrees
def searchxyz(starxyz,points,radius):
	xyz = spatial.cKDTree(starxyz)
	pts = spatial.cKDTree(points)
	return pts.query_ball_tree(xyz,2*abs(math.sin(math.radians(radius)/2)))
	
#xyz_dist takes xyz of two points
#returns distance in degrees
def xyz_dist(xyz1,xyz2):
	return math.degrees(math.asin(np.linalg.norm(np.cross(xyz1,xyz2))))
