#ifndef STARS_H
#define STARS_H

#include "config.h"

#include <assert.h> //assert()
#include <limits.h> //INT_MAX
#include <algorithm> //sort, nth_element
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <iterator>     // std::next


//Sadness inducing bit twiddling stuff

/**
* @brief hash function which interleaves values so that nearby stars will have nearby hashes
* this approach also  also allows you to truncate the hash when less precision is needed
* 
* for the magic numbers used for interleaving, see:
* https://stackoverflow.com/questions/1024754/how-to-compute-a-3d-morton-number-interleave-the-bits-of-3-ints
*/
size_t xyz_hash(float x, float y,float z){
	#define INTERLEAVE_THREE(X)\
		if (X > 0x1fffff) X=0x1fffff;\
		if (X < 0) X=0;\
		X = (X | X <<32) & 0x001f00000000ffff;\
		X = (X | X <<16) & 0x001f0000ff0000ff;\
		X = (X | X << 8) & 0x100f00f00f00f00f;\
		X = (X | X << 4) & 0x10c30c30c30c30c3;\
		X = (X | X << 2) & 0x1249249249249249;
	
	long h0=0x100000*(x+1.0);
	long h1=0x100000*(y+1.0);
	long h2=0x100000*(z+1.0);
	INTERLEAVE_THREE(h0)
	INTERLEAVE_THREE(h1)
	INTERLEAVE_THREE(h2)
	return h0<<2|h1<<1|h2;
	
	#undef INTERLEAVE_THREE
}
/**
* @brief returns a mask which can be used to clamp the hash to a specified precision
* @param radians mask a range of this size, rounded up to the nearest power of two
*/
size_t xyz_hash_mask(float radians) {
	long h0 = 0x100000*radians;
	if (h0 > 0xfffff) h0=0xfffff;
	if (h0 < 0) h0=0;
	size_t mask=-1;
	for (;h0;h0>>=1) mask<<=3;
	return mask;
}

struct star {
	float x;
	float y;
	float z;
	float flux;
	/**user defined id (ie hipparcos id, -1)*/
	int id;
	float px;
	float py;
	int unreliable;
	/** how many stars were inserted before this one? */
	int star_idx;
	float sigma_sq;
	size_t hash_val;
	/**
	* @brief add star from catalog 
	*
	* @param x ECI  'coming out of image'
	* @param y ECI 
	* @param z ECI z
	* @param flux Pixel brightness
	* @param id User defined id
	*/
	star(){};
	star(float x_, float y_, float z_, float flux_, int id_) {
		x=x_;
		y=y_;
		z=z_;
		flux=flux_;
		id=id_;
		
		px=y/(x*PIXX_TANGENT);
		py=z/(x*PIXY_TANGENT);

		unreliable=0;
		star_idx=-1;
		sigma_sq=POS_VARIANCE;
		hash_val=xyz_hash(x,y,z);
	}
	
	/**
	* @brief add star from image
	*
	* @param px Pixel x minus camera center
	* @param py Pixel y minus camera center
	* @param flux Pixel brightness
	* @param id  User defined id
	*/
	star(float px_, float py_, float flux_, int id_) {
		px=px_;
		py=py_;
		flux=flux_;
		id=id_;

		float j=(PIXX_TANGENT*px); /* j=(y/x) */
		float k=(PIXY_TANGENT*py); /* k=z/x */
		x=1./sqrt(j*j+k*k+1);
		y=j*x;
		z=k*x;

		unreliable=0;
		star_idx=-1;
		sigma_sq=IMAGE_VARIANCE/flux;
		hash_val=xyz_hash(x,y,z);
	}

	#define OP operator==
	bool OP(const star& s) const {return hash_val==s.hash_val;}
	#undef OP
	
	/**
	* @brief numerically stable method to calculate distance between stars
	* @param s Star
	* @return Angular seperation in arcsec
	*/
	#define OP operator*
	float OP(const star& s) const {
		float a=x*s.y - s.x*y;
		float b=x*s.z - s.x*z;
		float c=y*s.z - s.y*z;
		return (3600*180.0/PI)*asin(sqrt(a*a+b*b+c*c));
	}
	#undef OP
	#define OP operator size_t
	OP() const {return hash_val;}
	#undef OP
	/**
	* @brief Print debug info
	* @param s Label
	*/
	void DBG_(const char *s) {
		DBG_PRINT("%s\t",s);
		DBG_PRINT("x=%f ", x);
		DBG_PRINT("y=%f ", y);
		DBG_PRINT("z=%f ", z);
		DBG_PRINT("flux=%f ", flux);
		DBG_PRINT("star_idx=%d ", star_idx);
		DBG_PRINT("id=%d ", id);
		DBG_PRINT("unreliable=%d ", unreliable);
		DBG_PRINT("sigma_sq=%f ", sigma_sq);
		DBG_PRINT("px=%f ", px);
		DBG_PRINT("py=%f\n", py);
	}
};

bool star_gt_x(const star &s1, const star &s2) {return s1.x > s2.x;}
bool star_gt_y(const star &s1, const star &s2) {return s1.y > s2.y;}
bool star_gt_z(const star &s1, const star &s2) {return s1.z > s2.z;}
bool star_gt_flux(const star &s1, const star &s2) {return s1.flux > s2.flux;}
bool star_lt_x(const star &s1, const star &s2) {return s1.x < s2.x;}
bool star_lt_y(const star &s1, const star &s2) {return s1.y < s2.y;}
bool star_lt_z(const star &s1, const star &s2) {return s1.z < s2.z;}
bool star_lt_flux(const star &s1, const star &s2) {return s1.flux < s2.flux;}

struct star_db {
private:
	std::unordered_map<size_t,star> hash_map;
	std::set<size_t> hash_set;
	std::map<int,size_t> star_idx_map;
	std::multimap<float,size_t> flux_map;
	size_t sz;

	/**
	 * @brief Transcribe a portion of the db between first and last
	 * 
	 * @param first an iterator to a map. key can be anything, but value must be a star hash 
	 * @param last copy up to but not including the element pointed to by this iterator
	 */
	template<class T> star_db* copy(T first,T last) {
		star_db* s = new star_db;
		for (;first!=last;first++) (*s)+=hash_map[first->second];
		return s;
	}
public:
	//TODO
	float max_variance;
	star_db() {
		DBG_STAR_DB_COUNT++;
		DBG_PRINT("DBG_STAR_DB_COUNT++ %d\n",DBG_STAR_DB_COUNT);
		max_variance=0.0;
		sz=0;
	}
	~star_db() {
		DBG_STAR_DB_COUNT--;
		DBG_PRINT("DBG_STAR_DB_COUNT-- %d\n",DBG_STAR_DB_COUNT);
	}
	size_t size() {return sz;}
	///Philosophically inspired by python sets
	#define OP operator+=
	star_db* OP(const star* s) {
		if (count(s)==0) {
			if (max_variance<s->sigma_sq) max_variance=s->sigma_sq;
			star temp=s[0];
			temp.star_idx=size();
			hash_map[temp]=temp;
			hash_set.emplace(temp);
			flux_map.emplace(temp.flux,temp);
			star_idx_map.emplace(temp.star_idx,temp);
			sz++;
		}
		return this;
	}
	//TODO - faster to use insert?
	star_db* OP(star_db* s) {
		for (size_t i=0;i<s->size();i++) (*this)+=s->get_star(i);
		return this;
	}
	star_db* OP(const star& s) { return *this+=&s;}
	#undef OP
	#define OP operator-
	star_db* OP(star_db* s) {
		star_db* r = new star_db;
		for (auto it = star_idx_map.cbegin(); it != star_idx_map.cend(); ++it) {
			if (s->hash_map.count(it->second)==0) *r+=hash_map[it->second];
		}
		return r;
	}
	#undef OP
	#define OP operator&
	star_db* OP(star_db* s) {
		star_db* r = new star_db;
		for (auto it = star_idx_map.cbegin(); it != star_idx_map.cend(); ++it) {
			if (s->hash_map.count(it->second)>0) *r+=hash_map[it->second];
		}
		return r;
	}
	#undef OP
	/**
	 * @brief returns stars in the order they were added
	 * 
	 * @param idx the index of the star 
	 */
	star* get_star(int idx) {return size()>0?&(hash_map[star_idx_map[idx]]):NULL;}
	/**
	* @brief make a copy of the star db
	*/
	star_db* copy() {return copy(star_idx_map.cbegin(),star_idx_map.cend());}
	/**
	* @brief make a copy of the n brightest elements in the star db
	*/
	star_db* copy_n_brightest(size_t n) {return copy(flux_map.crbegin(),std::next(flux_map.crbegin(),std::min(n,size())));}
	/**
	* @brief return stars in the bounding volume around the specified star
	* @param r minimum radius of the bounding volume (radians).
	* 
	*/
	//TODO: group together by hash_lb,ub, n_brightest_search (maybe return a list?) 
	template<class T> void search(T &hs, float x,float y,float z, float r, float min_flux) {
		size_t mask=xyz_hash_mask(r);
		for (int8_t dx=-1;dx<=1;dx++) for (int8_t dy=-1;dy<=1;dy++) for (int8_t dz=-1;dz<=1;dz++) {
			size_t h=xyz_hash(x+dx*r,y+dy*r,z+dz*r);
			auto first = hash_set.lower_bound(h&mask);
			auto last = hash_set.upper_bound(h|(~mask));
			star *s=&hash_map[*first];
			for (;first!=last;first++) if (min_flux <= s->flux && s->x*s->x+s->y*s->y+s->z*s->z<=r*r) hs.insert(*first);
		}
	}


	/**
	* @brief Load stars from hip_main.dat
	*
	* @param catalog path to hip_main.dat
	* @param year Update star positions to the specified year
	*/
	void load_catalog(const char* catalog, float year) {
		FILE *stream = fopen(catalog, "r");
		if (stream == NULL) exit(EXIT_FAILURE);
		max_variance=POS_VARIANCE;
		
		ssize_t read;
		char *line = NULL;
		size_t len = 0;
		
		char* hip_record[78];
		while ((read = getline(&line, &len, stream)) != -1){
				float yeardiff=year-1991.25;
				hip_record[0]=strtok(line,"|");
				for (size_t j=1;j<sizeof(hip_record)/sizeof(hip_record[0]);j++) hip_record[j] = strtok(NULL,"|");
			
				float MAG=atof(hip_record[5]);
				float DEC=yeardiff*atof(hip_record[13])/3600000.0 + atof(hip_record[9]);
				float cosdec=cos(PI*DEC/180.0);
				float RA=yeardiff*atof(hip_record[12])/(cosdec*3600000.0) + atof(hip_record[8]);

				star s = star(cos(PI*RA/180.0)*cosdec,sin(PI*RA/180.0)*cosdec,sin(PI*DEC/180.0),BASE_FLUX*powf(10.0,-MAG/2.5),atoi(hip_record[1]));
				s.unreliable=((atoi(hip_record[29])==0||atoi(hip_record[29])==1)&&atoi(hip_record[6])!=3)?0:1;
				(*this)+=s;
		}
		free(line);
		fclose(stream);
		
	}

	size_t count(const star* s) {return hash_map.count(*s);}
	size_t count(star_db* s) {
		size_t n=0;
		for (size_t i=0;i<s->size();i++) n+=count(s->get_star(i));
		return n;
	}
	

	void DBG_(const char *s) {
		DBG_PRINT("%s\n",s);
		DBG_PRINT("star_db at %lu contains %lu elements\n",(size_t)this,size());
		DBG_PRINT("max_variance=%f\n",max_variance);
		for (size_t i=0; i<size(); i++) {
			DBG_PRINT("%lu:\t",i);
			get_star(i)->DBG_("star");
		}
	}
};

struct star_fov {
private:
	int *mask;
	star_db *stars;
	int *collision;
	int collision_size;
	float db_max_variance;
	float *s_px;
	float *s_py;
	/**
	* @brief Get the id of the best match to the specified coordinates 
	* Used to resolve collisions where the coordinates falls into the region of overlap between two stars
	* Adds a bit of complexity in exchange for being able to break ties at
	* the subpixel level, which can sometimes make a difference 
	* 
	* @param id - the id from the image map
	* Any id <-1 is interpreted as an index in the collision buffer (starts at -2)
	*
	* @param px Pixel x minus camera center
	* @param py Pixel y minus camera center
	*
	* @return The id of whichever star is the best match to the coordinates in question 
	*/
	int resolve_id(int id,float px,float py) {
		if (id>=-1) return id;
		//Any id <-1 is interpreted as an index in the collision buffer (starts at -2)
		id=-id;
		int id1=resolve_id(collision[id-2],px,py);
		int id2=resolve_id(collision[id-1],px,py);
		return (get_score(id1,px,py)>get_score(id2,px,py))?id1:id2;
	}
public:
	
	/**
	* @brief TODO
	*
	* @param id
	* @param px
	* @param py
	*
	* @return 
	*/
	float get_score(int id,float px,float py) {
		float sigma_sq,maxdist_sq;
		sigma_sq=stars->max_variance+db_max_variance;
		maxdist_sq=-sigma_sq*(log(sigma_sq)+MATCH_VALUE);
		
		float dx=px-s_px[id];
		float dy=py-s_py[id];
		if (dx<-0.5) dx+=1.0;/* use whichever corner of the pixel gives the best score */
		if (dy<-0.5) dy+=1.0;
		return (maxdist_sq-(dx*dx+dy*dy))/(2*sigma_sq);
	}
			
	/**
	* @brief TODO
	*
	* @param id
	* @param px
	* @param py
	* @param sigma_sq
	* @param maxdist_sq
	*
	* @return 
	*/
	float get_score(int id,float px,float py,float sigma_sq,float maxdist_sq) {
		float dx=px-s_px[id];
		float dy=py-s_py[id];
		if (dx<-0.5) dx+=1.0;/* use whichever corner of the pixel gives the best score */
		if (dy<-0.5) dy+=1.0;
		return (maxdist_sq-(dx*dx+dy*dy))/(2*sigma_sq);
	}
	
	int get_id(float px, float py) {
		int nx=(int)(px+IMG_X/2.0f);
		if (nx==-1) nx++;
		else if (nx==IMG_X) nx--;
	
		int ny=(int)(py+IMG_Y/2.0f);
		if (ny==-1) ny++;
		else if (ny==IMG_Y) ny--;
		int id=-1;
		if (nx>=0&&nx<IMG_X&&ny>=0&&ny<IMG_Y) id=mask[nx+ny*IMG_X];
		return resolve_id(id,px,py);
	}
	/**
	* @brief TODO
	* @param s
	* @param db_max_variance_
	*/
	star_fov(star_db* s, float db_max_variance_) {
		DBG_STAR_FOV_COUNT++;
		DBG_PRINT("DBG_STAR_FOV_COUNT++ %d\n",DBG_STAR_FOV_COUNT);
		db_max_variance=db_max_variance_;
		stars=s;
		collision=NULL;
		collision_size=0;
		s_px=(float*)malloc(stars->size()*sizeof(s_px[0]));
		s_py=(float*)malloc(stars->size()*sizeof(s_py[0]));
		mask=(int*)malloc(IMG_X*IMG_Y*sizeof(mask[0]));
		memset(mask, -1, IMG_X*IMG_Y*sizeof(mask[0]));
		/* generate image mask */
		for (size_t id=0;id<stars->size();id++){
			/* assume the dimmest possible star since we dont know the brightness of the other image */
			float sigma_sq,maxdist_sq;
			sigma_sq=stars->max_variance+db_max_variance;
			maxdist_sq=-sigma_sq*(log(sigma_sq)+MATCH_VALUE);
			float maxdist=sqrt(maxdist_sq);
			s_px[id]=stars->get_star(id)->px;
			s_py[id]=stars->get_star(id)->py;
			
			int xmin=s_px[id]-maxdist-1;
			int xmax=s_px[id]+maxdist+1;
			int ymin=s_py[id]-maxdist-1;
			int ymax=s_py[id]+maxdist+1;
			
			if(xmax>IMG_X/2) xmax=IMG_X/2;
			if(xmin<-IMG_X/2)xmin=-IMG_X/2;
			if(ymax>IMG_Y/2) ymax=IMG_Y/2;
			if(ymin<-IMG_Y/2)ymin=-IMG_Y/2;
			for(int i=xmin;i<xmax;i++) for (int j=ymin;j<ymax;j++) {

				float score=get_score(id, i,j, sigma_sq, maxdist_sq);
				
				int x=i+IMG_X/2;
				int y=j+IMG_Y/2;
			
				if (score>0) {
					/* has this pixel already been assigned to a different star? */
					int id2=mask[x+y*IMG_X];
					if (id2!=-1){
						collision_size+=2;
						mask[x+y*IMG_X]=-collision_size;
						collision=(int*)realloc(collision,collision_size*sizeof(collision[0]));
						collision[collision_size-2]=id;
						collision[collision_size-1]=id2;
					} else {
						mask[x+y*IMG_X]=id;
					}
				}
			}
		}
	}
	~star_fov() {
		DBG_STAR_FOV_COUNT--;
		DBG_PRINT("DBG_STAR_FOV_COUNT-- %d\n",DBG_STAR_FOV_COUNT);
		free(s_px);
		free(s_py);
		free(mask);
		free(collision);
	}
};

struct star_query {
	star *map;
	size_t map_size;
	size_t *kdresults;
private:
	uint8_t kdsorted;
	
	size_t kdresults_size;
	int8_t *kdmask;
	size_t kdresults_maxsize;
	star_db *stars;
		 
	 /** You may be looking at the most compact kd-tree in existence
	 *  It does not use pointers or indexes or leaf nodes or any extra memory
	 *  Instead the list is kdsorted in place using std::nth_element()
	 *  which is standard c++ implementation of quickselect.
	 * 
	 *  References:
	 *  Numerical Recipies (ISBN: 9780521884075)
	 *  https://en.wikipedia.org/wiki/Quickselect
	 *  https://stackoverflow.com/questions/17021379/balancing-kd-tree-which-approach-is-more-efficient
	 */
	
	#define KDSORT_NEXT(A,B)\
		int mid=(min+max)/2;\
		if (min+1<max) {\
			std::nth_element(map+min,map+mid,map+max,A);\
			if (mid-min>KDBUCKET_SIZE) B(min,mid);\
			else std::sort(map+min, map+mid,star_gt_flux);\
			if (max-(mid+1)>KDBUCKET_SIZE) B(mid+1,max);\
			else std::sort(map+(mid+1), map+max,star_gt_flux);\
		}
	void kdsort_x(int min, int max) {KDSORT_NEXT(star_lt_x,kdsort_y)}
	void kdsort_y(int min, int max) {KDSORT_NEXT(star_lt_y,kdsort_z)}
	void kdsort_z(int min, int max) {KDSORT_NEXT(star_lt_z,kdsort_x)}
	#undef KDSORT_NEXT
public:
	/**
	* @brief  TODO
	* @param s
	*/
	star_query(star_db *s) {
		DBG_STAR_QUERY_COUNT++;
		DBG_PRINT("DBG_STAR_QUERY_COUNT++ %d\n",DBG_STAR_QUERY_COUNT);
		stars=s;
		map_size=stars->size();
		kdsorted=0;
		kdresults_size=map_size;
		kdresults_maxsize=INT_MAX;

		kdmask=(int8_t*)malloc((map_size+1)*sizeof(kdmask[0]));
		kdresults=(size_t*)malloc((map_size+1)*sizeof(kdresults[0]));
		map=(star*)malloc(map_size*sizeof(map[0]));
		for (size_t i=0;i<map_size;i++){
			map[i]=stars->get_star(i)[0];
			kdresults[i]=i;
		}
		reset_kdmask();
	}
	~star_query() {
		DBG_STAR_QUERY_COUNT--;
		DBG_PRINT("DBG_STAR_QUERY_COUNT-- %d\n",DBG_STAR_QUERY_COUNT);
		free(map);
		free(kdresults);
		free(kdmask);
	}
	uint8_t is_kdsorted() {return kdsorted;}
	/**
	* @brief kdsort the list in question.
	*/
	void kdsort() {
		if (kdsorted==0) {
			kdsort_x(0,map_size);
			kdsorted=1;
		}
	}
	void sort() {
		std::sort(map, map+map_size,star_gt_flux);
		kdsorted=0;
	}
	size_t r_size() {return kdresults_size;}
	int8_t get_kdmask(size_t i) {return kdmask[i];}
	/**
	* @brief Clears kdmask, but does not reset kdresults. Slow.
	*/
	void reset_kdmask() {
		memset(kdmask,0,sizeof(kdmask[0])*stars->size());
	}
	
	/**
	* @brief Rewind kdresults, and set the corresponding kdmask to zero
	*/
	void clear_kdresults() {
		while (kdresults_size>0) {
			kdresults_size--;
			kdmask[kdresults[kdresults_size]]=0;
		}
	}

	/**
	* @brief TODO
	*
	* @param idx
	* @param x
	* @param y
	* @param z
	* @param r
	* @param min_flux
	*/
	void kdcheck(int idx, float x, float y, float z, float r, float min_flux){
		x-=map[idx].x;
		y-=map[idx].y;
		z-=map[idx].z;
		if (x-r <= 0 && 0 <= x+r && y-r <= 0 && 0 <= y+r &&z-r <= 0 && 0 <= z+r && min_flux <= map[idx].flux && kdmask[idx] == 0 && x*x+y*y+z*z<=r*r) {
			kdmask[idx]=1;
			/* Insertion sort into list from brightest to dimmest.
			 * 
			 * Note: Different sorting algorithms can result in different
			 * stars being selected in the case where there are two candidates of
			 * equal brightness for the last star. This has no effect on match quality
			 */
			 
			int n=kdresults_size++;
			
			float sm_flux=map[idx].flux;
			for (;n>0&&sm_flux>map[kdresults[n-1]].flux;n--) kdresults[n]=kdresults[n-1];
			kdresults[n]=idx;
			//if we go over the maximum, bump the dimmest star from the results
			if (kdresults_size>kdresults_maxsize) {
				kdresults_size=kdresults_maxsize;
				kdmask[kdresults[kdresults_size]]=0;
			}
		}
	}
	
	/**
	* @brief		search map for points within r pixels of x,y,z 
	* @details		put all results found into kdresults (sorted by brightness), mask via kdmask
	* @param x		cos(deg2rad*RA)*cos(deg2rad*DEC)
	* @param y		sin(deg2rad*RA)*cos(deg2rad*DEC)
	* @param z		sin(deg2rad*DEC)
	* @param r		search distance (pixels)
	* @param min_flux	minimum pixel brightness
	* @param min		start of bounding box
	* @param max		end of bounding box
	* @param dim		start dimension
	*/
	void kdsearch(float x, float y, float z, float r, float min_flux, int min, int max, int dim) {
		kdsort();
		float r_deg=r/3600.0;
		float r_rad=r_deg*PI/180.0;
		if (dim==0) kdsearch_x(x, y, z, 2*fabs(sin(r_rad/2.0)), min_flux,min,max);
		else if (dim==1) kdsearch_y(x, y, z, 2*fabs(sin(r_rad/2.0)), min_flux,min,max);
		else if (dim==2) kdsearch_z(x, y, z, 2*fabs(sin(r_rad/2.0)), min_flux,min,max);
		//std::unordered_set<size_t> s;
		//stars->search(s,x,y,z,r_rad,min_flux);
	}
	void kdsearch(float x, float y, float z, float r, float min_flux) {
		kdsearch(x, y, z, r, min_flux,0,stars->size(),0);
	}
	//use seperate functions for each diminsion so that the compiler can unroll the recursion
	#define KDSEARCH_NEXT(A,B,C,D)\
		int mid=(min+max)/2;\
		if (min<mid &&A <= map[mid].B) {\
			if (mid-min>KDBUCKET_SIZE) D(x,y,z,r,min_flux,min,mid);\
			else for (int i=min;i<mid&&min_flux<=map[i].flux;i++)kdcheck(i,x,y,z,r,min_flux);\
		}\
		if (mid<max) kdcheck(mid,x,y,z,r,min_flux);\
		if (kdresults_size==kdresults_maxsize) min_flux=map[kdresults[kdresults_size-1]].flux;\
		if (mid+1<max &&map[mid].B <= C) {\
			if (max-(mid+1)>KDBUCKET_SIZE) D(x,y,z,r,min_flux,mid+1,max);\
			else for (int i=mid+1;i<max&&min_flux<=map[i].flux;i++)kdcheck(i,x,y,z,r,min_flux);\
		}
	void kdsearch_x(const float x, const float y, const float z, const float r, float min_flux, int min, int max) {KDSEARCH_NEXT(x-r,x,x+r,kdsearch_y)}
	void kdsearch_y(const float x, const float y, const float z, const float r, float min_flux, int min, int max) {KDSEARCH_NEXT(y-r,y,y+r,kdsearch_z)}
	void kdsearch_z(const float x, const float y, const float z, const float r, float min_flux, int min, int max) {KDSEARCH_NEXT(z-r,z,z+r,kdsearch_x)}
	#undef KDSEARCH_NEXT

	/**
	* @brief set mask for stars that are too bright,too close, highly variable, or unreliable
	*/
	void kdmask_filter_catalog() {
		for (size_t i=0;i<stars->size();i++) {
			int8_t lastmask=kdmask[i];
			kdsearch(map[i].x,map[i].y,map[i].z,DOUBLE_STAR_PX*PIXSCALE,THRESH_FACTOR*IMAGE_VARIANCE);
			//TODO uncomment for real stars
			//if (kdresults_size>1||lastmask || map[i].flux<THRESH_FACTOR*IMAGE_VARIANCE||map[i].unreliable>0) {
			if (kdresults_size>1||lastmask || map[i].flux<THRESH_FACTOR*IMAGE_VARIANCE) {
				kdmask[i]=1;
				kdresults_size=0;
			} else {
				clear_kdresults();
			}
		}
	}
	/**
	* @brief Masks the dimmest stars in each area to produce a map with uniform density 
	* @param min_stars_per_fov Don't mask anything which could result in less than this many stars per field of view
	*/
	void kdmask_uniform_density(int min_stars_per_fov) {
		std::unordered_set<int> uniform_set;
		int kdresults_maxsize_old=kdresults_maxsize;
		kdresults_maxsize=min_stars_per_fov;
		for (size_t i=0;i<stars->size();i++) if (kdmask[i]==0) {
			kdsearch(map[i].x,map[i].y,map[i].z,MINFOV/2,THRESH_FACTOR*IMAGE_VARIANCE);
			for (size_t j=0;j<kdresults_size;j++) uniform_set.insert(kdresults[j]);
			clear_kdresults();
		}
		for (size_t i=0;i<stars->size();i++) kdmask[i]=1;
		std::unordered_set<int>::iterator it = uniform_set.begin();
		for (size_t i=0; i<uniform_set.size();i++,it++) kdmask[*it]=0;
		kdresults_maxsize=kdresults_maxsize_old;
	}
	/**
	* @brief Filter stardb based on mask
	* @return A new stardb containing only the stars which are not masked
	*/
	star_db* from_kdmask() {
		star_db* rd=new star_db;
		rd->max_variance=stars->max_variance;
		for (size_t i=0;i<stars->size();i++){
			if (kdmask[i]==0) *rd+=map[i];
		}
		return rd;
	}
	/**
	* @brief TODO
	* @return 
	*/
	star_db* from_kdresults() {
		star_db* rd=new star_db;
		rd->max_variance=stars->max_variance;
		for (size_t i=0;i<kdresults_size;i++){
			*rd+=map[kdresults[i]];
		}
		return rd;
	}
	
	void DBG_(const char *s) {
		DBG_PRINT("%s\n",s);
		DBG_PRINT("kdsorted=%d\n",kdsorted);
		DBG_PRINT("kdmask at %lu\n",(size_t)kdmask);
		DBG_PRINT("kdresults at %lu\n",(size_t)kdresults);
		DBG_PRINT("kdresults_size=%lu\n",kdresults_size);
		DBG_PRINT("kdresults_maxsize=%lu\n",kdresults_maxsize);
		if (kdresults_size>0){
			int i=0;
			DBG_PRINT("kdmask[%d]=%d\n",i,kdmask[i]);
			DBG_PRINT("kdresults[%d]=%lu\n",i,kdresults[i]);
			map[kdresults[i]].DBG_("STARS");
			DBG_PRINT(".\n.\n");
			i=kdresults_size-1;
			DBG_PRINT("kdmask[%d]=%d\n",i,kdmask[i]);
			DBG_PRINT("kdresults[%d]=%lu\n",i,kdresults[i]);
			map[kdresults[i]].DBG_("STARS");
		}
	}
};
#endif
