#ifndef STARS_H
#define STARS_H

#include "calibration.h"

#include <assert.h> //assert()
#include <limits.h> //INT_MAX
#include <algorithm> //sort, nth_element

struct star {
	float x;
	float y;
	float z;
	float photons;
	int32_t star_idx;
	int32_t id;
	int32_t unreliable;

	float sigma_sq;
	float px;
	float py;
	/* numerically stable method to calculate distance between stars */
	float dist_arcsec(const star& s) const {
		float a=x*s.y - s.x*y;
		float b=x*s.z - s.x*z;
		float c=y*s.z - s.y*z;
		return (3600*180.0/PI)*asin(sqrt(a*a+b*b+c*c));
	}
};

bool star_gt_x(const star &s1, const star &s2) {return s1.x > s2.x;}
bool star_gt_y(const star &s1, const star &s2) {return s1.y > s2.y;}
bool star_gt_z(const star &s1, const star &s2) {return s1.z > s2.z;}
bool star_gt_photons(const star &s1, const star &s2) {return s1.photons > s2.photons;}
bool star_lt_x(const star &s1, const star &s2) {return s1.x < s2.x;}
bool star_lt_y(const star &s1, const star &s2) {return s1.y < s2.y;}
bool star_lt_z(const star &s1, const star &s2) {return s1.z < s2.z;}
bool star_lt_photons(const star &s1, const star &s2) {return s1.photons < s2.photons;}



struct star_db {
	star *map;
	int map_size;
	float max_variance;
	int kdsorted;
	
	star_db() {
		map=NULL;
		map_size=0;
		kdsorted=0;

		
		max_variance=0.0;
	}
	~star_db() {
		free(map);
	}
	
	void add_star(float x, float y, float z, float mag, int id) {
		assert(kdsorted==0);
		int n=map_size++;
		if (n%16==0) map=(star*)realloc(map,(map_size+16)*sizeof(map[0]));

		map[n].x=x;
		map[n].y=y;
		map[n].z=z;
		map[n].id=id;
		map[n].unreliable=0;
		map[n].photons=PHOTONS*powf(10.0,-mag/2.5);/* TODO: change this to pixel value */
		map[n].px=y/(x*PIXX_TANGENT);
		map[n].py=z/(x*PIXY_TANGENT);
		map[n].sigma_sq=POS_VARIANCE;
		map[n].star_idx=n;
	}
	void add_star(float px, float py, float mag) {
		assert(kdsorted==0);
		float j=PIXX_TANGENT*px; /* j=(y/x) */
		float k=PIXY_TANGENT*py; /* k=z/x */
		float x=1./sqrt(j*j+k*k+1);
		float y=j*x;
		float z=k*x;
		int n=map_size;
		add_star(x,y,z,mag,-1);
		map[n].sigma_sq=IMAGE_VARIANCE/map[n].photons;
		if (max_variance<map[n].sigma_sq) max_variance=map[n].sigma_sq;
	}
	
	//TODO: should be part of config struct
	void load_catalog() {
		FILE *stream = fopen("catalog.dat", "r");
		if (stream == NULL) exit(EXIT_FAILURE);
		map_size=0;
		while(!feof(stream)) if(fgetc(stream) == '\n') map_size++;
		rewind(stream);
		max_variance=POS_VARIANCE;
		
		map = (star*)realloc(map,map_size*sizeof(map[0]));
		ssize_t read;
		char *line = NULL;
		size_t len = 0;
		
		for(int i=0;i<map_size;i++){
			if ((read = getline(&line, &len, stream)) != -1) {
				map[i].id=atoi(strtok(line," "));
				map[i].star_idx=i;
				float mag=atof(strtok(NULL," "));
				map[i].photons=PHOTONS*powf(10.0,-mag/2.5);
				
				map[i].x=atof(strtok(NULL," "));
				map[i].y=atof(strtok(NULL," "));
				map[i].z=atof(strtok(NULL," "));
				map[i].unreliable=atoi(strtok(NULL," "));
				
				map[i].sigma_sq=POS_VARIANCE;
			}
		}
		free(line);
		fclose(stream);
		
	}
	
	int* get_img_mask(float db_max_variance) {
		int *img_mask=(int*)malloc(IMG_X*IMG_Y*sizeof(img_mask[0]));
		memset(img_mask, -1, IMG_X*IMG_Y*sizeof(img_mask[0]));
		/* generate image mask */
		for (int id=0;id<map_size;id++){
			/* assume the dimmest possible star since we dont know the brightness of the other image */
			float sigma_sq=map[id].sigma_sq+db_max_variance;
			float maxdist_sq=-sigma_sq*(log(sigma_sq)+MATCH_VALUE);
			float maxdist=sqrt(maxdist_sq);
			int xmin=map[id].px-maxdist-1;
			int xmax=map[id].px+maxdist+1;
			int ymin=map[id].py-maxdist-1;
			int ymax=map[id].py+maxdist+1;
			
			if(xmax>IMG_X/2) xmax=IMG_X/2;
			if(xmin<-IMG_X/2)xmin=-IMG_X/2;
			if(ymax>IMG_Y/2) ymax=IMG_Y/2;
			if(ymin<-IMG_Y/2)ymin=-IMG_Y/2;
			for(int i=xmin;i<xmax;i++) for (int j=ymin;j<ymax;j++) {
				float a=((float)i-map[id].px);
				if (a<-0.5) a+=1.0;/* use whichever corner of the pixel gives the best score */
				float b=((float)j-map[id].py);
				if (b<-0.5) b+=1.0;
				
				int x=i+IMG_X/2;
				int y=j+IMG_Y/2;
				float score=(maxdist_sq-(a*a+b*b))/(2*sigma_sq);
			
				if (score>0) {
					/* has this pixel already been assigned to a different star? */
					int id2=img_mask[x+y*IMG_X];
					if (id2!=-1){
						float sigma_sq2=map[id2].sigma_sq+db_max_variance;
						float maxdist_sq2=-sigma_sq2*(log(sigma_sq2)+MATCH_VALUE);
						float px2=map[id2].px;
						float py2=map[id2].py;
						float a2=((float)x-px2-IMG_X/2);
						if (a2<-0.5) a2+=1.0;/* use whichever corner of the pixel gives the best score */
						float b2=((float)y-py2-IMG_Y/2);
						if (b2<-0.5) b2+=1.0;
						float score2 = (maxdist_sq2-(a2*a2+b2*b2))/(2*sigma_sq2);
						if (score>score2){
							img_mask[x+y*IMG_X]=id;
						}
					} else {
						img_mask[x+y*IMG_X]=id;
					}
				}
			}
		}
		return img_mask;
	}
};

struct star_query {
	int8_t *kdmask;

	int *kdresults;
	int kdresults_size;
	int kdresults_maxsize;
	
	int kdbucket_size;
	
	star_db *stars;
	star_query(star_db *s) {
		stars=s;
		
		kdresults=(int*)malloc(stars->map_size*sizeof(kdresults[0]));
		kdmask=(int8_t*)malloc(stars->map_size*sizeof(kdmask[0]));
		kdresults_size=0;
		kdresults_maxsize=INT_MAX;
		
		//TODO: put this in config file, come up with a way
		//to automatically figure out the optimal value
		kdbucket_size=(DEG_X*DEG_Y*3.5);
		
		if (stars->kdsorted==0) kdsort();
		
		reset_kdresults();
		reset_kdmask();
	}
	~star_query() {
		free(kdresults);
		free(kdmask);
	}
	/* You may be looking at the most compact kd-tree in existence
	 * It does not use pointers or indexes or leaf nodes or any extra memory
	 * Instead the list is kdsorted in place using std::nth_element()
	 * which is standard c++ implementation of quickselect.
	 * 
	 * References:
	 * Numerical Recipies (ISBN: 9780521884075)
	 * https://en.wikipedia.org/wiki/Quickselect
	 * https://stackoverflow.com/questions/17021379/balancing-kd-tree-which-approach-is-more-efficient
	 * 
	 */

	#define KDSORT_NEXT(A,B)\
		int mid=(min+max)/2;\
		if (min+1<max) {\
			std::nth_element(stars->map+min,stars->map+mid,stars->map+max,A);\
			if (mid-min>kdbucket_size) B(min,mid);\
			else std::sort(stars->map+min, stars->map+mid,star_gt_photons);\
			if (max-(mid+1)>kdbucket_size) B(mid+1,max);\
			else std::sort(stars->map+(mid+1), stars->map+max,star_gt_photons);\
		}

	void kdsort() {kdsort_x(0,stars->map_size);}
	void kdsort_x(int min, int max) {KDSORT_NEXT(star_lt_x,kdsort_y)}
	void kdsort_y(int min, int max) {KDSORT_NEXT(star_lt_y,kdsort_z)}
	void kdsort_z(int min, int max) {KDSORT_NEXT(star_lt_z,kdsort_x)}
	#undef KDSORT_NEXT
	 
	/* Reset kdresults without clearing kdmask
	 * Can be used to do clever things like reusing a mask generated
	 * by a previous search */
	void reset_kdresults() {kdresults_size=0;}
	
	/* Clears kdmask, but does not reset kdresults. Slow.*/
	void reset_kdmask() {memset(kdmask,0,sizeof(kdmask[0])*stars->map_size);}
	
	/* Undoes the effects of all kdsearches since the last time
	 * reset_kdresults() was called. Much faster than reset_kdmask()*/
	 
	void undo_kdsearch() {
		while (kdresults_size>0) {
			kdresults_size--;
			kdmask[kdresults[kdresults_size]]--;
		}
	}
	 
	void kdcheck(int idx, float x, float y, float z, float r, float min_photons){
		x-=stars->map[idx].x;
		y-=stars->map[idx].y;
		z-=stars->map[idx].z;
		if (x-r <= 0 && 0 <= x+r)
		if (y-r <= 0 && 0 <= y+r)
		if (z-r <= 0 && 0 <= z+r)
		if (min_photons <= stars->map[idx].photons)
		if (kdmask[idx] <= 0)
		if (x*x+y*y+z*z<=r*r) {
			kdmask[idx]=1;
			/* Insertion sort into list from brightest to dimmest.
			 * 
			 * Note: Different sorting algorithms can result in different
			 * stars being selected in the case where there are two candidates of
			 * equal brightness for the last star. Do not be alarmed -
			 * they all work fine
			 */
			 
			int n=kdresults_size++;
			for (;n>0&&stars->map[idx].photons>stars->map[kdresults[n-1]].photons;n--) kdresults[n]=kdresults[n-1];
			kdresults[n]=idx;
			//if we go over the maximum, bump the dimmest star from the results
			if (kdresults_size>kdresults_maxsize) {
				kdresults_size=kdresults_maxsize;
				kdmask[kdresults[kdresults_size]]--;
			}
		}
	}
	/**
	 * @fn		kdsearch(float x, float y, float z, float r, float min_photons)
	 * @param	x	cos(deg2rad*RA)*cos(deg2rad*DEC)
	 * @param	y	sin(deg2rad*RA)*cos(deg2rad*DEC)
	 * @param	z	sin(deg2rad*DEC)
	 * @param	r	search distance (pixels)
	 * @param	min	start of bounding box
	 * @param	max	end of bounding box
	 * @param	dim	start dimension
	 * @brief	search map for points within r pixels of x,y,z 
	 * @details put all results found into kdresults (sorted by brightness), mask via kdmask
	 */
	void kdsearch(float x, float y, float z, float r, float min_photons, int min, int max, int dim) {
		float r_deg=r*PIXSCALE/3600.0;
		float r_rad=r_deg*PI/180.0;
		if (dim==0) kdsearch_x(x, y, z, 2*fabs(sin(r_rad/2.0)), min_photons,min,max);
		else if (dim==1) kdsearch_y(x, y, z, 2*fabs(sin(r_rad/2.0)), min_photons,min,max);
		else if (dim==2) kdsearch_z(x, y, z, 2*fabs(sin(r_rad/2.0)), min_photons,min,max);
	}
	void kdsearch(float x, float y, float z, float r, float min_photons) {kdsearch(x, y, z, r, min_photons,0,stars->map_size,0);}
	
	//use seperate functions for each diminsion so that the compiler can unroll the recursion
	#define KDSEARCH_NEXT(A,B,C,D)\
		int mid=(min+max)/2;\
		if (min<mid &&A <= stars->map[mid].B) {\
			if (mid-min>kdbucket_size) D(x,y,z,r,min_photons,min,mid);\
			else for (int i=min;i<mid&&min_photons<=stars->map[i].photons;i++)kdcheck(i,x,y,z,r,min_photons);\
		}\
		if (mid<max) kdcheck(mid,x,y,z,r,min_photons);\
		if (kdresults_size==kdresults_maxsize) min_photons=stars->map[kdresults[kdresults_size-1]].photons;\
		if (mid+1<max &&stars->map[mid].B <= C) {\
			if (max-(mid+1)>kdbucket_size) D(x,y,z,r,min_photons,mid+1,max);\
			else for (int i=mid+1;i<max&&min_photons<=stars->map[i].photons;i++)kdcheck(i,x,y,z,r,min_photons);\
		}
	void kdsearch_x(const float x, const float y, const float z, const float r, float min_photons, int min, int max) {KDSEARCH_NEXT(x-r,x,x+r,kdsearch_y)}
	void kdsearch_y(const float x, const float y, const float z, const float r, float min_photons, int min, int max) {KDSEARCH_NEXT(y-r,y,y+r,kdsearch_z)}
	void kdsearch_z(const float x, const float y, const float z, const float r, float min_photons, int min, int max) {KDSEARCH_NEXT(z-r,z,z+r,kdsearch_x)}
	#undef KDSEARCH_NEXT

	void kdmask_filter_catalog() {
		for (int i=0;i<stars->map_size;i++) {
			//set mask for stars that are too bright,too close, highly variable, or unreliable
			int8_t lastmask=kdmask[i];
			//TODO: uncomment these two lines once secondary star matching is working
			
			/* kdsearch(map[i].x,map[i].y,map[i].z,3.5,BRIGHT_THRESH);
			 * if (kdresults_size>1||lastmask || map[i].unreliable>0)*/
			if (stars->map[i].photons<BRIGHT_THRESH) {
				kdmask[i]=1;
				reset_kdresults();
				continue;
			}
			undo_kdsearch();
		}
	}

	/* Masks the dimmest stars in each area to produce a 
	 * map with uniform density */
	void kdmask_uniform_density(int min_stars_per_fov) {
		std::set<int> uniform_set;
		int kdresults_maxsize_old=kdresults_maxsize;//TODO: eliminate this variable once restructured?
		kdresults_maxsize=min_stars_per_fov;
		for (int i=0;i<stars->map_size;i++) if (kdmask[i]==0) {
			kdsearch(stars->map[i].x,stars->map[i].y,stars->map[i].z,TODO_MINFOV_D2,BRIGHT_THRESH);
			for (int j=0;j<kdresults_size;j++) uniform_set.insert(kdresults[j]);
			undo_kdsearch();
		}
		for (int i=0;i<stars->map_size;i++) kdmask[i]=1;
		std::set<int>::iterator it = uniform_set.begin();
		for (int i=0; i<uniform_set.size();i++,it++) kdmask[*it]=0;
		kdresults_maxsize=kdresults_maxsize_old;
	}
};
#endif
