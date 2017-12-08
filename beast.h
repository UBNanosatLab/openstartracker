#include <stdint.h>
#include <math.h>	   /* sqrt */
#include <float.h>
#include <limits.h>
#include <string.h>
#include <search.h>
#include <assert.h>

#include <algorithm>
#include <set>

#include "config.h"

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


struct constellation_db {
	constellation* map;
	star_db* stars;
	int map_size;
	constellation_db() {map=NULL;map_size=0;}
	~constellation_db() {free(map);}
	//TODO: save and load
	//void save (FILE *fp, int32_t node) {
	//	fwrite (&next_id, sizeof(int32_t), 1, fp);
	//	int result = fwrite (the_tree, sizeof(kdnode), next_id, fp);
	//	if (result != next_id) perror ("Error writing kdtree\n");
	//}
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

struct beast_db {
	star_db* stars;
	star_query* results;
	constellation_db* constellations;
	//If no star list is provided, load from catalog
	beast_db() {
		
		stars=new star_db;
		stars->load_catalog();
		results=new star_query(stars);
		
		results->kdmask_filter_catalog();
		results->kdmask_uniform_density(REQUIRED_STARS);
		std::set<constellation> c_set;
		for (int i=0;i<stars->map_size;i++) if (results->kdmask[i]==0) {
			results->kdsearch(stars->map[i].x,stars->map[i].y,stars->map[i].z,TODO_MAXFOV_D2*2,BRIGHT_THRESH);
			constellation c;
			for (int j=0;j<results->kdresults_size;j++) if (i!=results->kdresults[j] && stars->map[i].photons>=stars->map[results->kdresults[j]].photons){
				c.p=stars->map[i].dist_arcsec(stars->map[results->kdresults[j]]);
				c.s1=i;
				c.s2=results->kdresults[j];
				c_set.insert(c);
			}
			results->undo_kdsearch();
		}
		constellations = new constellation_db;
		constellations->stars=stars;
		constellations->map_size=c_set.size();
		constellations->map=(constellation*)malloc(constellations->map_size*sizeof(constellations->map[0]));
		std::set<constellation>::iterator it = c_set.begin();
		for (int idx=0; idx<constellations->map_size;idx++,it++) {
			constellations->map[idx]=*it;
			constellations->map[idx].idx=idx;
		}
		
		results->reset_kdmask();
		results->reset_kdresults();
		
	}
	
	//Otherwise, this is an image
	beast_db(star_db *s) {
		stars=s;
		results=NULL;
		
		std::sort(stars->map, stars->map+stars->map_size,star_gt_photons);
		int ns=stars->map_size;/* number of stars to check */
		if (ns>REQUIRED_STARS+MAX_FALSE_STARS) ns=REQUIRED_STARS+MAX_FALSE_STARS;
		
		constellations = new constellation_db;
		constellations->stars=stars;
		constellations->map_size=ns*(ns-1)/2;
		constellations->map=(constellation*)malloc(constellations->map_size*sizeof(constellations->map[0]));
		int idx=0;
		for (int j=1;j<ns;j++) for (int i=0;i<j;i++,idx++) {
			constellations->map[idx].p=stars->map[i].dist_arcsec(stars->map[j]);
			constellations->map[idx].s1=i;
			constellations->map[idx].s2=j;
			constellations->map[idx].idx=idx;
		}
		std::sort(constellations->map, constellations->map+constellations->map_size,constellation_lt_p);
	}
	~beast_db() {
		delete stars;
		delete results;
		delete constellations;
	}
};

std::set<char*> env_str;//
beast_db *DB;
void load_config() {
	
	//move calibration.txt to constellation_db
	//move dbstats.txt to beastdb
	/* load config */
	
	FILE *stream = fopen("calibration.txt", "r");
	if (stream == NULL) exit(EXIT_FAILURE);

	ssize_t read;
	char *line = NULL;
	size_t len = 0;
	while ((read = getline(&line, &len, stream)) != -1) {
		char * s=(char *)malloc(sizeof(char) * len);
		putenv(strcpy(s,line));
		env_str.insert(s);
	}
	free(line);
	fclose(stream);

	IMG_X=atoi(getenv("IMG_X"));
	IMG_Y=atoi(getenv("IMG_Y"));
	TODO_MAXFOV_D2=sqrt(IMG_X*IMG_X+IMG_Y*IMG_Y)/2;
	TODO_MINFOV_D2=IMG_Y/2;
	DEG_X=atof(getenv("DEG_X"));
	DEG_Y=atof(getenv("DEG_Y"));
	POS_ERR_SIGMA=atof(getenv("POS_ERR_SIGMA"));
	PIXSCALE=atof(getenv("PIXSCALE"));
	POS_VARIANCE=atof(getenv("POS_VARIANCE"));/* sigma_r^2 */
	IMAGE_VARIANCE=atof(getenv("IMAGE_VARIANCE"));/* lambda */
	BRIGHT_THRESH=atof(getenv("BRIGHT_THRESH"));
	EXPECTED_FALSE_STARS=atof(getenv("EXPECTED_FALSE_STARS"));/* pfalse */
	MAX_FALSE_STARS=atoi(getenv("MAX_FALSE_STARS"));/* >10 is slow */
	REQUIRED_STARS=atoi(getenv("REQUIRED_STARS"));/* >10 is slow */
	PHOTONS=atoi(getenv("PHOTONS"));
	MATCH_VALUE=4*log(EXPECTED_FALSE_STARS/(IMG_X*IMG_Y))+log(2*PI);/* base */

	PIXX_TANGENT=2*tan(DEG_X*PI/(180*2))/IMG_X;
	PIXY_TANGENT=2*tan(DEG_Y*PI/(180*2))/IMG_Y;
	
	DB=new beast_db;
}

struct constellation_match {
	float R11,R12,R13;
	float R21,R22,R23;
	float R31,R32,R33;
	
	struct constellation_score *c_scores;
	int c_scores_size;
	float *winner_scores;
	//winner_id_map[new]=old
	int32_t *winner_id_map;
	float p_match;
	/* weighted_triad results */

	/* see https://en.wikipedia.org/wiki/Triad_method */
	/* and http://nghiaho.com/?page_id=846 */
	/* returns match results */

	/* when compiled, this section contains roughly 430 floating point operations */
	/* according to https://www.karlrupp.net/2016/02/gemm-and-stream-results-on-intel-edison/ */
	/* we can perform >250 MFLOPS with doubles, and >500 MFLOPS with floats */


	void weighted_triad(star db_s1,star db_s2,star img_s1,star img_s2){
		/* v=A*w */
		float wa1=db_s1.x,wa2=db_s1.y,wa3=db_s1.z;
		float wb1=db_s2.x,wb2=db_s2.y,wb3=db_s2.z;
		float va1=img_s1.x,va2=img_s1.y,va3=img_s1.z;
		float vb1=img_s2.x,vb2=img_s2.y,vb3=img_s2.z;
		float wc1=wa2*wb3 - wa3*wb2;
		float wc2=wa3*wb1 - wa1*wb3;
		float wc3=wa1*wb2 - wa2*wb1;
		float wcnorm=sqrt(wc1*wc1+wc2*wc2+wc3*wc3);
		wc1/=wcnorm;
		wc2/=wcnorm;
		wc3/=wcnorm;

		float vc1=va2*vb3 - va3*vb2;
		float vc2=va3*vb1 - va1*vb3;
		float vc3=va1*vb2 - va2*vb1;
		float vcnorm=sqrt(vc1*vc1+vc2*vc2+vc3*vc3);
		vc1/=vcnorm;
		vc2/=vcnorm;
		vc3/=vcnorm;
		
		float vaXvc1=va2*vc3 - va3*vc2;
		float vaXvc2=va3*vc1 - va1*vc3;
		float vaXvc3=va1*vc2 - va2*vc1;

		float waXwc1=wa2*wc3 - wa3*wc2;
		float waXwc2=wa3*wc1 - wa1*wc3;
		float waXwc3=wa1*wc2 - wa2*wc1;
		
		/* some of these are unused */
		float A11=va1*wa1 + vaXvc1*waXwc1 + vc1*wc1;
		/* float A12=va1*wa2 + vaXvc1*waXwc2 + vc1*wc2; */
		/* float A13=va1*wa3 + vaXvc1*waXwc3 + vc1*wc3; */
		float A21=va2*wa1 + vaXvc2*waXwc1 + vc2*wc1;
		/* float A22=va2*wa2 + vaXvc2*waXwc2 + vc2*wc2; */
		/* float A23=va2*wa3 + vaXvc2*waXwc3 + vc2*wc3; */
		float A31=va3*wa1 + vaXvc3*waXwc1 + vc3*wc1;
		float A32=va3*wa2 + vaXvc3*waXwc2 + vc3*wc2;
		float A33=va3*wa3 + vaXvc3*waXwc3 + vc3*wc3;
		
		wc1=-wc1;
		wc2=-wc2;
		wc3=-wc3;
		
		vc1=-vc1;
		vc2=-vc2;
		vc3=-vc3;
		float vbXvc1=vb2*vc3 - vb3*vc2;
		float vbXvc2=vb3*vc1 - vb1*vc3;
		float vbXvc3=vb1*vc2 - vb2*vc1;
		
		float wbXwc1=wb2*wc3 - wb3*wc2;
		float wbXwc2=wb3*wc1 - wb1*wc3;
		float wbXwc3=wb1*wc2 - wb2*wc1;

		/* some of these are unused */
		float B11=vb1*wb1 + vbXvc1*wbXwc1 + vc1*wc1;
		/* float B12=vb1*wb2 + vbXvc1*wbXwc2 + vc1*wc2; */
		/* float B13=vb1*wb3 + vbXvc1*wbXwc3 + vc1*wc3; */
		float B21=vb2*wb1 + vbXvc2*wbXwc1 + vc2*wc1;
		/* float B22=vb2*wb2 + vbXvc2*wbXwc2 + vc2*wc2; */
		/* float B23=vb2*wb3 + vbXvc2*wbXwc3 + vc2*wc3; */
		float B31=vb3*wb1 + vbXvc3*wbXwc1 + vc3*wc1;
		float B32=vb3*wb2 + vbXvc3*wbXwc2 + vc3*wc2;
		float B33=vb3*wb3 + vbXvc3*wbXwc3 + vc3*wc3;
		
		/* use weights based on magnitude */
		/* weighted triad */
		float weightA=1.0/(db_s1.sigma_sq+img_s1.sigma_sq);
		float weightB=1.0/(db_s2.sigma_sq+img_s2.sigma_sq);

		float sumAB=weightA+weightB;
		weightA/=sumAB;
		weightB/=sumAB;
		
		float cz,sz,mz;
		float cy,sy,my;
		float cx,sx,mx;
		
		cz=weightA*A11+weightB*B11;
		sz=weightA*A21+weightB*B21;
		mz=sqrt(cz*cz+sz*sz);
		cz=cz/mz;
		sz=sz/mz;
		
		cy=weightA*sqrt(A32*A32+A33*A33)+weightB*sqrt(B32*B32+B33*B33);
		sy=-weightA*A31-weightB*B31;
		my=sqrt(cy*cy+sy*sy);
		cy=cy/my;
		sy=sy/my;
		
		cx=weightA*A33+weightB*B33;
		sx=weightA*A32+weightB*B32;
		mx=sqrt(cx*cx+sx*sx);
		cx=cx/mx;
		sx=sx/mx;
		
		R11=cy*cz;
		R12=cz*sx*sy - cx*sz;
		R13=sx*sz + cx*cz*sy;
		
		R21=cy*sz;
		R22=cx*cz + sx*sy*sz;
		R23=cx*sy*sz - cz*sx;
		
		R31=-sy;
		R32=cy*sx;
		R33=cx*cy;
	}
	void add_score(constellation *db_const, constellation_score *cs, int32_t *img_mask, beast_db *db, beast_db *img){
		cs->db_id1=db_const->s1;
		cs->db_id2=db_const->s2;
		cs->id_map=(int32_t *)malloc(sizeof(int32_t)*img->stars->map_size);
		cs->scores=(float *)malloc(sizeof(float)*img->stars->map_size);
		for (int i=0;i<img->stars->map_size;i++) cs->id_map[i]=-1;
		for (int i=0;i<img->stars->map_size;i++) cs->scores[i]=0.0;
		
		
		cs->totalscore=log(EXPECTED_FALSE_STARS/(IMG_X*IMG_Y))*(2*img->stars->map_size);
		db->results->kdsearch(R11,R12,R13,TODO_MAXFOV_D2,BRIGHT_THRESH);
		for(int32_t i=0;i<db->results->kdresults_size;i++){
			int32_t o=db->results->kdresults[i];
			star s=db->stars->map[o];
			float x=s.x*R11+s.y*R12+s.z*R13;
			float y=s.x*R21+s.y*R22+s.z*R23;
			float z=s.x*R31+s.y*R32+s.z*R33;
			float px=y/(x*PIXX_TANGENT);
			float py=z/(x*PIXY_TANGENT);
			int nx,ny;
			nx=(int)(px+IMG_X/2.0f);
			ny=(int)(py+IMG_Y/2.0f);
			int32_t n=-1;
			if (nx==-1) nx++;
			else if (nx==IMG_X) nx--;
			if (ny==-1) ny++;
			else if (ny==IMG_Y) ny--;
			if (nx>=0&&nx<IMG_X&&ny>=0&&ny<IMG_Y) n=img_mask[nx+ny*IMG_X];
			if (n!=-1) {
				float sigma_sq=img->stars->map[n].sigma_sq+db->stars->max_variance;
				float maxdist_sq=-sigma_sq*(log(sigma_sq)+MATCH_VALUE);
				float a=(px-img->stars->map[n].px);
				float b=(py-img->stars->map[n].py);
				float score = (maxdist_sq-(a*a+b*b))/(2*sigma_sq);
				/* only match the closest star */
				if (score>cs->scores[n]){
					cs->id_map[n]=o;
					cs->scores[n]=score;
				}
			}
		}
		db->results->undo_kdsearch();
		for(int n=0;n<img->stars->map_size;n++) {
			cs->totalscore+=cs->scores[n];
		}
		
	}
	constellation_match(beast_db *db, beast_db *img) {
		c_scores=NULL;
		c_scores_size=0;
		/* Do we have enough stars? */
		if (db->stars->map_size<2||img->stars->map_size<2) return;

		winner_id_map=(int *)malloc(sizeof(int)*img->stars->map_size);
		winner_scores=(float *)malloc(sizeof(float)*img->stars->map_size);
		
		int32_t *img_mask = img->stars->get_img_mask(db->stars->max_variance);
		
		for (int n=0;n<img->constellations->map_size;n++) {
			constellation lb=img->constellations->map[n];
			constellation ub=img->constellations->map[n];
			lb.p-=POS_ERR_SIGMA*PIXSCALE*sqrt(img->stars->map[lb.s1].sigma_sq+img->stars->map[lb.s2].sigma_sq+2*db->stars->max_variance);
			ub.p+=POS_ERR_SIGMA*PIXSCALE*sqrt(img->stars->map[ub.s1].sigma_sq+img->stars->map[ub.s2].sigma_sq+2*db->stars->max_variance);
			constellation *lower=std::lower_bound (db->constellations->map, db->constellations->map+db->constellations->map_size, lb,constellation_lt_p);	
			constellation *upper=std::upper_bound (db->constellations->map, db->constellations->map+db->constellations->map_size, ub,constellation_lt_p);
			//rewind by one
			upper--;
			
			//TODO: get rid of cscores 
			if (lower->idx<=upper->idx) c_scores=(struct constellation_score*)realloc(c_scores,sizeof(struct constellation_score)*(c_scores_size+(upper->idx-lower->idx+1)*2));
			for (int o=lower->idx;o<=upper->idx;o++) {
				int32_t db_idx1=db->constellations->map[o].s1;
				int32_t db_idx2=db->constellations->map[o].s2;
				int32_t img_idx1=img->constellations->map[n].s1;
				int32_t img_idx2=img->constellations->map[n].s2;
				
				star db_s1=db->stars->map[db_idx1];
				star db_s2=db->stars->map[db_idx2];
				star img_s1=img->stars->map[img_idx1];
				star img_s2=img->stars->map[img_idx2];
				
				/* try both orderings of stars */
				weighted_triad(db_s1,db_s2,img_s1,img_s2);
				c_scores[c_scores_size].img_id1=img_idx1;
				c_scores[c_scores_size].img_id2=img_idx2;
				add_score(&(db->constellations->map[o]),&c_scores[c_scores_size],img_mask,db,img);
				c_scores_size++;
				
				weighted_triad(db_s1,db_s2,img_s2,img_s1);
				c_scores[c_scores_size].img_id1=img_idx2;
				c_scores[c_scores_size].img_id2=img_idx1;
				add_score(&(db->constellations->map[o]),&c_scores[c_scores_size],img_mask,db,img);
				c_scores_size++;
			}
		}
		std::sort(c_scores, c_scores+c_scores_size);
		for (int i=0;i<img->stars->map_size;i++) {winner_id_map[i]=-1;winner_scores[i]=0.0f;}
		if (c_scores_size>0) {
			for(int n=0;n<img->stars->map_size;n++) {
				int o=c_scores[0].id_map[n];
				if (o!=-1){
					winner_scores[img->stars->map[n].star_idx]=c_scores[0].scores[n];
					winner_id_map[img->stars->map[n].star_idx]=db->stars->map[o].id;
				}
			}
			
		}
		//TODO: move to add_score
		/* add up probabilities of all matches, excluding those which */
		/* are equivalent to the best match (S1==s1,S2=s2) */
		p_match=1.0;
		if (c_scores_size>0) {
			float bestscore=c_scores[0].totalscore;
			int db_id1=c_scores[0].db_id1;
			int db_id2=c_scores[0].db_id2;
			int img_id1=c_scores[0].img_id1;
			int img_id2=c_scores[0].img_id2;
			/* set attitude matrix to best match */
			weighted_triad(db->stars->map[db_id1],db->stars->map[db_id2],img->stars->map[img_id1],img->stars->map[img_id2]);
			for(int i=1;i<c_scores_size;i++) {
				if (c_scores[i].id_map[img_id1]!=db_id1&&c_scores[i].id_map[img_id2]!=db_id2){
					p_match+=exp(c_scores[i].totalscore-bestscore);
				}
			}
			//Turns out baysian hypothesis testing was the best way
			//after all. Who would've guessed?
			p_match=1.0/p_match;
		} else {
			p_match=0.0;
		}
		free(img_mask);
	}
	~constellation_match() {
		for (int i=0;i<c_scores_size; i++) {
			free(c_scores[i].scores);
			free(c_scores[i].id_map);
		}
		free(c_scores);
		free(winner_id_map);
		free(winner_scores);
	}
};
