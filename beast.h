#include <stdint.h>
#include <math.h>	   /* sqrt */
#include <float.h>
#include <limits.h>
#include <string.h>
#include <search.h>
#include <assert.h>

#include <algorithm>
#include <set>

namespace beast {
	#define PI		   3.14159265358979323846  /* pi */
	#define TWOPI		6.28318530717958647693
	struct constellation {
		float p;
		int32_t s1;
		int32_t s2;
		int32_t idx;
		bool operator<(const constellation &c) const {
			if (p!=c.p) return p<c.p;
			else if (s1!=c.s1) return s1<c.s1;
			else return s2<c.s2;
		}
	};
	bool constellation_lt_s1(const constellation &c1, const constellation &c2) {return c1.s1 < c2.s1;}
	bool constellation_lt_s2(const constellation &c1, const constellation &c2) {return c1.s2 < c2.s2;}
	bool constellation_lt_p(const constellation &c1, const constellation &c2) {return c1.p < c2.p;}

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

	struct  constellation_score {
		float totalscore;
		int32_t db_id1,db_id2;
		int32_t img_id1,img_id2;
		int *id_map; /* Usage: id_map[newstar]=oldstar */
		float *scores;
		/* backwards to sort in decending order */
		bool operator< (const constellation_score* c) { return totalscore > c->totalscore; }
		bool operator< (const constellation_score c) { return totalscore > c.totalscore; }
	};
	int IMG_X,IMG_Y,MAX_FALSE_STARS,REQUIRED_STARS;
	float DEG_X,DEG_Y,PIXX_TANGENT,PIXY_TANGENT;
	float PIXSCALE,POS_ERR_SIGMA,POS_VARIANCE;
	float IMAGE_VARIANCE,EXPECTED_FALSE_STARS,MATCH_VALUE;
	float PHOTONS,BRIGHT_THRESH;
	float TODO_MAXFOV_D2;
	float TODO_MINFOV_D2;

	void load_config() {
		
		//TODO: use an array of strings for environment variables
		//move calibration.txt to imgdb
		//move dbstats.txt to beastdb
		/* load config */
		
		FILE *stream = fopen("calibration/calibration.txt", "r");
		if (stream == NULL) exit(EXIT_FAILURE);

		ssize_t read;
		char *line = NULL;
		size_t len = 0;
		while ((read = getline(&line, &len, stream)) != -1) {
			/* This is not a true memory leak, as these variables are
			 * accessable via getenv(). That being said, it might be useful
			 * to add a dummy structure to keep track of them when debugging 
			 * with valgrind */
			putenv(strcpy((char *)malloc(sizeof(char) * len),line));
		}
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
		/* TODO: need to use this */
		MAX_FALSE_STARS=atoi(getenv("MAX_FALSE_STARS"));/* >10 is slow */
		REQUIRED_STARS=atoi(getenv("REQUIRED_STARS"));/* >10 is slow */
		PHOTONS=atoi(getenv("PHOTONS"));
		MATCH_VALUE=4*log(EXPECTED_FALSE_STARS/(IMG_X*IMG_Y))+log(2*PI);/* base */

		PIXX_TANGENT=2*tan(DEG_X*PI/(180*2))/IMG_X;
		PIXY_TANGENT=2*tan(DEG_Y*PI/(180*2))/IMG_Y;
	}

	struct stardb {
		constellation* db_map;
		int db_map_size;
		star *star_map;
		int star_map_size;
		
		int *kdresults;
		int8_t *kdmask;
		int kdresults_size;
		int kdresults_maxsize;

		float max_variance;

		stardb() {
			db_map=NULL;
			star_map=NULL;
			db_map_size=0;
			star_map_size=0;

			kdmask=NULL;
			kdresults=NULL;
			kdresults_size=0;
			kdresults_maxsize=INT_MAX;
			
			max_variance=0.0;
		}
		~stardb() {
			free(db_map);
			free(star_map);
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
		 */	
			
		#define KDSORT_NEXT(A,B)\
			int mid=(min+max)/2;\
			if (min+1<max) {\
				std::nth_element(star_map+min,star_map+mid,star_map+max,A);\
				B(min,mid);B(mid+1,max);\
			}
		void kdsort() {kdsort_x(0,star_map_size);}
		void kdsort_x(int min, int max) {KDSORT_NEXT(star_lt_x,kdsort_y)}
		void kdsort_y(int min, int max) {KDSORT_NEXT(star_lt_y,kdsort_z)}
		void kdsort_z(int min, int max) {KDSORT_NEXT(star_lt_z,kdsort_photons)}
		void kdsort_photons(int min, int max) {KDSORT_NEXT(star_gt_photons,kdsort_x)}
		#undef KDSORT_NEXT
		 
		/* Reset kdresults without clearing kdmask
		 * Can be used to do clever things like reusing a mask generated
		 * by a previous search */
		void reset_kdresults() {kdresults_size=0;}
		
		/* Clears kdmask, but does not reset kdresults. Slow.*/
		void reset_kdmask() {memset(kdmask,0,sizeof(kdmask[0])*star_map_size);}
		
		/* Undoes the effects of all kdsearches since the last time
		 * reset_kdresults() was called. Much faster than reset_kdmask()*/
		 
		void undo_kdsearch() {
			while (kdresults_size>0) {
				kdresults_size--;
				kdmask[kdresults[kdresults_size]]--;
			}
		}
		 
		void kdcheck(int idx, float x, float y, float z, float r, float min_photons){
			x-=star_map[idx].x;
			y-=star_map[idx].y;
			z-=star_map[idx].z;
			if (x-r <= 0 && 0 <= x+r)
			if (y-r <= 0 && 0 <= y+r)
			if (z-r <= 0 && 0 <= z+r)
			if (min_photons <= star_map[idx].photons)
			if (kdmask[idx] <= 0)
			if (x*x+y*y+z*z<=r*r) {
				kdmask[idx]=1;
				/* Insertion sort into list from brightest to dimmest.
				 * Results come out of the kdtree partially sorted, 
				 * so this is guaranteed to be fast 
				 * 
				 * Note: Different sorting algorithms can result in different
				 * stars being selected in the case where there are two candidates of
				 * equal brightness for the last star. Do not be alarmed -
				 * they all work fine
				 */
				 
				int n=kdresults_size++;
				for (;n>0&&star_map[idx].photons>star_map[kdresults[n-1]].photons;n--) kdresults[n]=kdresults[n-1];
				kdresults[n]=idx;
				//if we go over the maximum, bump the dimmest star from the results
				if (kdresults_size>kdresults_maxsize) {
					kdresults_size=kdresults_maxsize;
					kdmask[kdresults[kdresults_size]]--;
				}
			}
		}
		#define KDSEARCH_NEXT(A,B,C,D)\
			int mid=(min+max)/2;\
			if (min<mid &&A <= star_map[mid].B) D(min,mid,x,y,z,r,min_photons);\
			if (mid<max) kdcheck(mid,x,y,z,r,min_photons);\
			if (kdresults_size==kdresults_maxsize) min_photons=star_map[kdresults[kdresults_size-1]].photons;\
			if (mid+1<max &&star_map[mid].B <= C) D(mid+1,max,x,y,z,r,min_photons);
		/**
		 * @fn		kdsearch(float x, float y, float z, float r, float min_photons)
		 * @param	x	cos(deg2rad*RA)*cos(deg2rad*DEC)
		 * @param	y	sin(deg2rad*RA)*cos(deg2rad*DEC)
		 * @param	z	sin(deg2rad*DEC)
		 * @param	r	search distance (pixels)
		 * @brief	search star_map for points within r pixels of x,y,z 
		 * @details put all results found into kdresults (sorted by brightness), mask via kdmask
		 */
		void kdsearch(float x, float y, float z, float r, float min_photons) {
			float r_deg=r*PIXSCALE/3600.0;
			float r_rad=r_deg*PI/180.0;
			kdsearch_x(0,star_map_size,x, y, z, 2*fabs(sin(r_rad/2.0)), min_photons);
		}
		void kdsearch_x(int min, int max, const float x, const float y, const float z, const float r, float min_photons) {KDSEARCH_NEXT(x-r,x,x+r,kdsearch_y)}
		void kdsearch_y(int min, int max, const float x, const float y, const float z, const float r, float min_photons) {KDSEARCH_NEXT(y-r,y,y+r,kdsearch_z)}
		void kdsearch_z(int min, int max, const float x, const float y, const float z, const float r, float min_photons) {KDSEARCH_NEXT(z-r,z,z+r,kdsearch_photons)}
		void kdsearch_photons(int min, int max, const float x, const float y, const float z, const float r, float min_photons) {KDSEARCH_NEXT(FLT_MAX*-1,photons*-1,min_photons*-1,kdsearch_x)}
		#undef KDSEARCH_NEXT
	};


	struct imgdb: public stardb {
		//TODO: is it actually necessary to have two versions of add_star()?
		void add_star(float x, float y, float z, float mag, int id) {
			int n=star_map_size++;
			if (n%16==0) star_map=(star*)realloc(star_map,(star_map_size+16)*sizeof(star_map[0]));

			star_map[n].x=x;
			star_map[n].y=y;
			star_map[n].z=z;
			star_map[n].id=id;
			star_map[n].unreliable=0;
			star_map[n].photons=PHOTONS*powf(10.0,-mag/2.5);/* TODO: change this to pixel value */
			star_map[n].px=y/(x*PIXX_TANGENT);
			star_map[n].py=z/(x*PIXY_TANGENT);
			star_map[n].sigma_sq=POS_VARIANCE;
			star_map[n].star_idx=n;
		}
		void add_star(float px, float py, float mag) {
			float j=PIXX_TANGENT*px; /* j=(y/x) */
			float k=PIXY_TANGENT*py; /* k=z/x */
			float x=1./sqrt(j*j+k*k+1);
			float y=j*x;
			float z=k*x;
			int n=star_map_size;
			add_star(x,y,z,mag,-1);
			star_map[n].sigma_sq=IMAGE_VARIANCE/star_map[n].photons;
			if (max_variance<star_map[n].sigma_sq) max_variance=star_map[n].sigma_sq;
		}
		
		void load_catalog() {
			FILE *stream = fopen("catalog.dat", "r");
			if (stream == NULL) exit(EXIT_FAILURE);
			star_map_size=0;
			while(!feof(stream)) if(fgetc(stream) == '\n') star_map_size++;
			rewind(stream);
			max_variance=POS_VARIANCE;
			
			star_map = (star*)realloc(star_map,star_map_size*sizeof(star_map[0]));
			ssize_t read;
			char *line = NULL;
			size_t len = 0;
			
			for(int i=0;i<star_map_size;i++){
				if ((read = getline(&line, &len, stream)) != -1) {
					star_map[i].id=atoi(strtok(line," "));
					star_map[i].star_idx=i;
					float mag=atof(strtok(NULL," "));
					star_map[i].photons=PHOTONS*powf(10.0,-mag/2.5);
					
					star_map[i].x=atof(strtok(NULL," "));
					star_map[i].y=atof(strtok(NULL," "));
					star_map[i].z=atof(strtok(NULL," "));
					star_map[i].unreliable=atoi(strtok(NULL," "));
					
					star_map[i].sigma_sq=POS_VARIANCE;
				}
			}
			fclose(stream);
		}
		void kdmask_filter_catalog() {
			for (int i=0;i<star_map_size;i++) {
				//set mask for stars that are too bright,too close, highly variable, or unreliable
				int8_t lastmask=kdmask[i];
				//TODO: uncomment these two lines once secondary star matching is working
				
				/* kdsearch(star_map[i].x,star_map[i].y,star_map[i].z,3.5,BRIGHT_THRESH);
				 * if (kdresults_size>1||lastmask || star_map[i].unreliable>0)*/
				if (star_map[i].photons<BRIGHT_THRESH) {
					kdmask[i]=1;
					reset_kdresults();
					continue;
				}
				undo_kdsearch();
			}
		}
		//TODO: split functions which require kdsorted stars into a seprate  class
		//or add flag + assert() to kdsearch
		
		/* Masks the dimmest stars in each area to produce a 
		 * star_map with uniform density */
		void kdmask_uniform_density(int min_stars_per_fov) {
			std::set<int> uniform_set;
			int kdresults_maxsize_old=kdresults_maxsize;//TODO: eliminate this variable once restructured?
			kdresults_maxsize=min_stars_per_fov;
			for (int i=0;i<star_map_size;i++) if (kdmask[i]==0) {
				kdsearch(star_map[i].x,star_map[i].y,star_map[i].z,TODO_MINFOV_D2,BRIGHT_THRESH);
				for (int j=0;j<kdresults_size;j++) uniform_set.insert(kdresults[j]);
				undo_kdsearch();
			}
			for (int i=0;i<star_map_size;i++) kdmask[i]=1;
			std::set<int>::iterator it = uniform_set.begin();
			for (int i=0; i<uniform_set.size();i++,it++) kdmask[*it]=0;
			kdresults_maxsize=kdresults_maxsize_old;
		}
		
		void gendb_catalog() {
			kdsort();
			kdresults=(int*)realloc(kdresults,star_map_size*sizeof(kdresults[0]));
			kdmask=(int8_t*)realloc(kdmask,star_map_size*sizeof(kdmask[0]));
			reset_kdresults();
			reset_kdmask();

			kdmask_filter_catalog();
			int ct=0;
			kdmask_uniform_density(REQUIRED_STARS);
			for (int i=0; i<star_map_size;i++) if (kdmask[i]==0) ct++;
			std::set<constellation> c_set;
			for (int i=0;i<star_map_size;i++) if (kdmask[i]==0) {
				kdsearch(star_map[i].x,star_map[i].y,star_map[i].z,TODO_MAXFOV_D2*2,BRIGHT_THRESH);
				constellation c;
				for (int j=0;j<kdresults_size;j++) if (i!=kdresults[j] && star_map[i].photons>=star_map[kdresults[j]].photons){
					c.p=star_map[i].dist_arcsec(star_map[kdresults[j]]);
					c.s1=i;
					c.s2=kdresults[j];
					c_set.insert(c);
				}
				undo_kdsearch();
			}
			db_map_size=c_set.size();
			db_map=(constellation*)realloc(db_map,db_map_size*sizeof(db_map[0]));
			std::set<constellation>::iterator it = c_set.begin();
			for (int idx=0; idx<db_map_size;idx++,it++) {
				db_map[idx]=*it;
				db_map[idx].idx=idx;
			}
		}
		
		void gendb_img() {
			std::sort(star_map, star_map+star_map_size,star_gt_photons);
			int ns=star_map_size;/* number of stars to check */
			if (ns>REQUIRED_STARS+MAX_FALSE_STARS) ns=REQUIRED_STARS+MAX_FALSE_STARS;
			db_map_size=ns*(ns-1)/2;
			db_map=(constellation*)realloc(db_map,db_map_size*sizeof(db_map[0]));
			int idx=0;
			for (int j=1;j<ns;j++) for (int i=0;i<j;i++,idx++) {
				db_map[idx].p=star_map[i].dist_arcsec(star_map[j]);
				db_map[idx].s1=i;
				db_map[idx].s2=j;
				db_map[idx].idx=idx;
			}

			std::sort(db_map, db_map+db_map_size,constellation_lt_p);

		}
		//TODO: save and load
		//void save (FILE *fp, int32_t node) {
		//	fwrite (&next_id, sizeof(int32_t), 1, fp);
		//	int result = fwrite (the_tree, sizeof(kdnode), next_id, fp);
		//	if (result != next_id) perror ("Error writing kdtree\n");
		//}
		
		int* get_mask(float db_max_variance);
	};

	imgdb* DB;
	void load_db() {
		load_config();
		DB=new imgdb;
		DB->load_catalog();
		DB->gendb_catalog();
		//TODO: comment this out - instead use a seperate list once match has been confirmed
		DB->reset_kdmask();
		DB->reset_kdresults();
	}


	void add_score(constellation *db_const, constellation_score *cs, int32_t *img_mask, stardb *db, stardb *img,float db_max_variance);
	void find(stardb *db, imgdb *img);


	/* weighted_triad results */
	/* TODO: do this with something other than global variables */

	float R11,R12,R13;
	float R21,R22,R23;
	float R31,R32,R33;

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
}
