#ifndef CONSTELLATIONS_H
#define CONSTELLATIONS_H

#include "stars.h"

struct constellation {
	float p;
	int32_t s1;
	int32_t s2;
	int32_t idx;
	void _print(const char *s) {
		DBG_PRINT("%s\t",s);
		DBG_PRINT("p=%f ",p);
		DBG_PRINT("s1=%d ",s1);
		DBG_PRINT("s2=%d ",s2);
		DBG_PRINT("idx=%d\n",idx);
	}
};

struct  constellation_pair {
	float totalscore;
	int32_t db_s1,db_s2;
	int32_t img_s1,img_s2;
	
	void flip() {
		int32_t t=img_s1;
		img_s1=img_s2;
		img_s2=t;
	}
	void _print(const char *s) {
		DBG_PRINT("%s\t",s);
		DBG_PRINT("totalscore=%f ",totalscore);
		DBG_PRINT("db_s1=%d ",db_s1);
		DBG_PRINT("db_s2=%d ",db_s2);
		DBG_PRINT("img_s1=%d ",img_s1);
		DBG_PRINT("img_s2=%d\n",img_s2);
	}
};

struct constellation_lt {
	bool operator() (const constellation &c1, const constellation &c2) {
		if (c1.p!=c2.p) return c1.p<c2.p;
		else if (c1.s1!=c2.s1) return c1.s1<c2.s1;
		else return c1.s2<c2.s2;
	}
};

bool constellation_lt_s1(const constellation &c1, const constellation &c2) {return c1.s1 < c2.s1;}
bool constellation_lt_s2(const constellation &c1, const constellation &c2) {return c1.s2 < c2.s2;}
bool constellation_lt_p(const constellation &c1, const constellation &c2) {return c1.p < c2.p;}

struct constellation_db {
	int map_size;
	constellation* map;
	
	star_db* orig_stars;
	star_query* orig_results;
	star_db* stars;
	star_query* results;
		
	//If no star list is provided, load from catalog
	constellation_db() {

		orig_stars=new star_db;
		orig_stars->load_catalog();
		
		orig_results=new star_query(orig_stars);
		orig_results->kdmask_filter_catalog();
		orig_results->kdmask_uniform_density(REQUIRED_STARS);
		stars=orig_results->from_kdmask();
		results=new star_query(stars);
		orig_results->reset_kdmask();
		
		results->kdmask_uniform_density(2+DB_REDUNDANCY);
		int nmask=0;
		for (int i=0;i<stars->map_size;i++) if (results->kdmask[i]==0) nmask++;
		std::set<constellation,constellation_lt> c_set;
		for (int i=0;i<stars->map_size;i++) if (results->kdmask[i]==0) {
			results->kdsearch(stars->map[i].x,stars->map[i].y,stars->map[i].z,MAXFOV_D2*2,BRIGHT_THRESH);
			constellation c;
			for (int j=0;j<results->kdresults_size;j++) if (i!=results->kdresults[j] && stars->map[i].photons>=stars->map[results->kdresults[j]].photons){
				c.p=stars->map[i].dist_arcsec(stars->map[results->kdresults[j]]);
				c.s1=i;
				c.s2=results->kdresults[j];
				c_set.insert(c);
			}
			results->clear_kdresults();
		}
		results->reset_kdmask();
		//preallocate map
		map_size=c_set.size();
		map=(constellation*)malloc(map_size*sizeof(map[0]));
		std::set<constellation>::iterator it = c_set.begin();
		for (int idx=0; idx<map_size;idx++,it++) {
			map[idx]=*it;
			map[idx].idx=idx;
		}
		
	}
	constellation_db(star_db *s) {constellation_img_init(s,MAX_FALSE_STARS+2);}
	constellation_db(star_db *s, int n_brightest) {constellation_img_init(s,n_brightest);}
	//Otherwise, this is an image
	void constellation_img_init(star_db *s,int n_brightest) {
		stars=s;
		results=new star_query(s);
		orig_stars=NULL;
		orig_results=NULL;
		
		std::sort(stars->map, stars->map+stars->map_size,star_gt_photons);
		int ns=stars->map_size;/* number of stars to check */
		if (ns>n_brightest) ns=n_brightest;
		
		stars=stars;
		map_size=ns*(ns-1)/2;
		map=(constellation*)malloc(map_size*sizeof(map[0]));
		
		int idx=0;
		for (int j=1;j<ns;j++) for (int i=0;i<j;i++,idx++) {
			map[idx].p=stars->map[i].dist_arcsec(stars->map[j]);
			map[idx].s1=i;
			map[idx].s2=j;
		}
		std::sort(map, map+map_size,constellation_lt_p);
		while (--idx>=0) map[idx].idx=idx;
	}
	~constellation_db() {
		free(map);
		delete results;
		delete stars;
		delete orig_results;
		delete orig_stars;
	}
	void _print(const char *s) {
		DBG_PRINT("%s\n",s);
		results->_print("RESULTS");
		if (orig_results!=NULL) {
			orig_results->_print("ORIG_RESULTS");
		}
		map[0]._print("CONSTELLATION");
		DBG_PRINT(".\n%d total\n.\n",map_size);
		map[map_size-1]._print("CONSTELLATION");
	}

};

#endif
