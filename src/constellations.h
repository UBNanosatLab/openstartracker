#ifndef CONSTELLATIONS_H
#define CONSTELLATIONS_H

#include "stars.h"

struct constellation {
	float p;
	int32_t s1;
	int32_t s2;
	int32_t idx;
	void DBG_(const char *s) {
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
	void DBG_(const char *s) {
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
	
	star_db* stars;
	star_query* results;
		
	constellation_db(star_db *s,int stars_per_fov, int from_image) {
		stars=s;
		results=new star_query(s);

		if (from_image) {
			map=NULL;
			std::sort(stars->map, stars->map+stars->map_size,star_gt_flux);
			int ns=stars->map_size;/* number of stars to check */
			if (ns>stars_per_fov) ns=stars_per_fov;//MAX_FALSE_STARS+2
			
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
		} else {
			results->kdmask_uniform_density(stars_per_fov);//2+DB_REDUNDANCY
			std::set<constellation,constellation_lt> c_set;
			for (int i=0;i<stars->map_size;i++) if (results->kdmask[i]==0) {
				results->kdsearch(stars->map[i].x,stars->map[i].y,stars->map[i].z,MAXFOV,THRESH_FACTOR*IMAGE_VARIANCE);
				constellation c;
				for (int j=0;j<results->kdresults_size;j++) if (i!=results->kdresults[j] && stars->map[i].flux>=stars->map[results->kdresults[j]].flux){
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
	}
	~constellation_db() {
		free(map);
		delete results;
		delete stars;
	}
	void DBG_(const char *s) {
		DBG_PRINT("%s\n",s);
		stars->DBG_("STARS");
		results->DBG_("RESULTS");
		for (int i=0; i<map_size; i++) {
			DBG_PRINT("%d:\t",i);
			map[i].DBG_("C");
		}
	}
};

#endif
