#ifndef CONSTELLATIONS_H
#define CONSTELLATIONS_H

#include "stars.h"

struct constellation {
	float p;
	int32_t s1;
	int32_t s2;
	int32_t idx;
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
		results=orig_results->reduce_stars();
		orig_results->reset_kdmask();
		
		results->kdmask_uniform_density(2+DB_REDUNDANCY);
		stars=results->stars;
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
			results->undo_kdsearch();
		}
		//TODO: comment these out, and add second starid stage
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
	
	//Otherwise, this is an image
	constellation_db(star_db *s) {
		stars=s;
		orig_stars=s;
		//TODO: change this to results = sq
		results=NULL;
		orig_results=NULL;
		
		std::sort(stars->map, stars->map+stars->map_size,star_gt_photons);
		int ns=stars->map_size;/* number of stars to check */
		//if (ns>REQUIRED_STARS+MAX_FALSE_STARS) ns=REQUIRED_STARS+MAX_FALSE_STARS;
		if (ns>MAX_FALSE_STARS+2) ns=MAX_FALSE_STARS+2;
		
		stars=stars;
		map_size=ns*(ns-1)/2;
		map=(constellation*)malloc(map_size*sizeof(map[0]));
		
		int idx=0;
		for (int j=1;j<ns;j++) for (int i=0;i<j;i++,idx++) {
			map[idx].p=stars->map[i].dist_arcsec(stars->map[j]);
			map[idx].s1=i;
			map[idx].s2=j;
			map[idx].idx=idx;
		}
		std::sort(map, map+map_size,constellation_lt_p);
	}
	~constellation_db() {
		free(map);
		delete results;
		delete stars;
		delete orig_results;
		orig_stars=NULL;
	}
};

#endif
