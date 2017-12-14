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
#endif
