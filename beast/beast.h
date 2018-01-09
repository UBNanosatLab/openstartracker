#ifndef BEAST_H
#define BEAST_H

#include "constellations.h"
#include <float.h>

struct  match_result {
	constellation_db *db,*img;

	constellation* db_const;

	star_fov *img_mask;
	constellation_pair match;
	int *map; /* Usage: map[imgstar]=dbstar */
	int map_size;
	
	//eci to body (body=R*eci)
	float R11,R12,R13;
	float R21,R22,R23;
	float R31,R32,R33;
	
	match_result(constellation_db *db_, constellation_db *img_, star_fov *img_mask_) {
		DBG_MATCH_RESULT_COUNT++;
		DBG_PRINT("DBG_MATCH_RESULT_COUNT++ %d\n",DBG_MATCH_RESULT_COUNT);
		db=db_;
		img=img_;
		img_mask=img_mask_;
		map_size=img->stars->map_size;
		map=(int *)malloc(sizeof(int)*map_size);
		match.totalscore=-FLT_MAX;
		
	}
	~match_result() {
		DBG_MATCH_RESULT_COUNT--;
		DBG_PRINT("DBG_MATCH_RESULT_COUNT-- %d\n",DBG_MATCH_RESULT_COUNT);
		free(map);
	}
	void init(constellation &db_const_, constellation &img_const_) {
		db_const=&db_const_;
		
		match.img_s1=img_const_.s1;
		match.img_s2=img_const_.s2;
		match.db_s1=db_const_.s1;
		match.db_s2=db_const_.s2;
	}
	
	void copy_over(match_result *c) {
		assert(c->db==db);
		assert(c->img==img);
		assert(c->img_mask==img_mask);
		
		c->match=match;
		c->db_const=db_const;
		
		c->R11=R11,c->R12=R12,c->R13=R13;
		c->R21=R21,c->R22=R22,c->R23=R23;
		c->R31=R31,c->R32=R32,c->R33=R33;
		
		memcpy(c->map, map, sizeof(int)*map_size);
	}
	int related(constellation_pair &m) {
		
		if (match.totalscore==-FLT_MAX || m.totalscore==-FLT_MAX) return 0;
		return (map[m.img_s1]==m.db_s1 && map[m.img_s2]==m.db_s2)?1:0;
	}
	void search() {if (db->stars->kdsorted==1) db->results->kdsearch(R11,R12,R13,MAXFOV/2,THRESH_FACTOR*IMAGE_VARIANCE);}
	void clear_search() {if (db->stars->kdsorted==1) db->results->clear_kdresults();}
	void compute_score() {
		//TODO: figure out where 2*map_size came from
		match.totalscore=log(1.0/(IMG_X*IMG_Y))*(2*map_size);
		float* scores=(float *)malloc(sizeof(float)*map_size);
		for (int i=0;i<map_size;i++) {
			map[i]=-1;
			scores[i]=0.0;
		}
		for(int i=0;i<db->results->kdresults_size;i++) {
			int o=db->results->kdresults[i];
			star *s=&(db->stars->map[o]);
			float x=s->x*R11+s->y*R12+s->z*R13;
			float y=s->x*R21+s->y*R22+s->z*R23;
			float z=s->x*R31+s->y*R32+s->z*R33;
			float px=IMG_ROTATION*y/(x*PIXX_TANGENT);
			float py=IMG_ROTATION*z/(x*PIXY_TANGENT);
			
			int nx=(int)(px+IMG_X/2.0f);
			if (nx==-1) nx++;
			else if (nx==IMG_X) nx--;

			int ny=(int)(py+IMG_Y/2.0f);
			if (ny==-1) ny++;
			else if (ny==IMG_Y) ny--;
			int n=-1;
			if (nx>=0&&nx<IMG_X&&ny>=0&&ny<IMG_Y) n=img_mask->mask[nx+ny*IMG_X];
			if (n<-1) n=img_mask->resolve_id(n,px,py);
			if (n>=0) {
				float score = img_mask->get_score(n,px,py);
				if (score>scores[n]){/* only match the closest star */
					map[n]=o;
					scores[n]=score;
				}
			}
		}
		for(int n=0;n<map_size;n++) {
			match.totalscore+=scores[n];
		}
		free(scores);
	}
	//return matching stars from db, in order of star_idx
	star_db* from_match() {
		if (match.totalscore==-FLT_MAX) return NULL;
		
		star_db* s = new star_db;
		s->map_size=map_size;
		s->max_variance=db->stars->max_variance;
		s->map=(star *)malloc(sizeof(s->map[0])*map_size);
		memset(s->map,-1,sizeof(s->map[0])*map_size);
		for(int n=0;n<map_size;n++) {
			//catalog matching
			int o=map[n];
			if (o!=-1) {
				int img_star_idx=img->stars->map[n].star_idx;
				s->map[img_star_idx]=db->stars->map[o];
			}
		}
		return s;
	}
	
	/* weighted_triad results */

	/* see https://en.wikipedia.org/wiki/Triad_method */
	/* and http://nghiaho.com/?page_id=846 */
	/* returns match results */

	/* when compiled, this section contains roughly 430 floating point operations */
	/* according to https://www.karlrupp.net/2016/02/gemm-and-stream-results-on-intel-edison/ */
	/* we can perform >250 MFLOPS with doubles, and >500 MFLOPS with floats */
	
	void weighted_triad() {
		star db_s1=db->stars->map[match.db_s1];
		star db_s2=db->stars->map[match.db_s2];
		star img_s1=img->stars->map[match.img_s1];
		star img_s2=img->stars->map[match.img_s2];
		
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
	void DBG_(const char *s) {
		DBG_PRINT("%s\n",s);
		DBG_PRINT("%f\t%f\t%f\n", R11,R12,R13);
		DBG_PRINT("%f\t%f\t%f\n", R21,R22,R23);
		DBG_PRINT("%f\t%f\t%f\n", R31,R32,R33);
		
		db->DBG_("DB");
		img->DBG_("IMG");
		DBG_PRINT("map_size=%d\n", map_size);
		for (int i=0; i<map_size; i++) {
			DBG_PRINT("map[%d]=%d\n",i,map[i]);
		}
	}
	void print_ori() {
		fprintf(stderr,"DEC=%f\n",asin(R13)* 180 / PI);
		fprintf(stderr,"RA=%f\n",atan2(R12,R11)* 180 / PI);
		fprintf(stderr,"ORIENTATION=%f\n",fmod(360-atan2(R23,R33)* 180 / PI ,360));
	}
};

struct db_match {
	float p_match;
	
	match_result *winner;
	constellation_pair *c_pairs;
	int c_pairs_size;
	star_fov *img_mask;
	
	#define ADD_SCORE\
		m->compute_score();\
		if (m->match.totalscore>winner->match.totalscore) {\
			if (winner->match.totalscore!=-FLT_MAX) c_pairs[c_pairs_size++]=winner->match;\
			m->copy_over(winner);\
		} else c_pairs[c_pairs_size++]=m->match;
		
	db_match(constellation_db *db, constellation_db *img) {
		DBG_DB_MATCH_COUNT++;
		DBG_PRINT("DBG_DB_MATCH_COUNT++ %d\n",DBG_DB_MATCH_COUNT);
		winner=NULL;
		img_mask=NULL;
		c_pairs=NULL;
		c_pairs_size=0;
		p_match=0.0;

		if (db->stars->map_size<3||img->stars->map_size<3) return;
		img_mask = new star_fov(img->stars,db->stars->max_variance);
		
		//find stars
		match_result *m=new match_result(db, img, img_mask);
		winner=new match_result(db, img, img_mask);
		for (int n=0;n<img->map_size;n++) {
			constellation lb=img->map[n];
			constellation ub=img->map[n];
			lb.p-=POS_ERR_SIGMA*PIXSCALE*sqrt(img->stars->map[lb.s1].sigma_sq+img->stars->map[lb.s2].sigma_sq+2*db->stars->max_variance);
			ub.p+=POS_ERR_SIGMA*PIXSCALE*sqrt(img->stars->map[ub.s1].sigma_sq+img->stars->map[ub.s2].sigma_sq+2*db->stars->max_variance);
			constellation *lower=std::lower_bound (db->map, db->map+db->map_size, lb,constellation_lt_p);	
			constellation *upper=std::upper_bound (db->map, db->map+db->map_size, ub,constellation_lt_p);
			//rewind upper & do sanity checks
			if (db->map>=upper--) continue;
			if (db->map+db->map_size<=lower) continue;
			if (lower->idx<=upper->idx) {
				c_pairs=(struct constellation_pair*)realloc(c_pairs,sizeof(struct constellation_pair)*(c_pairs_size+(upper->idx-lower->idx+1)*2));
			}
			for (int o=lower->idx;o<=upper->idx;o++) {
				m->init(db->map[o],img->map[n]);
				m->weighted_triad();
				m->search();
				ADD_SCORE
				/* try both orderings of stars */
				m->match.flip();
				m->weighted_triad();
				ADD_SCORE
				m->clear_search();
			}
		}
		delete m;
		
		//calculate map
		if (winner->match.totalscore!=-FLT_MAX) { //Did we even match?
			//calculate p_match
			p_match=1.0;
			for (int idx=0; idx<c_pairs_size;idx++) {
				if (!winner->related(c_pairs[idx])){
					p_match+=exp(c_pairs[idx].totalscore-winner->match.totalscore);
				}
			}
			p_match=1.0/p_match;
		}
	}
	
	~db_match() {
		DBG_DB_MATCH_COUNT--;
		DBG_PRINT("DBG_DB_MATCH_COUNT-- %d\n",DBG_DB_MATCH_COUNT);
		delete winner;
		delete img_mask;
		free(c_pairs);
	}
};
#endif
