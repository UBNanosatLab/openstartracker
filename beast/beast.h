#ifndef BEAST_H
#define BEAST_H

#include "constellations.h"
#include <float.h>

struct  match_result {
//TODO: private:
	constellation_pair match;
	
	//eci to body (body=R*eci)
	float R11,R12,R13;
	float R21,R22,R23;
	float R31,R32,R33;
private:
	star_fov *img_mask;
	int *map; /* Usage: map[imgstar]=dbstar */
	size_t map_size;
	constellation* db_const;
	constellation_db *db,*img;
public:
	
	/**
	* @brief TODO
	* @param db_
	* @param img_
	* @param img_mask_
	*/
	match_result(constellation_db *db_, constellation_db *img_, star_fov *img_mask_) {
		DBG_MATCH_RESULT_COUNT++;
		DBG_PRINT("DBG_MATCH_RESULT_COUNT++ %d\n",DBG_MATCH_RESULT_COUNT);
		db=db_;
		img=img_;
		img_mask=img_mask_;
		map_size=img->stars->size();
		map=(int *)malloc(sizeof(map[0])*map_size);
		match.totalscore=-FLT_MAX;
		
	}
	~match_result() {
		DBG_MATCH_RESULT_COUNT--;
		DBG_PRINT("DBG_MATCH_RESULT_COUNT-- %d\n",DBG_MATCH_RESULT_COUNT);
		free(map);
	}
	size_t size() {return map_size;};
	/**
	* @brief TODO
	* @param db_const_
	* @param img_const_
	*/
	void init(constellation &db_const_, constellation &img_const_) {
		db_const=&db_const_;
		
		match.img_s1=img_const_.s1;
		match.img_s2=img_const_.s2;
		match.db_s1=db_const_.s1;
		match.db_s2=db_const_.s2;
	}
	
	/**
	* @brief TODO
	* @param c
	*/
	void copy_over(match_result *c) {
		assert(c->db==db);
		assert(c->img==img);
		assert(c->img_mask==img_mask);
		
		c->match=match;
		c->db_const=db_const;
		
		c->R11=R11,c->R12=R12,c->R13=R13;
		c->R21=R21,c->R22=R22,c->R23=R23;
		c->R31=R31,c->R32=R32,c->R33=R33;
		
		memcpy(c->map, map, sizeof(map[0])*map_size);
	}
	/**
	* @brief TODO
	* @param m
	* @return 
	*/
	int related(constellation_pair &m) {
		
		if (match.totalscore==-FLT_MAX || m.totalscore==-FLT_MAX) return 0;
		return (map[m.img_s1]==m.db_s1 && map[m.img_s2]==m.db_s2)?1:0;
	}
	/**
	* @brief TODO
	*/
	void search() {if (db->results->is_kdsorted()) db->results->kdsearch(R11,R21,R31,MAXFOV/2,THRESH_FACTOR*IMAGE_VARIANCE);}
	/**
	* @brief TODO
	*/
	void clear_search() {if (db->results->is_kdsorted()) db->results->clear_kdresults();}
	/**
	* @brief TODO
	*/
	void compute_score() {
		//TODO: figure out where 2*map_size came from
		match.totalscore=log(1.0/(IMG_X*IMG_Y))*(2*map_size);
		float* scores=(float *)malloc(sizeof(float)*map_size);
		for (size_t i=0;i<map_size;i++) {
			map[i]=-1;
			scores[i]=0.0;
		}
		for(size_t i=0;i<db->results->r_size();i++) {
			star *s=&(db->results->map[db->results->kdresults[i]]);
			int o=s->star_idx;
			float x=s->x*R11+s->y*R21+s->z*R31;
			float y=s->x*R12+s->y*R22+s->z*R32;
			float z=s->x*R13+s->y*R23+s->z*R33;
			float px=y/(x*PIXX_TANGENT);
			float py=z/(x*PIXY_TANGENT);
			
			int n=img_mask->get_id(px,py);
			if (n>=0) {
				float score = img_mask->get_score(n,px,py);
				if (score>scores[n]){/* only match the closest star */
					map[n]=o;
					scores[n]=score;
				}
			}
		}
		for(size_t n=0;n<map_size;n++) {
			match.totalscore+=scores[n];
		}
		free(scores);
	}
	/**
	* @return matching stars from db, in order of star_idx
	*/
	star_db* from_match() {
		if (match.totalscore==-FLT_MAX) return NULL;
		
		star_db* s = img->stars->copy();
		s->max_variance=db->stars->max_variance;
		for(size_t n=0;n<map_size;n++) {
			//catalog matching
			if (map[n]!=-1) {
				s->get_star(img->stars->get_star(n)->star_idx)[0]=db->stars->get_star(map[n])[0];
			} else {
				s->get_star(img->stars->get_star(n)->star_idx)->id=-1;
			}
		}
		return s;
	}

	/**
	* @brief weighted_triad results
	* see https://en.wikipedia.org/wiki/Triad_method 
	* and http://nghiaho.com/?page_id=846
	* 
	* when compiled, this section contains roughly 430 floating point operations
	* according to https://www.karlrupp.net/2016/02/gemm-and-stream-results-on-intel-edison
	* we can perform >250 MFLOPS with doubles, and >500 MFLOPS with floats
	*/
	void weighted_triad() {
		star *db_s1=db->stars->get_star(match.db_s1);
		star *db_s2=db->stars->get_star(match.db_s2);
		star *img_s1=img->stars->get_star(match.img_s1);
		star *img_s2=img->stars->get_star(match.img_s2);
		
		/* v=A*w */
		float wa1=db_s1->x,wa2=db_s1->y,wa3=db_s1->z;
		float wb1=db_s2->x,wb2=db_s2->y,wb3=db_s2->z;
		float va1=img_s1->x,va2=img_s1->y,va3=img_s1->z;
		float vb1=img_s2->x,vb2=img_s2->y,vb3=img_s2->z;
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
		float weightA=1.0/(db_s1->sigma_sq+img_s1->sigma_sq);
		float weightB=1.0/(db_s2->sigma_sq+img_s2->sigma_sq);

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
		R21=cz*sx*sy - cx*sz;
		R31=sx*sz + cx*cz*sy;
		
		R12=cy*sz;
		R22=cx*cz + sx*sy*sz;
		R32=cx*sy*sz - cz*sx;
		
		R13=-sy;
		R23=cy*sx;
		R33=cx*cy;
	}
	void DBG_(const char *s) {
		DBG_PRINT("%s\n",s);
		DBG_PRINT("%f\t%f\t%f\n", R11,R12,R13);
		DBG_PRINT("%f\t%f\t%f\n", R21,R22,R23);
		DBG_PRINT("%f\t%f\t%f\n", R31,R32,R33);
		
		db->DBG_("DB");
		img->DBG_("IMG");
		DBG_PRINT("map_size=%lu\n", map_size);
		for (size_t i=0; i<map_size; i++) {
			DBG_PRINT("map[%lu]=%d\n",i,map[i]);
		}
	}
	/**
	* @brief TODO
	*/
	void print_ori() {
		fprintf(stderr,"DEC=%f\n",fmod(360+asin(R31)* 180 / PI,360));
		fprintf(stderr,"RA=%f\n",fmod(360+atan2(R21,R11)* 180 / PI,360));
		fprintf(stderr,"ORIENTATION=%f\n",-atan2(R32,R33)* 180 / PI);
	}
};

struct db_match {
private:
	constellation_pair *c_pairs;
	size_t c_pairs_size;
	star_fov *img_mask;
public:
	float p_match;
	match_result *winner;
	
		
	/**
	* @brief TODO
	* @param db
	* @param img
	*/
	db_match(constellation_db *db, constellation_db *img) {
		DBG_DB_MATCH_COUNT++;
		DBG_PRINT("DBG_DB_MATCH_COUNT++ %d\n",DBG_DB_MATCH_COUNT);
		winner=NULL;
		img_mask=NULL;
		c_pairs=NULL;
		c_pairs_size=0;
		p_match=0.0;
//TODO - set size to 4 in python
		if (db->stars->size()<3||img->stars->size()<3) return;
		img_mask = new star_fov(img->stars,db->stars->max_variance);
		
		//find stars
		match_result *m=new match_result(db, img, img_mask);
		winner=new match_result(db, img, img_mask);
		for (size_t n=0;n<img->map_size;n++) {
			constellation lb=img->map[n];
			constellation ub=img->map[n];
			lb.p-=POS_ERR_SIGMA*PIXSCALE*sqrt(img->stars->get_star(lb.s1)->sigma_sq+img->stars->get_star(lb.s2)->sigma_sq+2*db->stars->max_variance);
			ub.p+=POS_ERR_SIGMA*PIXSCALE*sqrt(img->stars->get_star(ub.s1)->sigma_sq+img->stars->get_star(ub.s2)->sigma_sq+2*db->stars->max_variance);
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

				#define ADD_SCORE\
					m->compute_score();\
					if (m->match.totalscore>winner->match.totalscore) {\
						if (winner->match.totalscore!=-FLT_MAX) c_pairs[c_pairs_size++]=winner->match;\
						m->copy_over(winner);\
					} else c_pairs[c_pairs_size++]=m->match;

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
			for (size_t idx=0; idx<c_pairs_size;idx++) {
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
