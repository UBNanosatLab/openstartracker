;//#include <math.h>	   /* sqrt */
//#include <stdlib.h>
//
#include "beast.h"

namespace beast {

	struct constellation_score *c_scores;
	int c_scores_size;
	float *winner_scores;
	//winner_id_map[new]=old
	int32_t *winner_id_map;
	float p_match;
	void add_score(constellation *db_const, constellation_score *cs, int32_t *img_mask, constellation_db *db, constellation_db *img){
		cs->db_id1=db_const->s1;
		cs->db_id2=db_const->s2;
		cs->id_map=(int32_t *)malloc(sizeof(int32_t)*img->stars->map_size);
		cs->scores=(float *)malloc(sizeof(float)*img->stars->map_size);
		for (int i=0;i<img->stars->map_size;i++) cs->id_map[i]=-1;
		for (int i=0;i<img->stars->map_size;i++) cs->scores[i]=0.0;
		
		//cs->totalscore=log(EXPECTED_FALSE_STARS/(IMG_X*IMG_Y))*(2*img->stars->map_size);
		cs->totalscore=0;
		//for(int32_t i=db_const->firststar;i<smax;i++){
		//	int32_t o=db->db_stars[i];
		//db->stars->kdsearch(R11,R12,R13,TODO_MAXFOV_D2,BRIGHT_THRESH,db_const->kdbox_min,db_const->kdbox_max,db_const->kdbox_dim);
		db->stars->kdsearch(R11,R12,R13,TODO_MAXFOV_D2,BRIGHT_THRESH);
		for(int32_t i=0;i<db->stars->kdresults_size;i++){
			int32_t o=db->stars->kdresults[i];
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
		DB->stars->undo_kdsearch();
		for(int n=0;n<img->stars->map_size;n++) {
			cs->totalscore+=cs->scores[n];
		}
		
	}
	void find(constellation_db *db, constellation_db *img) {
		c_scores=NULL;
		c_scores_size=0;
		/* Do we have enough stars? */
		if (db->stars->map_size<2||img->stars->map_size<2) return;

		winner_id_map=(int *)malloc(sizeof(int)*img->stars->map_size);
		winner_scores=(float *)malloc(sizeof(float)*img->stars->map_size);
		
		int32_t *img_mask = img->stars->get_img_mask(db->stars->max_variance);
		
		for (int n=0;n<img->map_size;n++) {
			constellation lb=img->map[n];
			constellation ub=img->map[n];
			lb.p-=POS_ERR_SIGMA*PIXSCALE*sqrt(img->stars->map[lb.s1].sigma_sq+img->stars->map[lb.s2].sigma_sq+2*db->stars->max_variance);
			ub.p+=POS_ERR_SIGMA*PIXSCALE*sqrt(img->stars->map[ub.s1].sigma_sq+img->stars->map[ub.s2].sigma_sq+2*db->stars->max_variance);
			constellation *lower=std::lower_bound (db->map, db->map+db->map_size, lb,constellation_lt_p);	
			constellation *upper=std::upper_bound (db->map, db->map+db->map_size, ub,constellation_lt_p);
			//rewind by one
			upper--;
			if (lower->idx<=upper->idx) c_scores=(struct constellation_score*)realloc(c_scores,sizeof(struct constellation_score)*(c_scores_size+(upper->idx-lower->idx+1)*2));
			for (int o=lower->idx;o<=upper->idx;o++) {
				int32_t db_idx1=db->map[o].s1;
				int32_t db_idx2=db->map[o].s2;
				int32_t img_idx1=img->map[n].s1;
				int32_t img_idx2=img->map[n].s2;
				
				star db_s1=db->stars->map[db_idx1];
				star db_s2=db->stars->map[db_idx2];
				star img_s1=img->stars->map[img_idx1];
				star img_s2=img->stars->map[img_idx2];
				
				/* try both orderings of stars */
				weighted_triad(db_s1,db_s2,img_s1,img_s2);
				c_scores[c_scores_size].img_id1=img_idx1;
				c_scores[c_scores_size].img_id2=img_idx2;
				add_score(&(db->map[o]),&c_scores[c_scores_size],img_mask,db,img);
				c_scores_size++;
				
				weighted_triad(db_s1,db_s2,img_s2,img_s1);
				c_scores[c_scores_size].img_id1=img_idx2;
				c_scores[c_scores_size].img_id2=img_idx1;
				add_score(&(db->map[o]),&c_scores[c_scores_size],img_mask,db,img);
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
	//TODO: get rid of this
	void del_scores() {
		for (int i=0;i<c_scores_size; i++) {
			free(c_scores[i].scores);
			free(c_scores[i].id_map);
		}
		free(c_scores);
		free(winner_id_map);
		free(winner_scores);
	}
}
/***
 * Determines the HIP numbers of stars in a scene.
 * 
 * @param spikes The data of the spikes in the scene. x, y, and magnitude interleaved.
 * @param result The output array that should be filled with the HIP numbers.
 * @param length The number of spikes in the scene.
 */
void star_id(double spikes[], int result[], size_t length)
{
	beast::constellation_db* img=new beast::constellation_db;
	for(size_t i = 0; i < length; i++)
	{
		img->stars->add_star(spikes[3*i]-beast::IMG_X/2.0,-(spikes[3*i+1]-beast::IMG_Y/2.0),spikes[3*i+2]);
		result[i] = -1;
	}
	
	img->gendb_img();
	find(beast::DB,img);
	float p_match = beast::p_match;
	if (p_match>0.66) {
		for(size_t i = 0; i < length; i++) {
			result[i] = beast::winner_id_map[i];
		}
	}
	beast::del_scores();
	delete img;
}
