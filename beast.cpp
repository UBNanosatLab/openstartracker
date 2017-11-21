;//#include <math.h>	   /* sqrt */
//#include <stdlib.h>
//
#include "beast.h"

namespace beast {
	int* imgdb::get_mask(float db_max_variance) {
		int *img_mask=(int*)malloc(IMG_X*IMG_Y*sizeof(img_mask[0]));
		memset(img_mask, -1, IMG_X*IMG_Y*sizeof(img_mask[0]));
		/* generate image mask */
		for (int id=0;id<star_map_size;id++){
			/* assume the dimmest possible star since we dont know the brightness of the other image */
			float sigma_sq=star_map[id].sigma_sq+db_max_variance;
			float maxdist_sq=-sigma_sq*(log(sigma_sq)+MATCH_VALUE);
			float maxdist=sqrt(maxdist_sq);
			int xmin=star_map[id].px-maxdist-1;
			int xmax=star_map[id].px+maxdist+1;
			int ymin=star_map[id].py-maxdist-1;
			int ymax=star_map[id].py+maxdist+1;
			
			if(xmax>IMG_X/2) xmax=IMG_X/2;
			if(xmin<-IMG_X/2)xmin=-IMG_X/2;
			if(ymax>IMG_Y/2) ymax=IMG_Y/2;
			if(ymin<-IMG_Y/2)ymin=-IMG_Y/2;
			for(int i=xmin;i<xmax;i++) for (int j=ymin;j<ymax;j++) {
				float a=((float)i-star_map[id].px);
				if (a<-0.5) a+=1.0;/* use whichever corner of the pixel gives the best score */
				float b=((float)j-star_map[id].py);
				if (b<-0.5) b+=1.0;
				
				int x=i+IMG_X/2;
				int y=j+IMG_Y/2;
				float score=(maxdist_sq-(a*a+b*b))/(2*sigma_sq);
			
				if (score>0) {
					/* has this pixel already been assigned to a different star? */
					int id2=img_mask[x+y*IMG_X];
					if (id2!=-1){
						float sigma_sq2=star_map[id2].sigma_sq+db_max_variance;
						float maxdist_sq2=-sigma_sq2*(log(sigma_sq2)+MATCH_VALUE);
						float px2=star_map[id2].px;
						float py2=star_map[id2].py;
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
	
	struct constellation_score *c_scores;
	int c_scores_size;
	float *winner_scores;
	//winner_id_map[new]=old
	int32_t *winner_id_map;
	float p_match;
	void add_score(constellation *db_const, constellation_score *cs, int32_t *img_mask, stardb *db, stardb *img){
		cs->db_id1=db_const->s1;
		cs->db_id2=db_const->s2;
		cs->id_map=(int32_t *)malloc(sizeof(int32_t)*img->star_map_size);
		cs->scores=(float *)malloc(sizeof(float)*img->star_map_size);
		for (int i=0;i<img->star_map_size;i++) cs->id_map[i]=-1;
		for (int i=0;i<img->star_map_size;i++) cs->scores[i]=0.0;
	
		cs->totalscore=log(EXPECTED_FALSE_STARS/(IMG_X*IMG_Y))*(2*img->star_map_size);
		//for(int32_t i=db_const->firststar;i<smax;i++){
		//	int32_t o=db->db_stars[i];
		db->kdsearch(R11,R12,R13,TODO_MAXFOV_D2,BRIGHT_THRESH);
		for(int32_t i=0;i<db->kdresults_size;i++){//
			int32_t o=db->kdresults[i];
			star s=db->star_map[o];
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
				float sigma_sq=img->star_map[n].sigma_sq+db->max_variance;
				float maxdist_sq=-sigma_sq*(log(sigma_sq)+MATCH_VALUE);
				float a=(px-img->star_map[n].px);
				float b=(py-img->star_map[n].py);
				float score = (maxdist_sq-(a*a+b*b))/(2*sigma_sq);
				/* only match the closest star */
				if (score>cs->scores[n]){
					cs->id_map[n]=o;
					cs->scores[n]=score;
				}
			}
		}
		DB->undo_kdsearch();
		for(int n=0;n<img->star_map_size;n++) {
			cs->totalscore+=cs->scores[n];
		}
		
	}
	void find(stardb *db, imgdb *img) {
		c_scores=NULL;
		c_scores_size=0;
		/* Do we have enough stars? */
		if (db->star_map_size<2||img->star_map_size<2) return;

		winner_id_map=(int *)malloc(sizeof(int)*img->star_map_size);
		winner_scores=(float *)malloc(sizeof(float)*img->star_map_size);
		
		/* allocate space for image mask */
		int32_t *img_mask = img->get_mask(db->max_variance);
		for (int n=0;n<img->db_map_size;n++) {
			constellation lb=img->db_map[n];
			constellation ub=img->db_map[n];
			lb.p-=POS_ERR_SIGMA*PIXSCALE*sqrt(img->star_map[lb.s1].sigma_sq+img->star_map[lb.s2].sigma_sq+2*db->max_variance);
			ub.p+=POS_ERR_SIGMA*PIXSCALE*sqrt(img->star_map[ub.s1].sigma_sq+img->star_map[ub.s2].sigma_sq+2*db->max_variance);
			
			constellation *lower=std::lower_bound (db->db_map, db->db_map+db->db_map_size, lb,constellation_lt_p);	
			constellation *upper=std::upper_bound (db->db_map, db->db_map+db->db_map_size, ub,constellation_lt_p); 
			//if (lower->p<lb.p) lower++;
			upper--;
			//if (upper->idx==db->db_map_size) 
			int lower_idx=lower->idx;
			int upper_idx=upper->idx;
			if (lower_idx<=upper_idx) c_scores=(struct constellation_score*)realloc(c_scores,sizeof(struct constellation_score)*(c_scores_size+(upper_idx-lower_idx+1)*2));
			for (int o=lower_idx;o<=upper_idx;o++) {
				int32_t db_idx1=db->db_map[o].s1;
				int32_t db_idx2=db->db_map[o].s2;
				int32_t img_idx1=img->db_map[n].s1;
				int32_t img_idx2=img->db_map[n].s2;
				
				star db_s1=db->star_map[db_idx1];
				star db_s2=db->star_map[db_idx2];
				star img_s1=img->star_map[img_idx1];
				star img_s2=img->star_map[img_idx2];
				
				/* try both orderings of stars */
				weighted_triad(db_s1,db_s2,img_s1,img_s2);
				c_scores[c_scores_size].img_id1=img_idx1;
				c_scores[c_scores_size].img_id2=img_idx2;
				add_score(&(db->db_map[o]),&c_scores[c_scores_size],img_mask,db,img);
				c_scores_size++;
				
				weighted_triad(db_s1,db_s2,img_s2,img_s1);
				c_scores[c_scores_size].img_id1=img_idx2;
				c_scores[c_scores_size].img_id2=img_idx1;
				add_score(&(db->db_map[o]),&c_scores[c_scores_size],img_mask,db,img);
				c_scores_size++;
			}
		}
		std::sort(c_scores, c_scores+c_scores_size);
		for (int i=0;i<img->star_map_size;i++) {winner_id_map[i]=-1;winner_scores[i]=0.0f;}
		if (c_scores_size>0) {
			for(int n=0;n<img->star_map_size;n++) {
				int o=c_scores[0].id_map[n];
				if (o!=-1){
					winner_scores[img->star_map[n].star_idx]=c_scores[0].scores[n];
					winner_id_map[img->star_map[n].star_idx]=db->star_map[o].id;
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
			weighted_triad(db->star_map[db_id1],db->star_map[db_id2],img->star_map[img_id1],img->star_map[img_id2]);
			for(int i=1;i<c_scores_size;i++) {
				if (c_scores[i].id_map[img_id1]!=db_id1&&c_scores[i].id_map[img_id2]!=db_id2){
					p_match+=exp(c_scores[i].totalscore-bestscore);
				}
			}
			/* TODO: WRONG BAD WRONG BAD WRONG BAD */
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
	beast::imgdb* img=new beast::imgdb;
	for(size_t i = 0; i < length; i++)
	{
		img->add_star(spikes[3*i]-beast::IMG_X/2.0,-(spikes[3*i+1]-beast::IMG_Y/2.0),spikes[3*i+2]);
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
