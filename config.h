#define PI		   3.14159265358979323846  /* pi */
#define TWOPI		6.28318530717958647693

int IMG_X,IMG_Y,MAX_FALSE_STARS,REQUIRED_STARS;
float DEG_X,DEG_Y,PIXX_TANGENT,PIXY_TANGENT;
float PIXSCALE,POS_ERR_SIGMA,POS_VARIANCE;
float IMAGE_VARIANCE,EXPECTED_FALSE_STARS,MATCH_VALUE;
float PHOTONS,BRIGHT_THRESH;
float TODO_MAXFOV_D2;
float TODO_MINFOV_D2;

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

 
