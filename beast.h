#include <stdint.h>
#include <math.h>       /* sqrt */

#define PI           3.14159265358979323846  /* pi */
#define TWOPI        6.28318530717958647693

struct constellation {
	float p;
	int32_t s1;
	int32_t s2;
	int32_t numstars;
	int32_t firststar;
	int32_t idx;
	bool operator< (const constellation* c) const { return p < c->p; }
	bool operator< (const constellation c) const { return p < c.p; }
};

struct star {
	float x;
	float y;
	float z;
	float photons;
	int32_t star_idx;
	int32_t photon_idx;
	int32_t id;

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
struct  constellation_score {
	float totalscore;
	int32_t db_id1,db_id2;
	int32_t img_id1,img_id2;
	int *id_map; /* Usage: id_map[newstar]=oldstar */
	float *scores;
	bool operator< (const constellation_score* c) { return totalscore > c->totalscore; }
	bool operator< (const constellation_score c) { return totalscore > c.totalscore; }
};
int IMG_X,IMG_Y,MAX_FALSE_STARS;
float DEG_X,DEG_Y,PIXX_TANGENT,PIXY_TANGENT;
float PIXSCALE,POS_ERR_SIGMA,POS_VARIANCE;
float IMAGE_VARIANCE,EXPECTED_FALSE_STARS,MATCH_VALUE;
float PHOTONS;

int NUMCONST,NUMSTARS,STARTABLE;

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
	DEG_X=atof(getenv("DEG_X"));
	DEG_Y=atof(getenv("DEG_Y"));
	POS_ERR_SIGMA=atof(getenv("POS_ERR_SIGMA"));
	PIXSCALE=atof(getenv("PIXSCALE"));
	POS_VARIANCE=atof(getenv("POS_VARIANCE"));/* sigma_r^2 */
	IMAGE_VARIANCE=atof(getenv("IMAGE_VARIANCE"));/* lambda */
	/* BRIGHT_THRESH=atof(getenv("BRIGHT_THRESH")); */ /* mmin */
	EXPECTED_FALSE_STARS=atof(getenv("EXPECTED_FALSE_STARS"));/* pfalse */
	/* need to do something with this for lost in space */
	MAX_FALSE_STARS=atoi(getenv("MAX_FALSE_STARS"));/* >10 is slow */
	PHOTONS=atoi(getenv("PHOTONS"));
	MATCH_VALUE=4*log(EXPECTED_FALSE_STARS/(IMG_X*IMG_Y))+log(2*PI);/* base */

	PIXX_TANGENT=2*tan(DEG_X*PI/(180*2))/IMG_X;
	PIXY_TANGENT=2*tan(DEG_Y*PI/(180*2))/IMG_Y;

	stream = fopen("calibration/dbsize.txt", "r");
	if (stream == NULL) exit(EXIT_FAILURE);
	while ((read = getline(&line, &len, stream)) != -1) {
		putenv(strcpy((char *)malloc(sizeof(char) * len),line));
	}
	fclose(stream);
	
	NUMCONST=atoi(getenv("NUMCONST"));
	NUMSTARS=atoi(getenv("NUMSTARS"));
	STARTABLE=atoi(getenv("STARTABLE"));
}

struct stardb {
	struct constellation* db_map;
	int *db_stars;
	struct star *star_map;
	int db_map_size;
	int db_stars_size;
	int star_map_size;
	float max_variance;

	stardb() {
		db_map=NULL;
		db_stars=NULL;
		star_map=NULL;
		db_map_size=0;
		db_stars_size=0;
		star_map_size=0;
		
		max_variance=0.0;
	}
	~stardb() {
		free(db_map);
		free(db_stars);
		free(star_map);
	}
};
struct imgdb: public stardb {	
	void add_star(float px, float py, float mag);
	void add_star(float x, float y, float z, float mag);
	void gendb();
	int* get_mask(float db_max_variance);
};
struct beastdb: public stardb {	
	int fd;
	beastdb() {
		
		db_map_size=NUMCONST;
		db_stars_size=NUMSTARS;
		star_map_size=STARTABLE;
		max_variance=POS_VARIANCE;
		
		const char *filepath = "beastdb.bin";
		fd = open(filepath, O_RDONLY);
		if (fd == -1) exit(EXIT_FAILURE);

		/* Now the file is ready to be mmapped. */
		db_map = (constellation*)mmap(NULL, db_map_size*sizeof(constellation) + star_map_size*sizeof(int), PROT_READ, MAP_SHARED | MAP_POPULATE, fd, 0);
		if (db_map == MAP_FAILED) exit(EXIT_FAILURE);
		
		db_stars=(int*)(&db_map[db_map_size]);
		
		//TODO: do this via mmap using read() and write() methods
		FILE *stream = fopen("calibration/stars.txt", "r");
		if (stream == NULL) exit(EXIT_FAILURE);

		ssize_t read;
		char *line = NULL;
		size_t len = 0;
		star_map = (star*)malloc(db_stars_size*sizeof(star));
		for(int i=0;i<db_stars_size;i++){
			if ((read = getline(&line, &len, stream)) != -1) {
				star_map[i].id=atoi(strtok(line," "));
				star_map[i].star_idx=i;
				float mag=atof(strtok(NULL," "));
				star_map[i].photons=PHOTONS*powf(10.0,-mag/2.5);
				
				star_map[i].x=atof(strtok(NULL," "));
				star_map[i].y=atof(strtok(NULL," "));
				star_map[i].z=atof(strtok(NULL," "));
				
				star_map[i].sigma_sq=POS_VARIANCE;
			}
		}
		fclose(stream);
	}
	~beastdb() {
		if (munmap(db_map, db_map_size*sizeof(constellation) + star_map_size*sizeof(int)) == -1) exit(EXIT_FAILURE);
		close(fd);
		db_map=NULL;
		db_stars=NULL;
	}
};

beastdb* DB;
void load_db() {
	load_config();
	DB=new beastdb;
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

