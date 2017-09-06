#include <math.h>       /* sqrt */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>	/* For O_RDWR */
#include <sys/mman.h>
#include <assert.h>
#include <stdint.h>


#define PI           3.14159265358979323846  /* pi */
#define TWOPI        6.28318530717958647693

struct constellation {
	float p;
	int32_t s1;
	int32_t s2;
	int32_t numstars;
	int32_t staridx;
	int32_t last_l;
	int32_t last_m;
	int32_t last_h;
	int32_t mmi_l;
	int32_t mmi_h;
};

struct star {
	float x;
	float y;
	float z;
	float mag;
	int starnum;
	int magnum;
	int id;

	float sigma_sq;
	float px;
	float py;
};

struct  constellation_score {
	float totalscore;
	int oldid1;
	int oldid2;
	unsigned char newid1;
	unsigned char newid2;
	
	//Usage: id_map[newstar]=oldstar
	int *id_map;
	float *scores;
};

//SWIG complains about these pointers. They're fine.
int *map;
struct constellation *constptr;
int *db_starids;
struct star *starptr;

int fd;
size_t dbsize;
int mapsize;


int PARAM,NUMCONST,NUMSTARS,STARTABLE;

int IMG_X,IMG_Y,MAX_FALSE_STARS;
float DEG_X,DEG_Y,PIXX_TANGENT,PIXY_TANGENT;
float ARC_ERR,POS_VARIANCE;
float IMAGE_VARIANCE,EXPECTED_FALSE_STARS,MATCH_VALUE;
float PHOTONS;

struct star newstars[256];
struct star oldstars[256];
uint8_t numnewstars=0;	
uint8_t numoldstars=0;	
struct constellation_score *c_scores;
int c_scores_size;
float *winner_scores;
//winner_id_map[new]=old
int *winner_id_map;

int addednewstars;
int addedoldstars;

unsigned char *img_mask;
//rotation matrix
float R11,R12,R13;
float R21,R22,R23;
float R31,R32,R33;

void load_db() {
	//load config
	FILE *stream;
	char *line = NULL;
	winner_scores=NULL;
	winner_id_map=NULL;
	c_scores_size=0;
	
	size_t len = 0;
	ssize_t read;

	stream = fopen("calibration/calibration.txt", "r");
	if (stream == NULL) exit(EXIT_FAILURE);
	while ((read = getline(&line, &len, stream)) != -1) {
		//Don't listen to valgrind. This is fine. Everything is fine.
		putenv(strcpy((char *)malloc(sizeof(char) * len),line));
	}
	fclose(stream);
	
	stream = fopen("calibration/dbsize.txt", "r");
	if (stream == NULL) exit(EXIT_FAILURE);
	while ((read = getline(&line, &len, stream)) != -1) {
		//Don't listen to valgrind. This is fine. Everything is fine.
		putenv(strcpy((char *)malloc(sizeof(char) * len),line));
	}
	fclose(stream);
	
	IMG_X=atoi(getenv("IMG_X"));
	IMG_Y=atoi(getenv("IMG_Y"));
	DEG_X=atof(getenv("DEG_X"));
	DEG_Y=atof(getenv("DEG_Y"));
	ARC_ERR=atof(getenv("ARC_ERR"));
	POS_VARIANCE=atof(getenv("POS_VARIANCE"));//sigma_r^2
	IMAGE_VARIANCE=atof(getenv("IMAGE_VARIANCE"));//lambda
	//BRIGHT_THRESH=atof(getenv("BRIGHT_THRESH"));//mmin
	EXPECTED_FALSE_STARS=atof(getenv("EXPECTED_FALSE_STARS"));//pfalse
	//need to do something with this for lost in space
	MAX_FALSE_STARS=atoi(getenv("MAX_FALSE_STARS"));//>10 is slow
	PHOTONS=atoi(getenv("PHOTONS"));//>10 is slow
	MATCH_VALUE=4*log(EXPECTED_FALSE_STARS/(IMG_X*IMG_Y))+log(2*PI);//base
	
	PIXX_TANGENT=2*tan(DEG_X*PI/(180*2))/IMG_X;
	PIXY_TANGENT=2*tan(DEG_Y*PI/(180*2))/IMG_Y;
	
	PARAM=atoi(getenv("PARAM"));
	NUMCONST=atoi(getenv("NUMCONST"));
	NUMSTARS=atoi(getenv("NUMSTARS"));
	STARTABLE=atoi(getenv("STARTABLE"));
	
	mapsize=PARAM;
	int s_offset=mapsize*sizeof(int)+ NUMCONST*sizeof(struct constellation);
	dbsize = s_offset + STARTABLE*sizeof(int);
	
	/* Open a file for writing.
	 *  - Creating the file if it doesn't exist.
	 *  - Truncating it to 0 size if it already exists. (not really needed)
	 *
	 * Note: "O_WRONLY" mode is not sufficient when mmaping.
	 */

	const char *filepath = "beastdb.bin";
	fd = open(filepath, O_RDONLY);
	if (fd == -1) {
		perror("Error opening file for writing");
		exit(EXIT_FAILURE);
	}

	// Now the file is ready to be mmapped.

	map = (int*)mmap(NULL, dbsize, PROT_READ, MAP_SHARED | MAP_POPULATE, fd, 0);
	if (map == MAP_FAILED)
	{
		close(fd);
		perror("Error mmapping the file");
		exit(EXIT_FAILURE);
	}
	constptr=(struct constellation*)(&map[mapsize]);
	db_starids=(int*)(&map[s_offset/sizeof(int)]);
	
	stream = fopen("calibration/stars.txt", "r");
	if (stream == NULL) exit(EXIT_FAILURE);
	starptr = (struct star*)malloc(NUMSTARS*sizeof(struct star));
	for(int i=0;i<NUMSTARS;i++){
		if ((read = getline(&line, &len, stream)) != -1) {
			starptr[i].id=atoi(strtok(line," "));
			starptr[i].starnum=i;
			starptr[i].mag=atof(strtok(NULL," "));
			starptr[i].mag=PHOTONS*powf(10.0,-starptr[i].mag/2.5);;
			starptr[i].x=atof(strtok(NULL," "));
			starptr[i].y=atof(strtok(NULL," "));
			starptr[i].z=atof(strtok(NULL," "));
		}
	}
}

void unload_db() {
	// Don't forget to free the mmapped memory
	if (munmap(map, dbsize) == -1) {
			close(fd);
			perror("Error un-mmapping the file");
			exit(EXIT_FAILURE);
	}
	// Un-mmaping doesn't close the file, so we still need to do that.
	close(fd);
	free(starptr);
}

void c_scores_clear() {
	int i;
	for (i=0;i<c_scores_size; i++) {
		free(c_scores[i].scores);
		free(c_scores[i].id_map);
	}
	free(c_scores);
	c_scores=NULL;
	c_scores_size=0;
}
int compare_mag (const void *a, const void *b) {
	struct star *s1 =(struct star *)a;
	struct star *s2 =(struct star *)b;
	if (s1->mag > s2->mag) return 1;
	else if (s1->mag < s2->mag) return -1;
	else return 0;
}
int compare_totalscore (const void *a, const void *b) {
	struct constellation_score *cs1 =(struct constellation_score *)a;
	struct constellation_score *cs2 =(struct constellation_score *)b;
	if (cs1->totalscore > cs2->totalscore) return -1;
	else if (cs1->totalscore < cs2->totalscore) return 1;
	else return 0;
}

void add_star(float px, float py, float mag) {
	struct star s;

	float j=PIXX_TANGENT*px; //j=(y/x)
	float k=PIXY_TANGENT*py; //k=z/x
	s.x=1./sqrt(j*j+k*k+1);
	s.y=j*s.x;
	s.z=k*s.x;
	s.mag=PHOTONS*powf(10.0,-mag/2.5);
	s.id=-1;
	
	s.sigma_sq=IMAGE_VARIANCE/mag;
	s.px=px;
	s.py=py;

	//insert into list sorted by magnitude
	s.magnum=numnewstars++;
	if (s.magnum==0) addednewstars=0;
	s.starnum=addednewstars;
	addednewstars++;
	while (s.magnum>0&&compare_mag(&s,&newstars[s.magnum-1])>0){
		memcpy(&newstars[s.magnum],&newstars[s.magnum-1],sizeof(struct star));
		newstars[s.magnum].magnum++;
		s.magnum--;
	}
	memcpy(&newstars[s.magnum],&s,sizeof(struct star));
	if(numnewstars==255) numnewstars=254;
}
/* return sin(theta) where theta is the angle between vectors */
float dist3(float x1,float x2,float y1,float y2,float z1,float z2) {
	float a=x1*y2 - x2*y1;
	float b=x1*z2 - x2*z1;
	float c=y1*z2 - y2*z1;
	return sqrt(a*a+b*b+c*c);
}


//see https://en.wikipedia.org/wiki/Triad_method
//and http://nghiaho.com/?page_id=846
//returns match results

//Optimization for LIS:

//when compiled, this section contains roughly 430 floating point operations
//according to https://www.karlrupp.net/2016/02/gemm-and-stream-results-on-intel-edison/
//we can perform >250 MFLOPS with floats, and >500 MFLOPS with floats
//assuming 250 stars per pair * 100 pairs per image ~ .05 seconds

//tips to improve speed: 
//reduce max stars (x2-4)
//replace floats with floats(where 7 digits are enough) (x2-4)
//Perform triad method using both stars as pilot stars, and 
//set the rotation matrix to the weighted interpolation of the two

void weighted_triad(struct star old_s1,struct star old_s2,struct star new_s1,struct star new_s2,float variance){
	//v=A*w
	float wa1=old_s1.x,wa2=old_s1.y,wa3=old_s1.z;
	float wb1=old_s2.x,wb2=old_s2.y,wb3=old_s2.z;
	float va1=new_s1.x,va2=new_s1.y,va3=new_s1.z;
	float vb1=new_s2.x,vb2=new_s2.y,vb3=new_s2.z;
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
	
	//some of these are unused
	float A11=va1*wa1 + vaXvc1*waXwc1 + vc1*wc1;
	//float A12=va1*wa2 + vaXvc1*waXwc2 + vc1*wc2;
	//float A13=va1*wa3 + vaXvc1*waXwc3 + vc1*wc3;
	float A21=va2*wa1 + vaXvc2*waXwc1 + vc2*wc1;
	//float A22=va2*wa2 + vaXvc2*waXwc2 + vc2*wc2;
	//float A23=va2*wa3 + vaXvc2*waXwc3 + vc2*wc3;
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

	//some of these are unused
	float B11=vb1*wb1 + vbXvc1*wbXwc1 + vc1*wc1;
	//float B12=vb1*wb2 + vbXvc1*wbXwc2 + vc1*wc2;
	//float B13=vb1*wb3 + vbXvc1*wbXwc3 + vc1*wc3;
	float B21=vb2*wb1 + vbXvc2*wbXwc1 + vc2*wc1;
	//float B22=vb2*wb2 + vbXvc2*wbXwc2 + vc2*wc2;
	//float B23=vb2*wb3 + vbXvc2*wbXwc3 + vc2*wc3;
	float B31=vb3*wb1 + vbXvc3*wbXwc1 + vc3*wc1;
	float B32=vb3*wb2 + vbXvc3*wbXwc2 + vc3*wc2;
	float B33=vb3*wb3 + vbXvc3*wbXwc3 + vc3*wc3;
	
	//use weights based on magnitude
	//weighted triad
	float weightA=1.0/(variance+IMAGE_VARIANCE/new_s1.mag);
	float weightB=1.0/(variance+IMAGE_VARIANCE/new_s2.mag);

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
void set_mask(int x,int y, int id, float score, float variance){
	if (score>0) {
		//has this pixel already been assigned to a different star?
		unsigned char id2=img_mask[(x+IMG_X/2)+(y+IMG_Y/2)*IMG_X];
		if (id2!=255){
			float sigma_sq2=newstars[id2].sigma_sq+variance;
			float maxdist_sq2=-sigma_sq2*(log(sigma_sq2)+MATCH_VALUE);
			float px2=newstars[id2].px;
			float py2=newstars[id2].py;
			float a2=(x+.5-px2);
			float b2=(y+.5-py2);
			float score2 = (maxdist_sq2-(a2*a2+b2*b2))/(2*sigma_sq2);
			if (score>score2){
				img_mask[(x+IMG_X/2)+(y+IMG_Y/2)*IMG_X]=id;
			}
		} else {
			img_mask[(x+IMG_X/2)+(y+IMG_Y/2)*IMG_X]=id;
		}
	}
}

void add_score(struct constellation db_const,unsigned char newid1,unsigned char newid2){
	struct constellation_score cs;
	cs.oldid1=db_const.s1;
	cs.oldid2=db_const.s2;
	cs.newid1=newid1;
	cs.newid2=newid2;
	
	int i;
	cs.id_map=(int *)malloc(sizeof(int)*numnewstars);
	for (i=0;i<numnewstars;i++) cs.id_map[i]=-1;
	cs.scores=(float *)malloc(sizeof(float)*numnewstars);
	for (i=0;i<numnewstars;i++) cs.scores[i]=0.0;
	
	cs.totalscore=log(EXPECTED_FALSE_STARS/(IMG_X*IMG_Y))*(2*numnewstars);
	int smax=db_const.staridx+db_const.numstars;
	for(int i=db_const.staridx;i<smax;i++){
		int o=db_starids[i];
		float x=starptr[o].x*R11+starptr[o].y*R12+starptr[o].z*R13;
		float y=starptr[o].x*R21+starptr[o].y*R22+starptr[o].z*R23;
		float z=starptr[o].x*R31+starptr[o].y*R32+starptr[o].z*R33;
		float px=y/(x*PIXX_TANGENT);
		float py=z/(x*PIXY_TANGENT);
		unsigned char n=255;
		if (fabs(2*px)<IMG_X && fabs(2*py)<IMG_Y)
			n=img_mask[(int)(px+IMG_X/2.0f)+(int)(py+IMG_Y/2.0f)*IMG_X];
		if (n!=255) {
			float sigma_sq=newstars[n].sigma_sq+POS_VARIANCE;
			float maxdist_sq=-sigma_sq*(log(sigma_sq)+MATCH_VALUE);
			float a=(px-newstars[n].px);
			float b=(py-newstars[n].py);
			float score = (maxdist_sq-(a*a+b*b))/(2*sigma_sq);
			//only match the closest star
			if (score>cs.scores[n]){
				cs.id_map[n]=o;
				cs.scores[n]=score;
			}
		}
	}
	for(int n=0;n<numnewstars;n++) cs.totalscore+=cs.scores[n];
	c_scores=(struct constellation_score*)realloc(c_scores,sizeof(struct constellation_score)*(c_scores_size+1));
	c_scores[c_scores_size].totalscore=cs.totalscore;
	c_scores[c_scores_size].oldid1=cs.oldid1;
	c_scores[c_scores_size].oldid2=cs.oldid2;
	c_scores[c_scores_size].newid1=cs.newid1;
	c_scores[c_scores_size].newid2=cs.newid2;
	c_scores[c_scores_size].id_map=cs.id_map;
	c_scores[c_scores_size].scores=cs.scores;
	//printf("totalscore=%f\n",cs.totalscore);
	c_scores_size++;
}

float search() {
	int i;
	free(winner_id_map);
	winner_id_map=(int *)malloc(sizeof(int)*addednewstars);
	for (i=0;i<addednewstars;i++) winner_id_map[i]=-1;
	free(winner_scores);
	winner_scores=(float *)malloc(sizeof(float)*addednewstars);
	for (i=0;i<addednewstars;i++) winner_scores[i]=0.0;
	
	//Do we have enough stars?
	if (NUMSTARS<2||numnewstars<2) return 0.0;

	//allocate space for image mask 
	img_mask = (unsigned char*)malloc(IMG_X*IMG_Y);
	memset(img_mask, 255, IMG_X*IMG_Y);
	
	//generate image mask
	for (int id=0;id<numnewstars;id++){
		//assume the dimmest possible star since we dont know the brightness of the other image
		float sigma_sq=newstars[id].sigma_sq+POS_VARIANCE;
		float maxdist_sq=-sigma_sq*(log(sigma_sq)+MATCH_VALUE);
		float maxdist=sqrt(maxdist_sq);
		int xmin=fmaxf(-IMG_X/2.0f,(newstars[id].px-maxdist));
		int xmax=fminf(IMG_X/2.0f,newstars[id].px+maxdist+1);
		int ymin=fmaxf(-IMG_Y/2.0f,(newstars[id].py-maxdist));
		int ymax=fminf(IMG_Y/2.0f,(newstars[id].py+maxdist+1));
		for(int x=xmin;x<xmax;x++) for (int y=ymin;y<ymax;y++) {
			float a=(x+.5-newstars[id].px);
			float b=(y+.5-newstars[id].py);
			set_mask(x,y,id,(maxdist_sq-(a*a+b*b))/(2*sigma_sq),POS_VARIANCE);
		}
	}
	//use weighted triad
	c_scores_clear();
	for (int i=1;i<numnewstars&&i<MAX_FALSE_STARS+2;i++) {
		for (int j=0;j<i;j++) {
			float p=(3600*180.0/PI)*asin(dist3(newstars[i].x,newstars[j].x,newstars[i].y,newstars[j].y,newstars[i].z,newstars[j].z));
			int mmi=(int)(p/ARC_ERR)%mapsize;
			if (mmi<0) mmi+=mapsize;
			for (int constidx=map[mmi];constidx!=-1;) {
				if (fabs(constptr[constidx].p-p)<ARC_ERR){
					weighted_triad(starptr[constptr[constidx].s1],starptr[constptr[constidx].s2],newstars[j],newstars[i],POS_VARIANCE);
					add_score(constptr[constidx],j,i);
					weighted_triad(starptr[constptr[constidx].s1],starptr[constptr[constidx].s2],newstars[i],newstars[j],POS_VARIANCE);
					add_score(constptr[constidx],i,j);
				}
				if (mmi==constptr[constidx].mmi_l) constidx=constptr[constidx].last_l;
				else if (mmi==constptr[constidx].mmi_h) constidx=constptr[constidx].last_h;
				else constidx=constptr[constidx].last_m;
			}
		}
	}
	qsort(c_scores, c_scores_size, sizeof(struct constellation_score), compare_totalscore);
	//add up probabilities of all matches, excluding those which
	//are equivalent to the best match (S1==s1,S2=s2)
	float p_match=1.0;
	if (c_scores_size>0) {
		for(int n=0;n<numnewstars;n++) {
			int o=c_scores[0].id_map[n];
			if (o!=-1){
				winner_scores[newstars[n].starnum]=c_scores[0].scores[n];
				winner_id_map[newstars[n].starnum]=starptr[o].id;
			}
		}
		
		float bestscore=c_scores[0].totalscore;
		int oldid1=c_scores[0].oldid1;
		int oldid2=c_scores[0].oldid2;
		unsigned char newid1=c_scores[0].newid1;
		unsigned char newid2=c_scores[0].newid2;
		//set attitude matrix to best match
		weighted_triad(starptr[oldid1],starptr[oldid2],newstars[newid1],newstars[newid2],POS_VARIANCE);
		for(int i=1;i<c_scores_size;i++) {
			if (c_scores[i].id_map[newid1]!=oldid1&&c_scores[i].id_map[newid2]!=oldid2){
				p_match+=exp(c_scores[i].totalscore-bestscore);
			}
		}
		p_match=1.0/p_match;

	} else {
		p_match=0.0;
	}
	free(img_mask);
	return p_match;
}

void flip() {
		addedoldstars=addednewstars;
		numoldstars=numnewstars;
		memcpy(oldstars,newstars,sizeof(newstars));
		addednewstars=0;
		numnewstars=0;
		memset(newstars,0,sizeof(newstars));
		c_scores_clear();
}

/***
 * Prepares the data from the data file for final use.
 * This function should not run expensive preprocessing steps.
 * It should just put the data in a structure to be
 * usable by the star_id function.
 * 
 * @param data The data loaded from the _data_ file.
 * @param data_size The size of the data loaded from the file.
 * @return The prepared data.
 */
void* prepare(void* data, long data_size)
{
	load_db();
	return data;
}

/***
 * Determines the HIP numbers of stars in a scene.
 * 
 * @param spikes The data of the spikes in the scene. x, y, and magnitude interleaved.
 * @param result The output array that should be filled with the HIP numbers.
 * @param length The number of spikes in the scene.
 * @param data The data prepared by the prepare function.
 */
void star_id(double spikes[], int result[], size_t length, void* data)
{
	flip();
	for(size_t i = 0; i < length; i++)
	{
		add_star(spikes[3*i]-IMG_X/2.0,-(spikes[3*i+1]-IMG_Y/2.0),spikes[3*i+2]);
		result[i] = -1;
	}
	float p_match = search();
	if (p_match>0.66) {
		for(size_t i = 0; i < length; i++) {
			result[i] = winner_id_map[i];
		}
	}
}
