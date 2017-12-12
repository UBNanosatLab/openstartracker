#ifndef CALIBRATION_H
#define CALIBRATION_H

#include <stdio.h>
#include <stdlib.h>//EXIT_FAILURE
#include <string.h>
#include <math.h>
#include <set>

#define PI		   3.14159265358979323846  /* pi */
#define TWOPI		6.28318530717958647693

int IMG_X,IMG_Y,MAX_FALSE_STARS,REQUIRED_STARS;
float DEG_X,DEG_Y,PIXX_TANGENT,PIXY_TANGENT;
float PIXSCALE,POS_ERR_SIGMA,POS_VARIANCE;
float IMAGE_VARIANCE,EXPECTED_FALSE_STARS,MATCH_VALUE;
float PHOTONS,BRIGHT_THRESH;
float TODO_MAXFOV_D2;
float TODO_MINFOV_D2;

std::set<char*> ENV_VARS; //keep track of these so that valgrind is happy

void load_config() {
	
	//move calibration.txt to constellation_db
	//move dbstats.txt to beastdb
	/* load config */
	
	FILE *stream = fopen("calibration.txt", "r");
	if (stream == NULL) exit(EXIT_FAILURE);

	ssize_t read;
	char *line = NULL;
	size_t len = 0;
	while ((read = getline(&line, &len, stream)) != -1) {
		char * s=(char *)malloc(sizeof(char) * len);
		putenv(strcpy(s,line));
		ENV_VARS.insert(s);
	}
	free(line);
	fclose(stream);

	IMG_X=atoi(getenv("IMG_X"));
	IMG_Y=atoi(getenv("IMG_Y"));
	TODO_MAXFOV_D2=sqrt(IMG_X*IMG_X+IMG_Y*IMG_Y)/2;
	TODO_MINFOV_D2=IMG_Y/2;
	DEG_X=atof(getenv("DEG_X"));
	DEG_Y=atof(getenv("DEG_Y"));
	POS_ERR_SIGMA=atof(getenv("POS_ERR_SIGMA"));
	PIXSCALE=atof(getenv("PIXSCALE"));
	POS_VARIANCE=atof(getenv("POS_VARIANCE"));/* sigma_r^2 */
	IMAGE_VARIANCE=atof(getenv("IMAGE_VARIANCE"));/* lambda */
	BRIGHT_THRESH=atof(getenv("BRIGHT_THRESH"));
	EXPECTED_FALSE_STARS=atof(getenv("EXPECTED_FALSE_STARS"));/* pfalse */
	MAX_FALSE_STARS=atoi(getenv("MAX_FALSE_STARS"));/* >10 is slow */
	REQUIRED_STARS=atoi(getenv("REQUIRED_STARS"));/* >10 is slow */
	PHOTONS=atoi(getenv("PHOTONS"));
	MATCH_VALUE=4*log(EXPECTED_FALSE_STARS/(IMG_X*IMG_Y))+log(2*PI);/* base */

	PIXX_TANGENT=2*tan(DEG_X*PI/(180*2))/IMG_X;
	PIXY_TANGENT=2*tan(DEG_Y*PI/(180*2))/IMG_Y;
}
#endif
