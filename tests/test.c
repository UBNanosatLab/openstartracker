/****************************************************************************
 * Copyright (c) 2016 Joerg H. Mueller <nexyon@gmail.com>                   *
 *                                                                          *
 * Permission is hereby granted, free of charge, to any person obtaining a  *
 * copy of this software and associated documentation files (the            *
 * "Software"), to deal in the Software without restriction, including      *
 * without limitation the rights to use, copy, modify, merge, publish,      *
 * distribute, distribute with modifications, sublicense, and/or sell       *
 * copies of the Software, and to permit persons to whom the Software is    *
 * furnished to do so, subject to the following conditions:                 *
 *                                                                          *
 * The above copyright notice and this permission notice shall be included  *
 * in all copies or substantial portions of the Software.                   *
 *                                                                          *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS  *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF               *
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.   *
 * IN NO EVENT SHALL THE ABOVE COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,   *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR    *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR    *
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.                               *
 *                                                                          *
 * Except as contained in this notice, the name(s) of the above copyright   *
 * holders shall not be used in advertising or otherwise to promote the     *
 * sale, use or other dealings in this Software without prior written       *
 * authorization.                                                           *
 ****************************************************************************/

#include "../src/beast.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_STARS 1000
#define VECTOR_LENGTH 3

#define EDISON_SPEED_FACTOR 30

star_db *S_DB;
star_query *S_RESULTS;
constellation_db *C_DB;

/***
 * Determines the HIP numbers of stars in a scene.
 * 
 * @param spikes The data of the spikes in the scene. x, y, and magnitude interleaved.
 * @param result The output array that should be filled with the HIP numbers.
 * @param length The number of spikes in the scene.
 */
int IDNUM;
void star_id(double spikes[], int result[], size_t length)
{	
	star_db* s=new star_db;
	for(size_t i = 0; i < length; i++)
	{
		s->add_star(spikes[3*i]-IMG_X/2.0,-(spikes[3*i+1]-IMG_Y/2.0),BASE_FLUX*powf(10.0,-spikes[3*i+2]/2.5),-1);
		result[i] = -1;
	}
	
	constellation_db * img=new constellation_db(s,MAX_FALSE_STARS+2,1);
	db_match* lis = new db_match(C_DB,img);
	if (lis->p_match>0.99) {
		float x=lis->winner->R11;
		float y=lis->winner->R12;
		float z=lis->winner->R13;
		
		//Tests relative matching, and fills in missing stars
		//and fill in missing stars
		S_RESULTS->kdsearch(x,y,z,MAXFOV/2,THRESH_FACTOR*IMAGE_VARIANCE);
		C_DB->results->kdsearch(x,y,z,MAXFOV/2,THRESH_FACTOR*IMAGE_VARIANCE);
		constellation_db* fov_db = new constellation_db(S_RESULTS->from_kdresults(),C_DB->results->kdresults_size,1);
		C_DB->results->clear_kdresults();
		S_RESULTS->clear_kdresults();
		
		db_match* fov_match = new db_match(fov_db,img);
		
		for(size_t i = 0; i < length; i++) result[i] = fov_match->map[i];
		delete fov_db;
		delete fov_match;
	}
	delete lis;
	delete img;
}


int main(int argc, char* argv[])
{
	if(argc < 3) {
		printf("./test input.csv calibration.txt\n");
		return 0;
	}
	
	load_config(argv[2]);
	S_DB=new star_db;
	S_DB->load_catalog("hip_main.dat",1991.25);
	
	S_RESULTS=new star_query(S_DB);
	S_RESULTS->kdmask_filter_catalog();
	S_RESULTS->kdmask_uniform_density(REQUIRED_STARS);
	
	C_DB=new constellation_db(S_RESULTS->from_kdmask(),2+DB_REDUNDANCY,0);
	S_RESULTS->reset_kdmask();
	
	// read scenes and run star identification
	FILE* file = fopen(argv[1], "r");

	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	double* data = (double*)malloc(sizeof(double) * VECTOR_LENGTH * MAX_STARS);
	int* results = (int*)malloc(sizeof(int) * MAX_STARS);

	clock_t time_sum = 0;
	int run_times = 0;
	for(IDNUM=0;(read = getline(&line, &len, file)) != -1;IDNUM++) {
		size_t i = 0;
		for(char* result = strtok(line, ","); result != NULL; i++) {
			data[i] = atof(result);
			result = strtok(NULL, ",");
		}
		int length = i / VECTOR_LENGTH;
		clock_t start = clock();
		star_id(data, results, length);
		clock_t end = clock();
		time_sum += end - start;
		for(i = 0; i < length; i++) printf("%d%c", results[i], i == length - 1 ? '\n' : ',');
		run_times++;
		fprintf(stderr,"Time on edison: %f\n", ((float)time_sum*EDISON_SPEED_FACTOR) / (CLOCKS_PER_SEC*run_times));
	}



	free(line);
	free(results);
	free(data);
	fclose(file);

	return 0;
}

