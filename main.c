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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_STARS 200
#define VECTOR_LENGTH 3

#include "identification.c"

int main(int argc, char* argv[])
{
	char measure_time = 0;
	
	int repetitions = 1;
	int run_times = 1;

	if(argc > 1)
	{
		measure_time = 1;
		repetitions = atoi(argv[1]);
		run_times = 3;

		if(argc > 2)
			run_times = atoi(argv[2]);
	}

	// read data file

	FILE* file = fopen("data", "rb");

	long user_size = 0;
	void* user_data = NULL;
	void* prepared_data = NULL;

	if(file)
	{
		fseek(file, 0, SEEK_END);
		user_size = ftell(file);
		rewind(file);

		user_data = malloc(user_size);
		user_size = fread(user_data, user_size, 1, file);

		prepared_data = prepare(user_data, user_size);

		fclose(file);
	}

	// read scenes and run star identification

	file = fopen("input.csv", "r");

	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	double* data = (double*)malloc(sizeof(double) * VECTOR_LENGTH * MAX_STARS);
	int* results = (int*)malloc(sizeof(int) * MAX_STARS);

	clock_t best = 0;
	for(int run_time = 0; run_time < run_times; run_time++)
	{
		clock_t time_sum = 0;

		for(int repetition = 0; repetition < repetitions; repetition++)
		{

			while((read = getline(&line, &len, file)) != -1)
			{
				size_t i = 0;
		
				for(char* result = strtok(line, ","); result != NULL; i++)
				{
					data[i] = atof(result);
					result = strtok(NULL, ",");
				}
		
				int length = i / VECTOR_LENGTH;

				clock_t start = clock();
				star_id(data, results, length, prepared_data);

				clock_t end = clock();

				time_sum += end - start;

				if(!run_time && !repetition)
				{
					for(i = 0; i < length; i++)
					{
						printf("%d%c", results[i], i == length - 1 ? '\n' : ',');
					}
				}
			}

			rewind(file);
		}

		time_sum /= repetitions;

		if(best == 0)
			best = time_sum;
		else if(best > time_sum)
			best = time_sum;
	}

	if(measure_time)
		printf("Time elapsed: %f\n", ((float)best) / CLOCKS_PER_SEC);

	free(user_data);
	free(line);
	free(results);
	free(data);
	fclose(file);

	return 0;
}

