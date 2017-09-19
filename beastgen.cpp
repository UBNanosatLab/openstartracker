#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>	/* For O_RDWR */
#include <sys/mman.h>
#include <stdint.h>

#include "beast.h"
//TODO: merge into "write" function 
int main (int argc, char** argv) {
	//load config
	FILE *stream;
	char *line = NULL;
	size_t len = 0;
	stream = fopen("calibration/calibration.txt", "r");
	if (stream == NULL) exit(EXIT_FAILURE);
	while (getline(&line, &len, stream) != -1) {
		//Don't listen to valgrind. This is fine. Everything is fine.
		putenv(strcpy((char *)malloc(sizeof(char) * len),line));
	}
	fclose(stream);
	
	stream = fopen("calibration/dbsize.txt", "r");
	if (stream == NULL) exit(EXIT_FAILURE);
	while (getline(&line, &len, stream) != -1) {
		//Don't listen to valgrind. This is fine. Everything is fine.
		putenv(strcpy((char *)malloc(sizeof(char) * len),line));
	}
	fclose(stream);
	
	//constline
	stream=fopen("calibration/constellations.txt", "r");
	
	NUMCONST=atoi(getenv("NUMCONST"));
	STARTABLE=atoi(getenv("STARTABLE"));
	size_t dbsize = NUMCONST*sizeof(struct constellation) + STARTABLE*sizeof(int);
	
	/* Open a file for writing.
	 *  - Creating the file if it doesn't exist.
	 *  - Truncating it to 0 size if it already exists. (not really needed)
	 *
	 * Note: "O_WRONLY" mode is not sufficient when mmaping.
	 */

	const char *filepath = "beastdb.bin";
	int fd = open(filepath, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
	if (fd == -1) {
		perror("Error opening file for writing");
		exit(EXIT_FAILURE);
	}

	// Stretch the file size to the size of the (mmapped) array of char


	if (lseek(fd, dbsize-1, SEEK_SET) == -1) {
		close(fd);
		perror("Error calling lseek() to 'stretch' the file");
		exit(EXIT_FAILURE);
	}
	
	if (write(fd, "", 1) == -1) {
		close(fd);
		perror("Error writing last byte of the file");
		exit(EXIT_FAILURE);
	}
	// Now the file is ready to be mmapped.

	struct constellation *db_map = (struct constellation*)mmap(0, dbsize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (db_map == MAP_FAILED)
	{
		close(fd);
		perror("Error mmapping the file");
		exit(EXIT_FAILURE);
	}
	int *db_stars=(int*)(&db_map[NUMCONST]);
	memset(db_map, -1, dbsize);
	int curr_star=0;
	for (uint32_t idx=0;idx<NUMCONST;idx++){
		if (getline(&line, &len, stream) != -1) {
			db_map[idx].p=atof(strtok(line," "));
			db_map[idx].s1=atoi(strtok(NULL," "));
			db_map[idx].s2=atoi(strtok(NULL," "));
			db_map[idx].numstars=atoi(strtok(NULL," "));
			db_map[idx].firststar=curr_star;
			db_map[idx].idx=idx;
			for (int i=0;i<db_map[idx].numstars;i++) {
				db_stars[curr_star]=atoi(strtok(NULL," "));
				curr_star++;
			}
		}
	}
	
	// Write it now to disk
	if (msync(db_map, dbsize, MS_SYNC) == -1) {
		perror("Could not sync the file to disk");
	}
	
	// Don't forget to free the mmapped memory
	if (munmap(db_map, dbsize) == -1) {
		close(fd);
		perror("Error un-mmapping the file");
		exit(EXIT_FAILURE);
	}
	// Un-mmaping doesn't close the file, so we still need to do that.
	close(fd);
	fclose(stream);
}

