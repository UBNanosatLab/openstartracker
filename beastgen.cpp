#include <iostream>	 // cout
#include <string>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>	/* For O_RDWR */
#include <sys/mman.h>
#include "configuration.h"
#include <stdint.h>

struct constellation {
	float p;
	int32_t s1;
	int32_t s2;
	int32_t numstars;
	int32_t staridx;
	int32_t last;
};

int *map;
struct constellation *constptr;
int *star_ids;

void  add_entry(int mapidx,int curr_const) {
	int *staridx;
	for (staridx=&map[mapidx];*staridx!=-1;staridx=&(constptr[*staridx].last));
	if (staridx!=&(constptr[curr_const].last)) *staridx=curr_const;
}

int main (int argc, char** argv) {
	std::cout.precision(12);
	//load config
	configuration::data config1;
	std::ifstream cfgfile1("calibration/dbsize.txt");
	cfgfile1 >> config1;
	cfgfile1.close();
	
	configuration::data config2;
	std::ifstream cfgfile2("calibration/calibration.txt");
	cfgfile2 >> config2;
	cfgfile2.close();
	
	std::ifstream constline("calibration/constellations.txt");
	
	int PARAM=atoi(config1["PARAM"].c_str());
	int NUMCONST=atoi(config1["NUMCONST"].c_str());
	int STARTABLE=atoi(config1["STARTABLE"].c_str());
	float ARC_ERR=atof(config2["ARC_ERR"].c_str());
	int mapsize=PARAM;
	int s_offset=mapsize*sizeof(int)+ NUMCONST*sizeof(struct constellation);
	size_t dbsize = s_offset + STARTABLE*sizeof(int);
	
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

	map = (int*)mmap(0, dbsize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (map == MAP_FAILED)
	{
		close(fd);
		perror("Error mmapping the file");
		exit(EXIT_FAILURE);
	}
	constptr=(struct constellation*)(&map[mapsize]);
	star_ids=(int*)(&map[s_offset/sizeof(int)]);
	memset(map, -1, dbsize);
	int mmi_l,mmi_m,mmi_h;
	int curr_star=0;
	for (int curr_const=0;curr_const<NUMCONST;curr_const++){
		constline>>constptr[curr_const].p;
		constline>>constptr[curr_const].s1;
		constline>>constptr[curr_const].s2;
		constline>>constptr[curr_const].numstars;
		constptr[curr_const].staridx=curr_star;
		for (int i=0;i<constptr[curr_const].numstars;i++) {
			constline>>star_ids[curr_star];
			curr_star++;
		}
		//add entry to constellation table
		mmi_l=(int)(constptr[curr_const].p/ARC_ERR-1)%mapsize;
		mmi_m=(int)(constptr[curr_const].p/ARC_ERR)%mapsize;
		mmi_h=(int)(constptr[curr_const].p/ARC_ERR+1)%mapsize;
		if (mmi_l<0) mmi_l+=mapsize;
		if (mmi_m<0) mmi_m+=mapsize;
		if (mmi_h<0) mmi_h+=mapsize;
		
		add_entry(mmi_l,curr_const);
		add_entry(mmi_m,curr_const);
		add_entry(mmi_h,curr_const);
	}
	
	// Write it now to disk
	if (msync(map, dbsize, MS_SYNC) == -1) {
		perror("Could not sync the file to disk");
	}
	
	// Don't forget to free the mmapped memory
	if (munmap(map, dbsize) == -1) {
		close(fd);
		perror("Error un-mmapping the file");
		exit(EXIT_FAILURE);
	}
	// Un-mmaping doesn't close the file, so we still need to do that.
	close(fd);
	constline.close();
}

