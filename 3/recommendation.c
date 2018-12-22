#include <stdio.h>
#include <stdlib.h>
#include "functions.h"


int main(int argc, char** argv){
	
	FILE * f;
	int i,j;
	for(i = 0; i < argc;i++){
		if(!strcmp(argv[i],"-d")){
			f = fopen(argv[++i],"r");
		}
	}

	tweet * tweets;
	int no_of_samples;
	tweets = parse_data(f, &no_of_samples);
}