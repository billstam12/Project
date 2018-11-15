#include <stdio.h>
#include <stdlib.h>
#include "functions.h"


int main(int argc, char** argv){
	
	int metric  = 0; //0 euclidean, 1 cosine
	int k; // # clusters
	int no_of_functions = 4;
	int L = 5; // # hashtables

	point * data;
	FILE *f;
	FILE *c;
	//char out[50];
	int i,j;
	for(i = 0; i < argc;i++){
		if(!strcmp(argv[i],"-i")){
			f = fopen(argv[++i],"r");
		}
		if(!strcmp(argv[i],"-c")){
			c = fopen(argv[++i],"r");
		}
		if(!strcmp(argv[i], "-d")){
			metric = atoi(argv[++i]);
		}
		/*if(!strcmp(argv[i], "-o")){
			strcpy(out,argv[++i]);
		}*/
	}

	// Parse data files
	int no_of_samples;
	int no_of_dimensions;
	data = parse_data(f, &no_of_samples, &no_of_dimensions);
	parse_conf(c, &k, &no_of_functions, &L);
	//printf("Clusters:%d, functions:%d, hashtables:%d\n",k, no_of_functions, L);

	//Create centroids
	centroid * centroids;
	int type = 0; //0 random, 1 kmeans++

	centroids = init_centroids(k, data, no_of_samples, type);

	/* Compute the distance of each point from a centroid
	and assign each point to a centroid and each centroid to 
	the points it has */
	kmeans(centroids,data,no_of_samples,no_of_dimensions,k); //Loydds kmeans
	
}