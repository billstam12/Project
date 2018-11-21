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
	int type = 1; //0 random, 1 kmeans++

	centroids = init_centroids(k, data, no_of_samples, no_of_dimensions, type);
	
	/*
	for(i = 0; i < k; i++){
		printf("ID= %ld\n", centroids[i]->id);		
		printf("COUNT= %d\n",centroids[i]->count);
		print_coordinates_cent(centroids[i],  no_of_dimensions);
	}
	*/

	/* Compute the distance of each point from a centroid
	and assign each point to a centroid and each centroid to 
	the points it has */

	type = 1; //0 normal, 1 PAM
	int assignment = 0; //0 lloyds, 1 lsh, 2 hyperplane
	kmeans(centroids,data,no_of_samples,no_of_dimensions,k, assignment, metric, type); 

	//Evaluate centroids
	for(i = 0 ; i < no_of_samples; i++){
		int a_id = data[i]->centroid_id;
		int b_id = data[i]->centroid2_id;

		long double a = 0;
		long double b = 0;
		for(j = 0; j < centroids[a_id]->count; j++){
			a += euclidean_distance(data[i]->coordinates, centroids[a_id]->assigned_points[j]->coordinates, no_of_dimensions);
		}
		for(j = 0; j < centroids[b_id]->count; j++){
			b += euclidean_distance(data[i]->coordinates, centroids[b_id]->assigned_points[j]->coordinates, no_of_dimensions);
		}
		double max;
		if(a > b){
			max = a;
		}
		else{ max = b;}
		data[i]->silhouette= ((b-a)/max);
	}
}