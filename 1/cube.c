#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include <time.h>
#include "functions.h"


int main(int argc, char** argv){
	int seed = time(NULL);
	srand(seed);
	int i, j;
	int probes;
	int max = 10;
	int no_of_functions = 0;

	point * data;
	point * query_data;
	FILE *f;
	FILE *q;
	
	for(i = 0; i < argc;i++){
		if(!strcmp(argv[i],"-d")){
			f = fopen(argv[++i],"r");
		}
		if(!strcmp(argv[i],"-q")){
			q = fopen(argv[++i],"r");
		}
		if(!strcmp(argv[i], "-k")){
			no_of_functions = atoi(argv[++i]);
		}
		if(!strcmp(argv[i], "-M")){
			max = atoi(argv[++i]);
		}
		if(!strcmp(argv[i], "-probes")){
			probes = atoi(argv[++i]);
		}
	}
	int end = 0;
	while(end==0){	
		//Get dataset, its rows and its columns 		
		int no_samples = 0;
		int no_dimensions = 0;
		int metric;
		data = parse_data(f,&no_samples,&no_dimensions,&metric);
		int no_queries = 0;
		double radius = 0;
		query_data = parse_query(q, &no_queries, no_dimensions, &radius);	
		if(no_of_functions == 0){
			no_of_functions = log10(no_samples);
		}

		int table_size = pow(2,no_of_functions);

		hashtable hypercube = (hashtable)malloc(sizeof(struct hashtable));
		hashtable_init(&hypercube, table_size);
		double ** hyperplanes;
		hyperplanes = (double**)malloc(sizeof(double**));
		for(i = 0; i < no_of_functions; i++){
			hyperplanes[i] = create_random_vector(no_dimensions);
		}
		long long int index;
		for(i = 0; i < no_samples; i++){
			index = cosine_hash(&data[i], no_dimensions, no_of_functions, table_size, hyperplanes);
			hashtable_insert(&hypercube, data[i], index);
		}
		hypercube_query(data, query_data, hypercube, no_queries, no_samples, no_dimensions, radius, no_of_functions, table_size, probes, max, hyperplanes);
		
		/* Free Space */

		//Free random vectors
		free(hyperplanes);
		
		//free hashtables
		hashtable_free(&hypercube);

		//free data points
		for(i = 0; i<no_samples; i++){
			free(data[i]->coordinates);
			free(data[i]);
		}
		for(i = 0; i<no_queries; i++){
			free(query_data[i]->coordinates);
			free(query_data[i]);
		}
		free(data);
		free(query_data);
		
		
		// Probe user to continue
		char probe[1];
		printf("Would you like to continue? [y/n]\n");
		scanf("%s",probe);
		if(strcmp(probe,"y")==0){
			char new_query[50];
			char new_data[50];
			printf("New data? [y/n]\n");
			scanf("%s",probe);
			if(strcmp(probe,"y")==0){
				printf("Give a new data file!\n");
				scanf("%s", new_data);
				f = fopen(new_data,"r");
			}
			printf("New Query? [y/n]\n");
			scanf("%s",probe);
			if(strcmp(probe,"y")==0){
				printf("Give a new query file!\n");
				scanf("%s", new_query);
				q = fopen(new_query,"r");
			}
		}
		else{
			end = 1;
			printf("Finishing...\n");
		}
	}
	return 0;
}