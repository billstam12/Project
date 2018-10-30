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
	int window = 400;
	int L = 4;
	int no_of_functions = 4;


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
		if(!strcmp(argv[i], "-L")){
			L = atoi(argv[++i]);
		}
	}
	int end=0;
	//Get dataset, its rows and its columns 
	while(end==0){	
		int no_samples = 0;
		int no_dimensions = 0;
		int metric = 0; // 0 eucl, 1 cosine
		hashtable * hts;
		data = parse_data(f,&no_samples,&no_dimensions,&metric);
		int table_size = no_samples/4;
		int no_queries = 0;
		double radius = 0;
		query_data = parse_query(q, &no_queries, no_dimensions, &radius);

		/*  r, t and random_vectors are 2d, 2d, and 3d vectors each
			one dimension is for each of the hashtables and the other
			is for the data  the arrays contain
		*/
		
		double ** t; // t array
		t = (double**)malloc(L*sizeof(double*));
		int ** r; // r array
		r = (int**)malloc(L*sizeof(int*));
		double *** random_vectors; // random hyperplane vectors
		random_vectors = (double***)malloc(L*sizeof(double**));
		hts = (hashtable*)malloc(L * sizeof(hashtable));
	
		if(metric == 0){
			for(i = 0; i < L; i ++){
				hashtable_init(&hts[i], table_size);
			}
			for(i=0;i<no_samples;i++){
				data[i]->g_functions = (long long int**)malloc(sizeof(long long int)* L);
			}	
			for(i = 0; i < L; i++){
				random_vectors[i] = init(&hts[i], data ,  no_samples,  no_dimensions, window,  no_of_functions, table_size, L,  i, metric, &t[i], &r[i]);
			}

			/* Now we will give the queries to our hashtables */	
			euclidean_lsh_query(data,  query_data,  hts,  no_queries, no_samples, no_dimensions,  L,  radius,  window,  no_of_functions,  table_size,  t,  r,  random_vectors);
		}
		else{
			table_size = pow(2,no_of_functions);
			for(i = 0; i < L; i ++){
				hashtable_init(&hts[i], table_size);
			}
			for(i = 0; i < L; i++){
				random_vectors[i] = init(&hts[i], data ,  no_samples,  no_dimensions, window,  no_of_functions, table_size, L, i, metric, &t[i], &r[i]);
			}
			cosine_query(data,  query_data,  hts,  no_queries, no_samples, no_dimensions,  L,  radius,  window,  no_of_functions,  table_size,  t,  r,  random_vectors);
		
		}
		/* Free Space */

		//Free t, r and random vectors

		for(i = 0; i < no_of_functions; i++){
			free(t[i]);
			free(r[i]);
			free(random_vectors[i]);
		}
		free(t);
		free(r);
		free(random_vectors);
		
		//free hashtables
		for(i = 0; i < L; i++){
			hashtable_free(&hts[i]);
		}

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