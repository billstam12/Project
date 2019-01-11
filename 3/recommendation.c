#include <stdio.h>
#include <stdlib.h>
#include "functions.h"


int main(int argc, char** argv){
	
	FILE * f, * c;
	int i,j;
	int p = 5; // Nearest Neighbors
	for(i = 0; i < argc;i++){
		if(!strcmp(argv[i],"-d")){
			f = fopen(argv[++i],"r");
		}
		if(!strcmp(argv[i],"-p")){
			p = atoi(argv[++i]);
		}
	}

	users_and_tweets data;
	tweet * tweets;
	user * users;
	centroid * centroids;
	cluster * clusters;
	int no_of_samples, no_of_users, no_of_clusters;
	
	data = parse_data(f); //Data contains both users and tweets
	//Get all data from the struct
	no_of_samples = data->no_of_tweets;
	no_of_users = data->no_of_users;
	no_of_clusters = data->no_of_clusters;
	tweets = data->tweets;
	users = data->users;
	clusters = data->clusters;
	
	/* QUESTION A */
	/*Find the nearest neighbors of each user with lsh */
	//cosine_lsh_implementation(users, no_of_users, p , 0); //Done

	/* Find the nearest neighbors of each user using clustering */
	// Clustering method used will be the best from the previous exercise 
	// Looyds basic assignment with basic update, and k-means++ initialization
	int k = 4; 
	int type = 0;
	//centroids = init_centroids(k, users, no_of_users, 100, type);
	//kmeans(centroids, users, no_of_users, 100, k, p); 

	/* QUESTION B */
	//clustering_lsh_implementation(users, clusters, no_of_users, no_of_clusters, p);
	clustering_cluster_implementation(users, clusters, no_of_users, no_of_clusters, p);
	for(i = 0; i < no_of_samples; i++){
	    free_list(tweets[i]->coin_list);
		for(j = 0; j < tweets[i]->no_of_words; j++){
			free(tweets[i]->words[j]);
		}
		free(tweets[i]->words);
		free(tweets[i]);
	}
	free(tweets);
	return 0;
}