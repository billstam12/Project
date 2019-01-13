#include <stdio.h>
#include <stdlib.h>
#include "functions.h"


int main(int argc, char** argv){
	FILE * f, * c;
	int i,j;
	char output[50];
	int p = 5; // Nearest Neighbors
	for(i = 0; i < argc;i++){
		if(!strcmp(argv[i],"-d")){
			f = fopen(argv[++i],"r");
		}
		if(!strcmp(argv[i],"-p")){
			p = atoi(argv[++i]);
		}
		if(!strcmp(argv[i],"-o")){
			strcpy(output,argv[++i]);
		}
	}

	users_and_tweets data;
	tweet * tweets;
	user * users, * users2, *users3, *users4;
	char** dict; //coins dictionary
	cluster * clusters;
	int no_of_samples, no_of_users, no_of_clusters;
	int no_of_coins;
	data = parse_data(f); //Data contains both users and tweets
	//Get all data from the struct
	no_of_samples = data->no_of_tweets;
	no_of_users = data->no_of_users;
	no_of_clusters = data->no_of_clusters;
	no_of_coins = data->no_of_coins;
	tweets = data->tweets;
	users = data->users;
	users2 = copy_users(users, no_of_users);
	users3 = copy_users(users, no_of_users);
	users4 = copy_users(users, no_of_users);
	clusters = data->clusters;
	dict = data->dict;
	/* QUESTION A */
	/*Find the nearest neighbors of each user with lsh */
	cosine_lsh_implementation(users, no_of_users, p , 0, dict ,output); //Done

	/* Find the nearest neighbors of each user using clustering */
	// Clustering method used will be the best from the previous exercise 
	// Looyds basic assignment with basic update, and k-means++ initialization
	int k = 4; 
	int type = 0;
	cosine_cluster_implementation(users2, no_of_users, 100, k, p, type, dict, output);
	
	/* QUESTION B */
	clustering_lsh_implementation(users3, clusters, no_of_users, no_of_clusters, p, dict, output);
	clustering_cluster_implementation(users4, clusters, no_of_users, no_of_clusters, p, dict, output);

	/* Free all data */
	free_data(tweets, no_of_samples, users, no_of_users, clusters, dict, no_of_coins);
	free_users(users2, no_of_users);
	free_users(users3, no_of_users);
	free_users(users4, no_of_users);
	return 0;
}