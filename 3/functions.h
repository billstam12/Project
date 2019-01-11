#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include <assert.h>
#include <time.h>
#include "hashtable.h"
#include <stddef.h>

#ifdef _WIN32
	#include <windows.h>
#else
	#include <unistd.h>
#endif


typedef struct unique{
	int count;
	int* array;
	int* tweets_per_user;
}  unique;

typedef struct list_of_coins{
	char** array;
	int no_of_elements;
} *list_of_coins;

typedef struct dictionary_entry{
	char* word;
	double score;
} *dictionary_entry;
list_of_coins * parse_coins_query(FILE*, int * );
users_and_tweets  parse_data(FILE*);
cluster * parse_data_2(tweet* , int , int);
unique get_unique(int *, int *, int );
void calculate_scores(tweet * , int);
double euclidean_distance(score_array * , score_array * , int );
double sample_gaussian();
double sample_uniform(int );
double * create_random_vector(int );
int binary_string_to_num(char* );
double cosine_similarity(score_array * , score_array * , int );
int** cosine_lsh_implementation(user* , int , int , int);
double** init(hashtable* , user * , int , int , int, int , int , int , int , double ** ,  int ** );
int cosine_hash(user , int , int , int , double** );
int** cosine_lsh_query(user* , hashtable , int ,  int , int , int , int , int , int , double** , int** , double** , int);
int cmp(const void *, const void *);

centroid * init_centroids(int , user* , int , int , int );
void kmeans(centroid* , user* , int , int , int , int );
int check_convergence(centroid * , centroid* , int , int);
void copy_user(user* , user , int );
void copy_centroid(centroid* , centroid , int );
void compute_cluster(centroid * , user * , int , int , int );
void lloyds_assignment(centroid * , user* , int , int , int );
void basic_update(centroid *, user* , int , int );
void shilouette_evaluation(centroid* , user* , int , int, int );
void k_means_recommend(centroid* , user* , int , int , int , int );

void clustering_lsh_implementation(user* , cluster* , int , int , int);
void clustering_cluster_implementation(user* , cluster* , int , int , int);
