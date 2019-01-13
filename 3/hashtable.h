#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h> 
#include <time.h>



/* USEFUL STRUCTS FOR THE PROBLEM */
typedef struct node {
    int val;
    struct node * next;
} node_t;

typedef struct tweet2{
	double * coordinates;
	int id;
} *tweet2;

typedef struct tweet{
	int id;
	int user_id;
	int no_of_words;
	node_t * coin_list;
	double * coordinates;
	char ** words;  
	int cluster_id;
	double score;
} *tweet;

typedef struct cluster{
	tweet* tweets;
	int size;
} * cluster;

typedef struct score_array{
	double  value;
	double old_value; //used for when need to re-initialize
	int  id;
} * score_array;

typedef struct user{
	int id;
	int no_of_tweets;
	struct user * next;
	tweet * tweets;
	score_array * score_vector;

	//centroid stuff
	int centroid_id;
	int centroid2_id;
	double dist; //distance from centroids
	double dist_as_centroid;
	double silhouette;
} *user;

typedef struct centroid{
	long int id;
	int count;
	int prev_count;
	double dist;
	double silhouette_of_cluster;
	user* assigned_users;
	score_array * score_vector;
} *centroid;


typedef struct users_and_tweets{
	int no_of_users;
	user * users;
	int no_of_tweets;
	tweet * tweets;
	int no_of_clusters;
	cluster* clusters;
	int no_of_coins;
	char ** dict;
} *users_and_tweets;

/* HASHTABLE STUFF */
typedef struct nearest_neighbor{
	long int id;
	double distance;
	int table;
	struct nearest_neighbor *next;
} *nearest_neighbor;

typedef struct nn_list{
	nearest_neighbor first, last;
} *nn_list;


typedef struct bucket{
	user first, last;
} *bucket;

typedef struct hashtable{
	bucket *buckets;
	int cnt;
	int size;
} * hashtable;

/* LINKED LIST */
void print_list(node_t * );
void push(node_t * , int );
int not_in_list(node_t * , int );
void free_list(node_t *);

/* HASHTABLE*/

void hashtable_init(hashtable *, int );
int hashtable_size(hashtable);
void hashtable_insert(hashtable*, user, long long int);
void hashtable_print(hashtable);
void hashtable_free(hashtable *, int);
void bucket_print(hashtable , int);
