#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h> 
#include <time.h>

typedef struct nearest_neighbor{
	long int id;
	double distance;
	int table;
	struct nearest_neighbor *next;
} *nearest_neighbor;

typedef struct nn_list{
	nearest_neighbor first, last;
} *nn_list;

typedef struct point{
	long int id;
	double * coordinates;
	long long int ** g_functions;
	int centroid_id;
	int centroid2_id;
	double dist; //distance from centroids
	double dist_as_centroid;
	struct point * next;
	double silhouette;
} *point;

typedef struct centroid{
	long int id;
	int count;
	int prev_count;
	double dist;
	long long int ** g_functions;
	double * coordinates;
	double silhouette_of_cluster;
	point* assigned_points;
} *centroid;


typedef struct cluster{
	point* data;
	centroid* centroids;
} *cluster;


typedef struct bucket{
	point first, last;
} *bucket;

typedef struct hashtable{
	bucket buckets;
	int cnt;
	int size;
} * hashtable;

void hashtable_init(hashtable *, int );
int hashtable_size(hashtable);
void hashtable_insert(hashtable*, point, long long int);
void hashtable_print(hashtable);
void hashtable_free(hashtable *);
void bucket_print(hashtable , int);
void nn_list_insert(nn_list * , long int , double, int);
