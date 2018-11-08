#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
	struct point * next;
} *point;

typedef struct centroid{
	long int id;
	int count;
	double * coordinates;
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
