#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include <assert.h>
#include <time.h>
#include "hashtable.h"


point * parse_data(FILE* , int* , int* );
void parse_conf(FILE* , int* , int* , int* );
double sample_gaussian();
double sample_uniform(int );
double * create_random_vector(int );
int binary_string_to_num(char* );
void print_coordinates(point , int );
void print_coordinates_cent(centroid , int );
double euclidean_distance(double * , double * , int );
double cosine_similarity(double * , double * , int );
centroid * init_centroids(int , point* , int , int, int );
void compute_cluster(centroid * , point * , int, int, int*, int,  int,int, int, int );
int check_convergence(centroid * , centroid* , int , int , int);
void copy_centroids(centroid** , centroid *, int , int );
void copy_centroid(centroid* , centroid, int );
void copy_point(point* , point , int );
void kmeans(centroid* ,point* , int , int ,int , int, int, int, int, int );
int my_rand(double* , int , int);
void basic_update(centroid *, point* , int , int * );
void pam(centroid *, point* , int , int *);
void lloyds_assignment(centroid * , point* , int , int , int *);
void lsh_assignment(centroid* , point* , int , int , int, int, int, int *);
double** init(hashtable*, point * , int , int , int  , int , int , int, int, int, double ** ,  int ** );
void euclidean_lsh_query(point* , centroid* ,  hashtable* , int , int , int , int , double , int , int , int , double** , int** , double*** );
void cosine_query(point* , point* , char*, hashtable* , int , int , int , int , double , int , int , int , double** , int** , double*** );
int cosine_hash(point* , int , int , int , double** );
long long int euclidean_hash(point*, int , int ,  int , int , int , double* , int *, double** );
long long int euclidean_hash_centroid(centroid*, int , int ,  int , int , int , double* , int *, double** );
int compare_gfuncs(long long int * , long long int * , int );
void shilouette_evaluation(centroid* ,point* ,int ,int );
void *change_mem(void *, size_t , size_t );
int check_lsh_convergence(int *,int *, int);