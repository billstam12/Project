#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 
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
centroid * init_centroids(int , point* , int , int );
void compute_cluster(centroid * , point * , int, int, int* );
int check_convergence(centroid * , centroid* , int , int , int);
void copy_centroids(centroid** , centroid *, int , int );
void copy_centroid(centroid* , centroid, int );
void copy_point(point* , point , int );
void kmeans(centroid* ,point* , int , int ,int );