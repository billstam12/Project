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
centroid * init_centroids(int , int );
cluster compute_cluster(centroid * , point * , int, int, int );