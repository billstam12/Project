#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include <time.h>
#include "hashtable.h"


double sample_gaussian();
double sample_uniform(int );
double * create_random_vector(int );
int binary_string_to_num(char* );
void print_coordinates(point , int );
long long int euclidean_hash(point *,  int , int ,  int ,int, int,   double* , int* , double** );
int cosine_hash(point* , int , int , int , double** );
void hypercube_query(point* , point* , hashtable , int , int , int , double , int , int , int , int , double** );
double euclidean_distance(double * , double * , int );
int compare_gfuncs(long long int *, long long int*, int);
double** init(hashtable*, point * , int , int , int  , int , int , int, int, int, double ** ,  int ** );
double print_euclidean_results(FILE**, nn_list, double);
double print_cosine_results(FILE**, nn_list, double);
double print_hypercube_results(FILE **, nn_list , double );
void euclidean_lsh_query(point* , point* , hashtable* , int, int , int , int , double , int, int, int, double** , int** , double***);
void cosine_query(point* , point* , hashtable* , int , int , int , int , double , int , int , int , double** , int** , double*** );
point* parse_data(FILE* , int* , int* ,int*);
point* parse_query(FILE* ,int*, int, double*);