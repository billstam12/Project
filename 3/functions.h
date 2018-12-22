#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include <assert.h>
#include <time.h>
#include "hashtable.h"

typedef struct unique{
	int count;
	int* array;
	int* tweets_per_user;
}  unique;

tweet * parse_data(FILE*, int * );
unique get_unique(int *, int *, int );


