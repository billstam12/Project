#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h> 
#include <time.h>

typedef struct tweet{
	int id;
	int user_id;
	int no_of_words;
	char ** words;  
	int score;
} *tweet;



typedef struct user{
	long int id;
	tweet * tweets;
	
} *user;