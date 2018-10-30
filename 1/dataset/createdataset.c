#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// INPUT IS NUMBER OF DATA AND NUMBER OF DIMENSIONS

int main(int argc, char** argv){
	//Create new random dataset file
	// Default values
	int no_of_samples = 100;
	int no_of_dim = 10;
	int type = 0; //0 for dataset 1 for query
	int i, j;
	for(i = 0; i < argc;i++){
		if(!strcmp(argv[i],"-n")){
			no_of_samples = atoi(argv[++i]);	
		}
		else if(!strcmp(argv[i], "-d")){
			no_of_dim = atoi(argv[++i]);
		}
		else if(!strcmp(argv[i], "-t")){
			type = atoi(argv[++i]);
		}
	}

	float M = 50;
	float N = 150;
	srand(3);
	int *distribs = malloc(no_of_dim*sizeof(int));
	for(i=0; i < no_of_dim; i++){
		distribs[i] = M + (rand() / (RAND_MAX + 1.0)) * ( N - M + 1.0); // random int between M and N
		//printf("%d\n", distribs[i]);
	}
	long id;
	if(type == 0){
		FILE *f = fopen("dataset.csv", "w");
		
		if (f == NULL)
		{
			printf("Error opening file!\n");
			exit(1);
		}
		double r;
		
		
		
		// Metric to use 33% cosine, 33% euclid, 33% nothing

		int rn = rand() % 100;
		printf("%d", rn);
		if(rn < 33){
			fprintf(f, "@metric cosine\n");
		}
		else if(rn < 66){
			fprintf(f, "@metric euclidean\n");
		}
		else{
			fprintf(f, "\n");
		}
		id = 0;

		// Create dimension coordinates that are close to each other numerically
		for(i=0;i < no_of_samples; i++){
			for(j=0; j < no_of_dim-1; j++){
				if(j == 0){
					fprintf(f,"%ld\t", id);
					id ++;
				}
				r = (distribs[j] - 1) + (rand()/ (RAND_MAX  + 1.0f))*( - 50.0f);
				fprintf(f,"%f\t", r);
			}
			r = (distribs[no_of_dim-1] - 1) + (rand()/ (RAND_MAX  + 1.0f))*( - 5.0f);
			fprintf(f,"%f\n", r);		
		}
	}
	else if(type == 1){

		id = 10000;

		FILE *f = fopen("query.csv", "w");
		
		if (f == NULL)
		{
			printf("Error opening file!\n");
			exit(1);
		}
		double r;
		
		// Metric to use 33% cosine, 33% euclid, 33% nothing

		double radius = rand() % 100;
		
		fprintf(f, "%f\n",radius);
		long id = 0;
		// Create a no_of_dim dimension array, which every element x will be the [x-2, x]
		// distribution of each coordinate.
		

		// Create dimension coordinates that are close to each other numerically
		for(i=0;i < no_of_samples; i++){
			for(j=0; j < no_of_dim-1; j++){
				if(j == 0){
					fprintf(f,"%ld\t", id);
					id ++;
				}
				r = (distribs[j] - 1) + (rand()/ (RAND_MAX  + 1.0f))*( - 50.0f);
				fprintf(f,"%f\t", r);
			}
			r = (distribs[no_of_dim-1] - 1) + (rand()/ (RAND_MAX  + 1.0f))*( - 5.0f);
			fprintf(f,"%f\n", r);		
		}
	}
	else{
		fprintf(stderr, "Type must be 0 or 1! Exiting...\n");
    	exit(-1);
 
	}
	
	
	
	return 0;
}