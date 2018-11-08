#include "functions.h"







void parse_conf(FILE* f, int* no_clusters, int* no_functions, int* no_hashtables){
	rewind(f);
	char * line = NULL;
	size_t len = 0;
	int lines = 0;
	ssize_t read;
	while((read = getline(&line, &len, f)) !=  -1){
		lines ++;
	}
	rewind(f);
	int i= 0;
	while(i<lines && ((read = getline(&line, &len, f))!=0)) {
		line = strtok(line," ");
		if(!strcmp(line,"number_of_clusters:")){
			*no_clusters =  atoi(strtok(NULL,"\n"));
		}
		else if (!strcmp(line,"number_of_hash_functions:")){
			*no_functions = atoi(strtok(NULL,"\n"));
		}
		else if (!strcmp(line,"number_of_hash_tables:")){
			*no_hashtables = atoi(strtok(NULL,"\n"));
		}
		else{
			fprintf(stderr, "Error! Wrong configuration Format");
			exit(-1);
		}
		i++;
	}
	if((*no_clusters)==0){
		fprintf(stderr, "Error! Clusters cannot be 0\n");
		exit(-1);

	}

}

double euclidean_distance(double * c1, double * c2, int no_dimensions){
	int i;
	double dist = 0.0;
	for(i = 0; i < no_dimensions; i++){
		dist += pow((c1[i] - c2[i]), 2);
	}	
	return sqrt(dist);
}


double cosine_similarity(double * c1, double * c2, int no_dimensions){
	int i;
	double similarity = 0.0;
	double c1_dist = 0.0;
	double c2_dist = 0.0;
	for(i = 0; i < no_dimensions; i++){
		c1_dist += pow(c1[i], 2);
		c2_dist += pow(c2[i], 2);
	}
	c1_dist = sqrt(c1_dist); c2_dist = sqrt(c2_dist);

	for(i = 0; i < no_dimensions; i++){
		similarity += (c1[i] * c2[i]);
	}
	return similarity/(c1_dist*c2_dist);

}


centroid * init_centroids(int k, int no_dimensions){
	centroid *centroids = (centroid*)malloc(k*sizeof(centroid));
	int i;
	for(i = 0; i< k ; i++){
		centroid new_centroid;
		new_centroid = (centroid)malloc(sizeof(struct centroid));
		centroids[i] = new_centroid;
		centroids[i]->id = i;
		centroids[i]->count = 0;
		centroids[i]->assigned_points = (point*)malloc(sizeof(point));
		centroids[i]->coordinates = create_random_vector(no_dimensions);
	}
	return centroids;
}

cluster compute_cluster(centroid * centroids, point * data, int no_of_samples, int no_of_dimensions, int k){
	int i,j;
	double min;
	long int id;
	for(i = 0; i < no_of_samples; i++){
		min = euclidean_distance((centroids[0])->coordinates, (data[i])->coordinates, no_of_dimensions);
		id = (centroids[0])->id;
		double dist;

		for(j = 1; j < k; j++){
			dist = euclidean_distance((centroids[j])->coordinates, (data[i])->coordinates, no_of_dimensions);
			if(dist < min){
				min = dist;
				id = (centroids[j])-> id;
			}
		}
		(data[i])->centroid_id = id;

		(centroids[id])->count +=1;
		if((centroids[id])->count == 1){
			(centroids[id])->assigned_points[0] = data[i];
		}
		else{
			centroids[id]->assigned_points = (point*)realloc(centroids[id]->assigned_points,(centroids[id]->count)*sizeof(point));
			centroids[id]->assigned_points[(centroids[id]->count)-1] = data[i];
		}
		//printf("Point with ID:%ld was classified in %ld with distance %f\n",data[i]->id, id, dist);
	}
	cluster cl;
	cl = (cluster)malloc(sizeof(struct cluster));
	cl->data = data;
	cl->centroids = centroids;

	return cl;	
}


point * parse_data(FILE* f, int* no_samples, int* no_dimensions){
	rewind(f);
	char delimit[]="\t,";

	char * line = NULL;
	size_t len = 0;
	ssize_t read = getline(&line, &len, f);
	
	//Get no of samples
	int samps = 0;
	while((read = getline(&line, &len, f)) !=  -1){
		samps ++;
	}

	*no_samples = samps;
	rewind(f);
	read = getline(&line, &len, f); //Rewind and read first line

	//Get no of dimensions
	int flag = 1;
	int dims = 0;
	int i = 0;
	
	read = getline(&line, &len, f);
	if(flag == 1){
		flag = 0;
		while(i < strlen(line)){
			if(line[i] == '\t' || line[i] == ','){
				dims ++;
			}
			i++;
		}
	}
	*no_dimensions = dims;
	rewind(f);
	read = getline(&line, &len, f); //Rewind and read first line

	/*  First create an array of points so we can store 
		each of the points in the memory. Then, we will begin
		hashing each point from our table of point-pointers 
		to the hashtable.
		So initialize the array of point-pointers.
	*/
	point * data;
	data = (point *)malloc(samps * sizeof(point));

	int j;
	i = 0;
	while((read = getline(&line, &len, f)) !=  -1 ){
		point new_point;
		new_point = malloc(sizeof(struct point));
		new_point->coordinates = (double*) malloc(dims*sizeof(double));

		j = 0;		
		line = strtok(line,delimit);
		new_point-> id = atof(line);

		while((line = strtok(NULL,delimit)) != NULL){		
			new_point -> coordinates[j] = atof(line);
			data[i] = new_point;
			j++;
		}
		i++;
	}

	return data;
}

int binary_string_to_num(char* key){
	char* start = &key[0];
	int total = 0;
	while (*start){
		total *= 2;
		if (*start++ == '1') total += 1;
	}
	return total;
}

long long int real_modulo(long long int x, long long int n){
	return ((x % n + n ) %n );
}

double sample_gaussian() {
    double u = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double v = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double r = u * u + v * v;
    if (r == 0 || r > 1) return sample_gaussian();
    double c = sqrt(-2 * log(r) / r);
    return u * c;
}

double sample_uniform(int w) {
	return w * (double)rand() / ((double)RAND_MAX + 1); 
}

double * create_random_vector(int no_dimensions){
	float M = 50;
	float N = 150;
	int i;
	int *distribs = malloc(no_dimensions*sizeof(int));
	for(i=0; i < no_dimensions; i++){
		distribs[i] = M + (rand() / (RAND_MAX + 1.0)) * ( N - M + 1.0); // random int between M and N
		//printf("%d\n", distribs[i]);
	}
	double *r = malloc(no_dimensions*sizeof(double));
	for(i=0; i < no_dimensions; i++){
		r[i] = (distribs[i] - 1) + (rand()/ (RAND_MAX  + 1.0f))*( - 50.0f);
	}
	
	return r;
}
void print_coordinates(point p, int no_dimensions){
	int i;
	for( i = 0; i < no_dimensions; i++)
		printf("%f ", p->coordinates[i]);	
	printf("\n");
}

void print_coordinates_cent(centroid p, int no_dimensions){
	int i;
	for( i = 0; i < no_dimensions; i++)
		printf("%f ", p->coordinates[i]);	
	printf("\n");
}
