#include "functions.h"

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
	double * r;
	r = (double *)malloc(no_dimensions*sizeof(double));
	int i;
	for(i=0;i<no_dimensions;i++){
		r[i] = sample_gaussian();
	}

	return r;
}
void print_coordinates(point p, int no_dimensions){
	int i;
	for( i = 0; i < no_dimensions; i++)
		printf("%f ", p->coordinates[i]);	
	printf("\n");
}

int cosine_hash(point* p, int no_dimensions, int no_of_functions, int table_size, double** random_vectors_temp){
	int index = 0;
	int i, j;
	char * hashes;
	hashes = (char *)malloc(no_of_functions * (sizeof( char )+1));
	long double first_hash = 0.0;
	for(i = 0; i < no_of_functions; i++){
		for(j = 0; j < no_dimensions; j++){
			first_hash += ((*p) -> coordinates[j]) * random_vectors_temp[i][j];
		}
		if(first_hash >= 0){
			hashes[i] = '1';
		}
		else{
			hashes[i] = '0';
		}
	}
	hashes[no_of_functions] = '\0';
	index = binary_string_to_num(hashes);

	return index;
}

long long int euclidean_hash(point* p, int no_dimensions, int window,  int no_of_functions, int table_size, int table_id, double* t, int *r, double** random_vectors_temp){
	long long int index = 0;
	int i, j;
	long long int * first_hashes;
	first_hashes = (long long int *)malloc(no_of_functions * sizeof(long long int ));
	long double first_hash = 0.0;
	for(i = 0; i < no_of_functions; i++){
		for(j = 0; j < no_dimensions; j++){
			first_hash += ((((*p) -> coordinates[j] * random_vectors_temp[i][j]) + t[i])/window);
		}
		
		first_hashes[i] = (long long int)floor(first_hash);
	
	}
	(*p) -> g_functions[table_id] = first_hashes;
	// Now that we have the g1(x), g2(x)... gk(x), we will create the φ function = index
	for(i = 0; i < no_of_functions; i++){
		index += (r[i] * (first_hashes[i]));
	}

	if(index > 0)
		index = (index%table_size);
	else
		index = real_modulo(index,table_size);
	
	//printf("%lld\n",index);
	return index;
}

int compare_gfuncs(long long int * g1, long long int * g2, int no_of_functions){
	int i = 0;
	int b = 1;
	while( (b == 1) && (i < no_of_functions) ){
		b = (g1[i] == g2[i]);
		i++;
	}
	return b;
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

double** init(hashtable* ht, point *data , int no_samples, int no_dimensions, int window, int no_of_functions, int table_size,int L, int table_id, int metric, double ** t,  int ** r){
	// Initialize window t and number of functions
	int i,j;
	double *t_temp = (double*)malloc(no_of_functions*sizeof(double));
	for(i = 0; i < no_of_functions; i++){
		t_temp[i] = sample_uniform(window);
	}	
	
	(*t) = t_temp;
	int *r_temp = (int*)malloc(no_of_functions*sizeof(int));
	for(i = 0; i < no_of_functions; i++){
		r_temp[i] = rand() % 100;
	}

	//Create r-scalars array
	(*r) = r_temp;
	// Create matrix of k random vectors and transpose it so we can multiply with data
	double **random_vectors_temp = (double**)malloc(no_of_functions* sizeof(double*));
	for(i = 0; i < no_of_functions; i++){
		random_vectors_temp[i] = create_random_vector(no_dimensions);
	}
	
	long long int index;
	for(i = 0; i < no_samples; i++){
		if(metric == 0){
			index = euclidean_hash(&data[i], no_dimensions, window, no_of_functions, table_size, table_id, t_temp, r_temp, random_vectors_temp);
			/* initialize g functions here not in the hash */
			hashtable_insert(ht, data[i], index);
		}
		else{
			index = cosine_hash(&data[i], no_dimensions, no_of_functions, table_size,  random_vectors_temp);
			/* initialize g functions here not in the hash */
			hashtable_insert(ht, data[i], index);
		}
	}
	return random_vectors_temp;

}

void euclidean_lsh_query(point* data, point* query_data, hashtable* hts, int no_queries, int no_samples, int no_dimensions, int L, double radius, int window, int no_of_functions, int table_size, double** t, int** r, double*** random_vectors){
	int i, j;
	FILE * o;
	o = fopen("output/euclidean.txt", "w");
	for(j = 0; j < no_queries; j++){
		/* L S H */
		clock_t start, end;
	    double cpu_time_used;
	     
	    start = clock();
	    
		fprintf(o,"-------------------------------------------------------\n");
		fprintf(o,"Query: %ld\n", query_data[j]->id);
		fprintf(o,"R-near neighbors:\n");
		query_data[j]->g_functions = (long long int**)malloc(sizeof(long long int)* L);
		nn_list nearest_neighbors = (struct nn_list *)malloc(sizeof(struct nn_list));
		int count = 0;
		for(i = 0; i < L; i++){
			
			//printf("TABLE NO: %d\n", i);
			long long int index = euclidean_hash(&query_data[j], no_dimensions,  window,   no_of_functions,  table_size, i, t[i], r[i],  random_vectors[i]);			
			point tmp;
			tmp = hts[i]->buckets[index].first;
 			
			while(tmp != NULL){
				int g_func;
				g_func = compare_gfuncs(tmp->g_functions[i], query_data[j]->g_functions[i], no_of_functions);
				if(g_func == 1){
					count ++;
					double distance = euclidean_distance(tmp->coordinates, query_data[j]->coordinates, no_dimensions);
					if(distance < radius && distance!=0){
						nn_list_insert(&nearest_neighbors, tmp->id, distance, i);
					}
				}
				tmp = tmp -> next;
			}
			free(tmp);
		}
		end = clock();
	    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		double min_pr = print_euclidean_results(&o, nearest_neighbors, cpu_time_used);			
		//printf("COUNT=%d\n",count); //ABOUT 1/4th of the total

	    /* Deterministic */
	    double min = euclidean_distance(data[0]->coordinates, query_data[j]->coordinates, no_dimensions);
	    long int found_id = data[0]->id;
	   	start = clock();
	    for (i = 1; i < no_samples; i++){
	    	double dist;
	    	dist = euclidean_distance(data[i]->coordinates, query_data[j]->coordinates, no_dimensions);
	    	if(dist < min && dist != 0){
	    		min = dist;
	    		found_id = data[i]->id;
	    	}
	    }
	    end = clock();
	    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	    fprintf(o,"Nearest neighbor Deterministic: %ld\n", found_id);
		fprintf(o,"Nearest neighbor Deterministic distance: %f\n", min);
		fprintf(o,"Nearest neighbor Deterministic time: %f\n", cpu_time_used);
		fprintf(o,"Probabilistic Distance by Deterministic Distance: %f\n", min_pr/min);
	}
}

void cosine_query(point* data, point* query_data, hashtable* hts, int no_queries, int no_samples, int no_dimensions, int L, double radius, int window, int no_of_functions, int table_size, double** t, int** r, double*** random_vectors){
	int i, j;
	FILE * o;
	o = fopen("output/cosine.txt", "w");
	for(j = 0; j < no_queries; j++){
		/*C o s i n e */
		clock_t start, end;
	    double cpu_time_used;
	     
	    start = clock();
	    
		fprintf(o,"-------------------------------------------------------\n");
		fprintf(o,"Query: %ld\n", query_data[j]->id);
		fprintf(o,"R-near neighbors:\n");
		query_data[j]->g_functions = (long long int**)malloc(sizeof(long long int)* L);
		nn_list nearest_neighbors = (struct nn_list *)malloc(sizeof(struct nn_list));
		int count = 0;
		for(i = 0; i < L; i++){
			
			//printf("TABLE NO: %d\n", i);
			long long int index = cosine_hash(&query_data[j], no_dimensions, no_of_functions,  table_size, random_vectors[i]);			
			point tmp;
			tmp = hts[i]->buckets[index].first;
 			
			double distance, similarity;
			while(tmp != NULL){
				count++;	
	    		similarity = cosine_similarity(data[i]->coordinates, query_data[j]->coordinates, no_dimensions);
				distance = 1 - similarity;
				if(distance < radius && distance!=0){
					nn_list_insert(&nearest_neighbors, tmp->id, distance, i);
				}
				tmp = tmp -> next;
			}
			free(tmp);
		}
		end = clock();
	    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		double min_pr = print_cosine_results(&o, nearest_neighbors, cpu_time_used);			
		//printf("COUNT=%d\n",count); //ABOUT 1/2th of the total

	    /* Deterministic */
	    double min = euclidean_distance(data[0]->coordinates, query_data[j]->coordinates, no_dimensions);
	    long int found_id = data[0]->id;
	   	start = clock();
	    for (i = 1; i < no_samples; i++){
	    	double distance, similarity;
	    	similarity = cosine_similarity(data[i]->coordinates, query_data[j]->coordinates, no_dimensions);
	    	distance = 1- similarity;
	    	if(distance < min && distance != 0){
	    		min = distance;
	    		found_id = data[i]->id;
	    	}
	    }
	    end = clock();
	    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	    fprintf(o,"Nearest neighbor Deterministic: %ld\n", found_id);
		fprintf(o,"Nearest neighbor Deterministic distance: %f\n", min);
		fprintf(o,"Nearest neighbor Deterministic time: %f\n", cpu_time_used);
		fprintf(o,"Probabilistic Distance by Deterministic Distance: %f\n", min_pr/min);

	}
}

void hypercube_query(point* data, point* query_data, hashtable hypercube , int no_queries, int no_samples, int no_dimensions, double radius, int no_of_functions, int table_size, int probes, int max, double** random_vectors){
	int i, j;
	FILE * o;
	o = fopen("output/hypercube.txt", "w");
	for(j = 0; j < no_queries; j++){
		/* L S H */
		clock_t start, end;
	    double cpu_time_used;
	     
	    start = clock();
	    
		fprintf(o,"-------------------------------------------------------\n");
		fprintf(o,"Query: %ld\n", query_data[j]->id);
		fprintf(o,"R-near neighbors:\n");
		nn_list nearest_neighbors = (struct nn_list *)malloc(sizeof(struct nn_list));
			
		//printf("TABLE NO: %d\n", i);
		long long int index = cosine_hash(&query_data[j], no_dimensions, no_of_functions,  table_size, random_vectors);			
		point tmp;
		tmp = hypercube->buckets[index].first;
 		
 		int count = 0;
		while(tmp != NULL ){
			count ++;
			if(count <= max){
				double distance = euclidean_distance(tmp->coordinates, query_data[j]->coordinates, no_dimensions);
				if(distance < radius && distance!=0){
					nn_list_insert(&nearest_neighbors, tmp->id, distance, i);
				}
				
			}
			else break;
			tmp = tmp -> next;
		}
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		double min_pr = print_hypercube_results((&o), nearest_neighbors, cpu_time_used);			

		/* Deterministic */
		double min = euclidean_distance(data[0]->coordinates, query_data[j]->coordinates, no_dimensions);
		long int found_id = data[0]->id;
		start = clock();
		for (i = 1; i < no_samples; i++){
		  	double dist;
		   	dist = euclidean_distance(data[i]->coordinates, query_data[j]->coordinates, no_dimensions);
		   	if(dist < min && dist != 0){
		   		min = dist;
		   		found_id = data[i]->id;
		   	}
		}
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		fprintf(o,"Nearest neighbor Deterministic: %ld\n", found_id);
		fprintf(o,"Nearest neighbor Deterministic distance: %f\n", min);
		fprintf(o,"Nearest neighbor Deterministic time: %f\n", cpu_time_used);
		fprintf(o,"Probabilistic Distance by Deterministic Distance: %f\n", min_pr/min);
	}
}


double print_cosine_results(FILE **o, nn_list nn, double lsh_tm){
	fprintf((*o), "Cosine\n");
	nearest_neighbor tmp;
	tmp = nn->first;
	double min = 0;
	if(tmp == NULL){
		fprintf((*o), "Cosine LSH found no nearest neighbors for point!\n");
	}
	else{
		while(tmp != NULL){
			fprintf((*o), "Item: %ld ", tmp->id);
			fprintf((*o), "Distance: %f ", tmp->distance);
			fprintf((*o), "Table: %d\n",tmp->table);
			tmp = tmp -> next;
		}
		tmp = nn->first;
		// Find nearest neighbor
		min = tmp->distance;
		long int nn_id = tmp->id;
		while(tmp != NULL){
			if(tmp->distance < min){
				min = tmp->distance;
				nn_id = tmp->id;
			}
			tmp = tmp -> next;
		}
		fprintf((*o), "Nearest neighbor LSH(Cosine): %ld\n", nn_id);
		fprintf((*o), "Nearest neighbor LSH(Cosine) distance: %f\n", min);
		fprintf((*o), "Nearest neighbor LSH(Cosine) time: %f\n", lsh_tm);
	}
	fprintf((*o),"\n");
	return min;
}

double print_hypercube_results(FILE** o, nn_list nn, double lsh_tm){
	fprintf((*o), "Hypercube\n");
	nearest_neighbor tmp;
	tmp = nn->first;
	double min = 0;
	if(tmp == NULL){
		fprintf((*o), "Hypercube search found no nearest neighbors for point!\n");
	}
	else{
		while(tmp != NULL){
			fprintf((*o),"Item: %ld ", tmp->id);
			fprintf((*o),"Distance: %f\n", tmp->distance);
			tmp = tmp -> next;
		}
		tmp = nn->first;
		// Find nearest neighbor
		min = tmp->distance;
		long int nn_id = tmp->id;
		while(tmp != NULL){
			if(tmp->distance < min){
				min = tmp->distance;
				nn_id = tmp->id;
			}
			tmp = tmp -> next;
		}
		fprintf((*o),"Nearest neighbor Hypercube: %ld\n", nn_id);
		fprintf((*o),"Nearest neighbor Hypercube: %f\n", min);
		fprintf((*o),"Nearest neighbor Hypercube time: %f\n", lsh_tm);
	}
	fprintf((*o),"\n");
	return min;
}

double print_euclidean_results(FILE **o, nn_list nn, double lsh_tm){
	fprintf((*o), "Euclidean\n");
	nearest_neighbor tmp;
	tmp = nn->first;
	double min = 0;
	if(tmp == NULL){
		fprintf((*o), "Euclidean LSH found no nearest neighbors for point!\n");
	}
	else{
		while(tmp != NULL){
			fprintf((*o), "Item: %ld ", tmp->id);
			fprintf((*o), "Distance: %f ", tmp->distance);
			fprintf((*o), "Table: %d\n",tmp->table);
			tmp = tmp -> next;
		}
		tmp = nn->first;
		// Find nearest neighbor
		min = tmp->distance;
		long int nn_id = tmp->id;
		while(tmp != NULL){
			if(tmp->distance < min){
				min = tmp->distance;
				nn_id = tmp->id;
			}
			tmp = tmp -> next;
		}
		fprintf((*o), "Nearest neighbor LSH(Euclidean): %ld\n", nn_id);
		fprintf((*o), "Nearest neighbor LSH(Euclidean) distance: %f\n", min);
		fprintf((*o), "Nearest neighbor LSH(Euclidean) time: %f\n", lsh_tm);
	}
	fprintf((*o),"\n");
	return min;
}

point * parse_data(FILE* f, int* no_samples, int* no_dimensions, int* m){
	rewind(f);
	char * line = NULL;
	size_t len = 0;
	ssize_t read = getline(&line, &len, f);
	
	//Get metric
	char* metric = strtok(line, "\t");
	metric = strtok(NULL, "\n");
	if(!strcmp(metric,"cosine")){
		*m = 1;
	}

	else{ *m = 0; }
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
			if(line[i] == '\t'){
				dims ++;
			}
			i++;
		}
	}
	*no_dimensions = dims;
	rewind(f);
	read = getline(&line, &len, f); //Rewind and read first line

	
	/* Begin the hashing */

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
		line = strtok(line,"\t");
		new_point-> id = atof(line);

		while((line = strtok(NULL,"\t")) != NULL){		
			new_point -> coordinates[j] = atof(line);
			data[i] = new_point;
			j++;
		}
		i++;
	}

	// Now data contains the addresses of all of our points
	/*
	for(i = 0; i < samps; i++){
		for (int j = 0; j < dims;j++){
			printf("%f ",data[i] -> coordinates[j]);
		}
		printf("\n");
	}
	*/
	

	
	return data;
}

point * parse_query(FILE* f, int* no_queries, int no_dimensions, double* radius){
	rewind(f);
	char * line = NULL;
	size_t len = 0;
	ssize_t read = getline(&line, &len, f);
	
	//Get metric
	*radius = atoi(strtok(line, "\n"));

	//Get no of samples
	int samps = 0;
	while((read = getline(&line, &len, f)) !=  -1){
		samps ++;
	}

	*no_queries = samps;

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
			if(line[i] == '\t'){
				dims ++;
			}
			i++;
		}
	}
	if(dims != no_dimensions){
		fprintf(stderr, "Error! Query file components have different dimensions");
		exit(-1);
	}

	rewind(f);
	read = getline(&line, &len, f); //Rewind and read first line

	
	/* Begin the hashing */

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
		line = strtok(line,"\t");
		new_point-> id = atof(line);

		while((line = strtok(NULL,"\t")) != NULL){		
			new_point -> coordinates[j] = atof(line);
			data[i] = new_point;
			j++;
		}
		i++;
	}

	// Now data contains the addresses of all of our points
	/*
	for(i = 0; i < samps; i++){
		for (int j = 0; j < dims;j++){
			printf("%f ",data[i] -> coordinates[j]);
		}
		printf("\n");
	}
	*/

	return data;
}