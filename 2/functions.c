#include "functions.h"

/*DATA FUNCTIONS */
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
		new_point-> dist = 0.0;
		while((line = strtok(NULL,delimit)) != NULL){		
			new_point -> coordinates[j] = atof(line);
			data[i] = new_point;
			j++;
		}
		i++;
	}

	return data;
}

/* DISTANCE FUNCTIONS */
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

/*GENERAL FUNCTIONS */
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

void print_coordinates_cent(centroid p, int no_dimensions){
	int i;
	for( i = 0; i < no_dimensions; i++)
		printf("%f ", p->coordinates[i]);	
	printf("\n");
}

int my_rand(double* freq, int n, int id){
    // Create and fill prefix array
    double *prefix;
    prefix  = (double*)malloc(n*sizeof(double));
    int i;
    prefix[0] = freq[0];
    for (i = 1; i < n; ++i){
        prefix[i] = prefix[i - 1] + freq[i];
    }
 
    // prefix[n-1] is sum of all frequencies. Generate a random number
    // with value from 1 to this sum
    int r = (rand() % (int)prefix[n - 1]) + 1;
    // Find index of ceiling of r in prefix arrat
    return r;
}


void copy_centroid(centroid* c1, centroid c2, int no_of_dimensions){
	int i;
	(*c1)->id = c2->id;
	(*c1)->count = c2->count;
	(*c1)->coordinates = (double*)malloc(no_of_dimensions*sizeof(double));
	for(i = 0; i < no_of_dimensions; i++){
		(*c1)->coordinates[i] = c2->coordinates[i];
	}
	(*c1)->assigned_points = (point*)malloc((*c1)->count*sizeof(point));
	for(i = 0; i < (*c1)->count; i++){
		copy_point(&((*c1)->assigned_points[i]), c2->assigned_points[i], no_of_dimensions);
	}
}

void copy_point(point* p1, point p2, int no_of_dimensions){
	int i;
	(*p1)->id = p2->id;
	(*p1)->centroid_id = p2->centroid_id;
	printf("%ld\n", (*p1)->id);
	(*p1)->coordinates = (double*)malloc(no_of_dimensions*sizeof(double));
	for(i = 0; i < no_of_dimensions; i++){
		(*p1)->coordinates[i] = p2->coordinates[i];
	}
}


/* INITIALIZATION FUNCTIONS */


centroid * init_centroids(int k, point* data, int no_of_samples, int no_of_dimensions, int type){
	centroid *centroids = (centroid*)malloc(k*sizeof(centroid));
	int i,j,z;
	if(type==0){
		for(i = 0; i< k ; i++){
			centroid new_centroid;
			new_centroid = (centroid)malloc(sizeof(struct centroid));
			centroids[i] = new_centroid;
			centroids[i]->id = i;
			centroids[i]->count = 0;
			centroids[i]->prev_count = 1;
			centroids[i]->assigned_points = (point*)malloc(sizeof(point));
			centroids[i]->coordinates = data[rand()%no_of_samples]->coordinates;
		}
	}
	else{
		int id = (int)sample_uniform(no_of_samples);
		int found_centers = 0; //calculate how many centers we have found
		// Get first random centroid
		centroid new_centroid;
		new_centroid = (centroid)malloc(sizeof(struct centroid));
		centroids[0] = new_centroid;
		centroids[0]->count = 0;
		centroids[0]->prev_count = 1;
		centroids[0]->assigned_points = (point*)malloc(sizeof(point));
		centroids[0]->coordinates = data[id]->coordinates;
		found_centers = 1;

		
		//Compute Distances from each centroid
		for(z = 0; z < k; z++){
			double* freqs = (double*)malloc(sizeof(double)*no_of_samples);
			double* probs = (double*)malloc(sizeof(double)*no_of_samples);
			double freqs_sum = 0;
			for(j = 0; j < no_of_samples; j++){
				if(j!=id){
					double dist;
					for(i = 0; i < found_centers; i++){
						data[j]-> dist += euclidean_distance(centroids[i]->coordinates, data[j]->coordinates, no_of_dimensions);
					}
					freqs[j] = pow(data[j]->dist,2);
					freqs_sum += freqs[j];
				}
			}
			for( j = 0; j < no_of_samples; j++){
				if(j!=id){
					probs[j] = freqs[j]/(freqs_sum-freqs[j]);
				}
			}

			//Get random id based on probability distribution
			double max_prob=0;
			for( j = 0; j < no_of_samples; j++){
				if(max_prob < probs[j]){
					max_prob  = probs[j];
					id = j;
				}
			}
						
			centroid new_centroid;
			new_centroid = (centroid)malloc(sizeof(struct centroid));
			centroids[z+1] = new_centroid;
			centroids[z+1]->id = z+1;
			centroids[z+1]->count = 0;
			centroids[z+1]->prev_count = 1;
			centroids[z+1]->assigned_points = (point*)malloc(sizeof(point));
			centroids[z+1]->coordinates = data[id]->coordinates;
			found_centers++;
		}
	}
	centroids[0]->id = 0;
	return centroids;
}


/* ASSIGNMENT FUNCTIONS */
void euclidean_lsh_query(point* data, centroid* query_data, hashtable* hts, int no_queries, int no_samples, int no_dimensions, int L, double radius, int window, int no_of_functions, int table_size, double** t, int** r, double*** random_vectors){
	int i, j;
	for(j = 0; j < no_queries; j++){
		/* L S H */
		query_data[j]->g_functions = (long long int**)malloc(sizeof(long long int)* L);
		nn_list nearest_neighbors = (struct nn_list *)malloc(sizeof(struct nn_list));
		int count = 0;
		for(i = 0; i < L; i++){
			
			//printf("TABLE NO: %d\n", i);
			long long int index = euclidean_hash_centroid(&query_data[j], no_dimensions,  window,   no_of_functions,  table_size, i, t[i], r[i],  random_vectors[i]);			
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
	}
}

void lsh_assignment(centroid* centroids, point* data, int no_of_samples, int no_of_dimensions, int metric, int *k){
	int table_size = no_of_samples/12;
	int L = 4;
	int window = 30;
	int no_of_functions = 40;
	int i, j,z;
	double radius;
	/*Calculate radius */
	double min_dist= euclidean_distance(centroids[0]->coordinates, centroids[1]->coordinates, no_of_dimensions);
	for(i = 0; i < *k; i++){
		for(j = i+1; j < *k; j++){
			double dist = euclidean_distance(centroids[i]->coordinates, centroids[j]->coordinates, no_of_dimensions);
			if( dist < min_dist){
				min_dist = dist;
			}
		}
	}
	radius = min_dist /2;
	//printf("%f\n",radius);
	double ** t; // t array
	t = (double**)malloc(L*sizeof(double*));
	int ** r; // r array
	r = (int**)malloc(L*sizeof(int*));
	double *** random_vectors; // random hyperplane vectors
	random_vectors = (double***)malloc(L*sizeof(double**));
	hashtable * hts;
	hts = (hashtable*)malloc(L * sizeof(hashtable));

	if(metric == 0){
		for(i = 0; i < L; i ++){
			hashtable_init(&hts[i], table_size);
		}
		for(i=0;i<no_of_samples;i++){
			data[i]->g_functions = (long long int**)malloc(sizeof(long long int)* L);
		}	
		for(i = 0; i < L; i++){
			random_vectors[i] = init(&hts[i], data ,  no_of_samples,  no_of_dimensions, window,  no_of_functions, table_size, L,  i, metric, &t[i], &r[i]);
			hashtable_print(hts[i]);
		}

		/* Now we will give the queries to our hashtables */	
		euclidean_lsh_query(data,  centroids, hts,  *k, no_of_samples, no_of_dimensions,  L,  radius,  window,  no_of_functions,  table_size,  t,  r,  random_vectors);
	}
}

void *my_realloc(void *ptr, size_t originalLength, size_t newLength)
{
   // Note that because we cannot rely on implementation details of the standard library,
   // we can never grow a block in place like realloc() can. However, it is semantically
   // equivalent to allocate a new block of the appropriate size, copy the original data
   // into it, and return a pointer to that new block. In fact, realloc() is allowed to
   // do this, as well. So we lose a possible performance optimization (that is rarely
   // possible in practice anyway), but correctness is still ensured, and the caller
   // never need be the wiser.
   // 
   // Likewise, we cannot actually shrink a block of memory in-place, so we either
   // have to return the block unchanged (which is legal, because a block of memory
   // is always allowed to be *larger* than necessary), or allocate a new smaller
   // block, copy the portion of the original data that will fit, and return a
   // pointer to this new shrunken block. The latter would actually be slower,
   // so we'll avoid doing this extra work when possible in the current implementation.
   // (You might need to change this if memory pressure gets out of control.)

   if (newLength == 0)
   {
      free(ptr);
      return NULL;
   }
   else if (!ptr)
   {
      return malloc(newLength);
   }
   else if (newLength <= originalLength)
   {
      return ptr;
   }
   else
   {
      assert((ptr) && (newLength > originalLength));
      void *ptrNew = malloc(newLength);
      if (ptrNew)
      {
          memcpy(ptrNew, ptr, originalLength);
          free(ptr);
      }
      return ptrNew;
    }
}

void lloyds_assignment(centroid * centroids, point* data, int no_of_samples, int no_of_dimensions, int *k){
	int i,j,z;
	double min;
	int id, min2_id;

	for(i = 0; i < no_of_samples; i++){
		min = euclidean_distance((centroids[0])->coordinates, (data[i])->coordinates, no_of_dimensions);
		id = (centroids[0])->id;
		double dist;

		for(j = 1; j < (*k); j++){
			dist = euclidean_distance((centroids[j])->coordinates, (data[i])->coordinates, no_of_dimensions);
			if(dist < min){
				min2_id = id;
				min = dist;
				id = (centroids[j])-> id;
			}
		}
		(data[i])->centroid_id = id;
		(data[i])->centroid2_id = min2_id;

		(centroids[id])->count +=1;
		if((centroids[id])->count == 1){
			(centroids[id])->assigned_points[0] = data[i];
		}
		else{
			centroids[id]->assigned_points = (point*)my_realloc(centroids[id]->assigned_points,(centroids[id]->prev_count*sizeof(point)),(centroids[id]->count)*sizeof(point));//Temporal
			centroids[id]->assigned_points[(centroids[id]->count)-1] = data[i];
			centroids[id]->prev_count = centroids[id]->count;
		}
		//printf("Point with ID:%ld was classified in %ld with distance %f\n",data[i]->id, id, dist);
	}
}

/* UPDATE FUNCTIONS */
void basic_update(centroid *centroids, point* data, int no_of_dimensions, int *k){
// Compute cluster mean and find new centroid
	int i,j,z ;
	for(i = 0; i < (*k); i++){
		// If for some reason one cluster is empty, then remove it 
		if(centroids[i]->count==0){
			for(j=i;j<(*k);j++){
				centroids[j] = centroids[j+1];
			}
			i -= 1;
			*k -= 1;
		}
		else{
			double * mean = (double*)malloc(sizeof(double)*no_of_dimensions);
			for(j = 0; j < no_of_dimensions; j++){
				mean[j] = 0;
				for(z = 0; z < centroids[i]->count; z++){
					mean[j] += centroids[i]->assigned_points[z]->coordinates[j];
				}
				mean[j] /= centroids[i]->count;
			}
			centroids[i]->coordinates = mean;
		}
	}
}

void pam(centroid *centroids, point* data, int no_of_dimensions, int *k){
// Compute cluster cost function and find new centroid
	int i,j,z ;
	for(i = 0; i < (*k); i++){
		if(centroids[i]->count==0){
			for(j=i;j<(*k);j++){
				centroids[j] = centroids[j+1];
			}
			i -= 1;
			*k -= 1;
		}
		double min_cost = centroids[i]->dist;
		long int min_id;
		for( j = 0; j < centroids[i]->count; j++){
			centroids[i]->assigned_points[j]->dist_as_centroid = 0.0;
			for( z = 0; z < centroids[i]->count; z++ ){
				//Compute new cost
				if(z!=j){
					centroids[i]->assigned_points[j]->dist_as_centroid += euclidean_distance(centroids[i]->assigned_points[j]->coordinates, centroids[i]->assigned_points[z]->coordinates,no_of_dimensions);
				}
			}
			if(centroids[i]->assigned_points[j]->dist_as_centroid < min_cost){
				min_id =  j;
				min_cost = centroids[i]->assigned_points[j]->dist_as_centroid;
			}

		}
		//Update centroid with the new one
		for( j = 0; j< centroids[i]->count; j++){
			centroids[i]->coordinates[j] = centroids[i]->assigned_points[min_id]->coordinates[j];
		}
	}
}

/* K-MEANS */
void compute_cluster(centroid * centroids, point * data, int no_of_samples, int no_of_dimensions, int* k, int assignment, int metric, int type){
	int i,j,z;
	// Check where each point falls
	if(assignment == 0){
		lloyds_assignment(centroids, data,  no_of_samples,  no_of_dimensions, k);
	}
	else if(assignment == 1){
		lsh_assignment(centroids, data, no_of_samples, no_of_dimensions, metric, k);
	}
	for(i = 0; i < (*k); i++){
		centroids[i]->dist = 0.0;
		for(j = 0; j < centroids[i]->count; j++){
			centroids[i]->dist += centroids[i]->assigned_points[j]->dist;
		}
		
	}
	if(type == 0){
		basic_update(centroids, data, no_of_dimensions, k);
	}
	else{
		pam(centroids, data, no_of_dimensions, k);
	}
}


void kmeans(centroid* centroids,point* data, int no_of_samples, int no_of_dimensions,int k, int assignment, int metric, int type){
	int converged = 0;
	int i;
	while(converged == 0){
		// Save the centroids for later
		for( i = 0; i< k;i++){
			centroids[i] -> count = 0;
		}
		int k_prev = k;
		centroid * old_cent;
		old_cent = (centroid*)malloc(sizeof(centroid)*k);
		for( i = 0; i< k; i++){
			old_cent[i] = (centroid)malloc(sizeof(struct centroid));
			copy_centroid(&(old_cent[i]), centroids[i], no_of_dimensions );
		}
		compute_cluster(centroids, data, no_of_samples, no_of_dimensions, &k, assignment, metric, type);
		
		
		for(i = 0; i < k; i++){
			print_coordinates_cent(old_cent[i], no_of_dimensions);
			print_coordinates_cent(centroids[i], no_of_dimensions);		
			printf("\n");
		}
		
		
		//Check for convergence
		converged = check_convergence(old_cent, centroids, no_of_dimensions, k_prev, k);

		free(old_cent);
	}
}

int check_convergence(centroid * c1, centroid* c2, int no_of_dimensions, int k_prev, int k){
	int j = 0;
	int b = 1;
	int i = 0;
	
	for(i = 0;  i < k; i++){
		for( j = 0; j < no_of_dimensions; j++){
			if(c1[i]->coordinates[j] != c2[i]->coordinates[j]){
				return 0;
			}
		}
	}
	return 1;
}

/* HASHTABLE FUNCTIONS */

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


long long int euclidean_hash_centroid(centroid* p, int no_dimensions, int window,  int no_of_functions, int table_size, int table_id, double* t, int *r, double** random_vectors_temp){
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