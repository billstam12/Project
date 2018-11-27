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
			centroids[i]->prev_count = -1;
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
		centroids[0]->prev_count = 0;
		centroids[0]->assigned_points = (point*)malloc(sizeof(point));
		centroids[0]->coordinates = (double*)malloc(no_of_dimensions*sizeof(double));
		for(i = 0; i< no_of_dimensions;i++){
			centroids[0]->coordinates[i] = data[id]->coordinates[i];
		}
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
			centroids[z+1]->prev_count = 0;
			centroids[z+1]->assigned_points = (point*)malloc(sizeof(point));
			centroids[z+1]->coordinates = (double*)malloc(no_of_dimensions*sizeof(double));
			for(i = 0; i< no_of_dimensions;i++){
				centroids[z+1]->coordinates[i] = data[id]->coordinates[i];
			}
			found_centers++;
		}
	}
	centroids[0]->id = 0;
	return centroids;
}
/* ASSIGNMENT FUNCTIONS */

void lloyds_assignment(centroid * centroids, point* data, int no_of_samples, int no_of_dimensions, int *k){
	int i,j,z;
	double min;
	int id, min2_id;

	for(i = 0; i < no_of_samples; i++){
		if(data[i]->centroid_id == -1){
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
				centroids[id]->prev_count = centroids[id]->count;
			}
			else{
				centroids[id]->assigned_points = (point*)change_mem(centroids[id]->assigned_points,(centroids[id]->prev_count*sizeof(point)),(centroids[id]->count)*sizeof(point));//Temporal
				centroids[id]->assigned_points[(centroids[id]->count)-1] = data[i];
				centroids[id]->prev_count = centroids[id]->count;
			}
			//printf("Point with ID:%ld was classified in %ld with distance %f\n",data[i]->id, id, dist);
		}
	}
}



void lsh_assignment(centroid* centroids, point* data, int no_of_samples, int no_of_dimensions, int metric, int no_of_functions, int L, int *k){
	int table_size = no_of_samples/12;
	int window = 30;
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
			for(j = 0; j < L; j++){
				data[i] -> g_functions[j] = (long long int*)malloc(no_of_functions*sizeof(long long int));
			}
		}	
		for(i = 0; i < L; i++){
			random_vectors[i] = init(&hts[i], data ,  no_of_samples,  no_of_dimensions, window,  no_of_functions, table_size, L,  i, metric, &t[i], &r[i]);
			//hashtable_print(hts[i]);
		}

		/* Now we will give the queries to our hashtables */
	
		int converged = 0;
		int * old_num = (int*)malloc(sizeof(int)*(*k));

		for(i=0; i < *k; i++){
			old_num[i] = centroids[i]->count;
		}

		while(converged == 0){ //Convergence will be done when the number of assigned points for each cluster has remained the same
			
			euclidean_lsh_query(data,  centroids, hts,  *k, no_of_samples, no_of_dimensions,  L,  radius,  window,  no_of_functions,  table_size,  t,  r,  random_vectors);
			radius = radius *2;
			int * new_num = (int*)malloc(sizeof(int)*(*k));
			for(i=0; i < *k; i++){
				new_num[i] = centroids[i]->count;
			}
			converged = check_lsh_convergence(new_num,old_num, *k);
			for(i=0; i < *k; i++){
				old_num[i] = new_num[i];
			}
		}
		int count = 0;
		int count2 =0;
		for(i = 0; i< no_of_samples; i++){
			if(data[i]->centroid_id == -1){
				count ++;
			}
			else{
				count2 ++;
			}
		}
		printf("count=%d\n", count);
		printf("count2=%d\n", count2);

		// Now assign each unassigned point with lloyds
		lloyds_assignment(centroids, data, no_of_samples, no_of_dimensions, k);
	}


	// Free data
	for(i = 0; i< L;i++){
		hashtable_free(&hts[i]);
		free(t[i]);
		free(r[i]);
		for(j = 0; j < no_of_functions; j++){
			free(random_vectors[i][j]);
		}
		free(random_vectors[i]);
	}
	for(i=0; i<no_of_samples;i++){
		free(data[i]->g_functions);
	}
	free(hts);
	free(t);
	free(r);
	free(random_vectors);
}

void euclidean_lsh_query(point* data, centroid* query_data, hashtable* hts, int no_queries, int no_samples, int no_of_dimensions, int L, double radius, int window, int no_of_functions, int table_size, double** t, int** r, double*** random_vectors){
	int i, j, z;
	for(i=0;i< no_queries;i++){
		query_data[i]->g_functions = (long long int**)malloc(sizeof(long long int)* L);
		for(j = 0; j < L; j++){
			query_data[i] -> g_functions[j] = (long long int*)malloc(no_of_functions*sizeof(long long int));
		}
	}	
	for(j = 0; j < no_queries; j++){
		/* L S H */

		for(i = 0; i < L; i++){
			
			//printf("TABLE NO: %d\n", i);
			long long int index = euclidean_hash_centroid(query_data[j], no_of_dimensions,  window,   no_of_functions,  table_size, i, t[i], r[i],  random_vectors[i]);			
			point tmp;
			tmp = hts[i]->buckets[index].first;
 			
			while(tmp != NULL){
				int g_func;
				g_func = compare_gfuncs(tmp->g_functions[i], query_data[j]->g_functions[i], no_of_functions);
				if(g_func == 1){
					double distance = euclidean_distance(tmp->coordinates, query_data[j]->coordinates, no_of_dimensions);
					if(distance < radius && distance!=0){
						if(data[tmp->id]->centroid_id == -1){\
							query_data[j]->count++;
							data[tmp->id]->centroid_id = j;
							data[tmp->id]->centroid2_id = j;
							query_data[j]->prev_count = query_data[j]->count;
						}
						else{
							double d1 = euclidean_distance(data[tmp->id]->coordinates, query_data[j]->coordinates, no_of_dimensions);
							double d2 = euclidean_distance(data[tmp->id]->coordinates, query_data[data[tmp->id]->centroid_id]->coordinates, no_of_dimensions);
							if( d1 > d2 ){
								//Add to new cluster	
								query_data[j]->count++;
								query_data[j]->prev_count = query_data[j]->count;
								//Remove from the old
								query_data[data[tmp->id]->centroid_id]->count--;
								data[tmp->id]->centroid_id = j;
								data[tmp->id]->centroid2_id = data[tmp->id]->centroid_id;
							}
						}
					}
				}
				tmp = tmp -> next;
			}
			
			free(tmp);
		}
	}
	for(i=0;i< no_queries ;i++){
		if(query_data[i]->count!=0){
			query_data[i]->assigned_points = (point*)malloc((query_data[i]->count* sizeof(point)));
			int jj = 0; //position in table
			for(j = 0; j < no_samples; j++){
				if(data[j]->centroid_id == i){
					query_data[i]->assigned_points[jj] = data[j];
					jj++;
				}
			}
		}
	}
}



void *change_mem(void *ptr, size_t originalLength, size_t newLength){

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
			printf("Found empty centroid\n");
		}
		else{

			for(j = 0; j < no_of_dimensions; j++){
				centroids[i]->coordinates[j] = 0;
				for(z = 0; z < centroids[i]->count; z++){
					centroids[i]->coordinates[j] += centroids[i]->assigned_points[z]->coordinates[j];
				}
				centroids[i]->coordinates[j] /= centroids[i]->count;
			}
		}
	}
}

void pam(centroid *centroids, point* data, int no_of_dimensions, int *k){
// Compute cluster cost function and find new centroid
	int i,j,z ;
	for(i = 0; i < (*k); i++){
		centroids[i]->dist = 0.0;
		for(j = 0; j < centroids[i]->count; j++){
			centroids[i]->dist += centroids[i]->assigned_points[j]->dist;
		}		
	}
	for(i = 0; i < (*k); i++){
		double min_cost = centroids[i]->dist;
		long int min_id = 0;
		if(centroids[i]->count==0){
			for(j=i;j<(*k);j++){
				centroids[j] = centroids[j+1];
			}
			i -= 1;
			*k -= 1;
			printf("PAM found empty centroid\n");
		}
		else{
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
			//printf("%ld %f\n",min_id, min_cost);
			for( j = 0; j< no_of_dimensions; j++){
				centroids[i]->coordinates[j] = centroids[i]->assigned_points[min_id]->coordinates[j];
			}
		}
	}
}

/* K-MEANS */
void compute_cluster(centroid * centroids, point * data, int no_of_samples, int no_of_dimensions, int* k, int assignment, int metric,int no_of_functions, int L, int type){
	int i,j,z;
	// Check where each point falls
	for(i = 0; i < no_of_samples; i++){
		data[i]->centroid_id = -1;
		data[i]->centroid2_id = -1;
	}

	if(assignment == 0){
		lloyds_assignment(centroids, data,  no_of_samples,  no_of_dimensions, k);

	}
	else if(assignment == 1){
		lsh_assignment(centroids, data, no_of_samples, no_of_dimensions, metric, no_of_functions,  L, k);
	}
	
	if(type == 0){
		basic_update(centroids, data, no_of_dimensions, k);
	}
	else{
		pam(centroids, data, no_of_dimensions, k);
	}
	
	
}



void kmeans(centroid* centroids,point* data, int no_of_samples, int no_of_dimensions,int k, int assignment, int metric, int no_of_functions, int L, int type){
	int converged = 0;
	int i,j;
	while(converged == 0){
		// Save the centroids for later
		
		for( i = 0; i< k;i++){
			free(centroids[i]->assigned_points);
			centroids[i]->assigned_points = (point*)malloc(sizeof(point));
			centroids[i] -> count = 0;
			centroids[i] -> prev_count = 0;
		}
		
		int k_prev = k;
		centroid * old_cent;
		old_cent = (centroid*)malloc(sizeof(centroid)*k);
		for( i = 0; i< k; i++){
			old_cent[i] = (centroid)malloc(sizeof(struct centroid));
			copy_centroid(&(old_cent[i]), centroids[i], no_of_dimensions );
		}
		compute_cluster(centroids, data, no_of_samples, no_of_dimensions, &k, assignment, metric,  no_of_functions,  L, type);
		/*
		for(i = 0; i < k; i++){
			print_coordinates_cent(old_cent[i], no_of_dimensions);
			print_coordinates_cent(centroids[i], no_of_dimensions);		
			printf("\n");
		}*/
		printf("TOTAL SAMPLES = %d\n", no_of_samples);
		for(i=0; i < k; i++){
			printf("%d\n", centroids[i]->count);
		}
		//Check for convergence
		converged = check_convergence(old_cent, centroids, no_of_dimensions, k_prev, k);

		free(old_cent);
		//converged=-1;
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


void shilouette_evaluation(centroid* centroids,point* data,int no_of_samples,int no_of_dimensions){
	int i,j;
	for(i = 0 ; i < no_of_samples; i++){
		int a_id = data[i]->centroid_id;
		int b_id = data[i]->centroid2_id;

		long double a = 0;
		long double b = 0;
		for(j = 0; j < centroids[a_id]->count; j++){
			a += euclidean_distance(data[i]->coordinates, centroids[a_id]->assigned_points[j]->coordinates, no_of_dimensions);
		}
		for(j = 0; j < centroids[b_id]->count; j++){
			b += euclidean_distance(data[i]->coordinates, centroids[b_id]->assigned_points[j]->coordinates, no_of_dimensions);
		}
		double max;
		if(a > b){
			max = a;
		}
		else{ max = b;}
		data[i]->silhouette= ((b-a)/max);
	}
}



int check_lsh_convergence(int *new_num,int *old_num, int k){
	int i;
	int converged = 1;
	for(i = 0; i < k; i++){
		converged = (new_num[i]==old_num[i]);
	}
	return converged;
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
			index = euclidean_hash(data[i], no_dimensions, window, no_of_functions, table_size, table_id, t_temp, r_temp, random_vectors_temp);
			/* initialize g functions here not in the hash */
			hashtable_insert(ht, data[i], index);

		}
		else{
			index = cosine_hash(data[i], no_dimensions, no_of_functions, table_size,  random_vectors_temp);
			/* initialize g functions here not in the hash */
			hashtable_insert(ht, data[i], index);
		}
	}
	return random_vectors_temp;

}


int cosine_hash(point p, int no_dimensions, int no_of_functions, int table_size, double** random_vectors_temp){
	int index = 0;
	int i, j;
	char * hashes;
	hashes = (char *)malloc(no_of_functions * (sizeof( char )+1));
	long double first_hash = 0.0;
	for(i = 0; i < no_of_functions; i++){
		for(j = 0; j < no_dimensions; j++){
			first_hash += ((p) -> coordinates[j]) * random_vectors_temp[i][j];
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

long long int euclidean_hash(point p, int no_dimensions, int window,  int no_of_functions, int table_size, int table_id, double* t, int *r, double** random_vectors_temp){
	long long int index = 0;
	int i, j;
	long long int * first_hashes;
	first_hashes = (long long int *)malloc(no_of_functions * sizeof(long long int ));
	long double first_hash = 0.0;

	for(i = 0; i < no_of_functions; i++){
		for(j = 0; j < no_dimensions; j++){
			first_hash += ((((p) -> coordinates[j] * random_vectors_temp[i][j]) + t[i])/window);
		}
		
		first_hashes[i] = (long long int)floor(first_hash);
	}
	for(i = 0 ; i< no_of_functions; i++){
		(p) -> g_functions[table_id][i] = first_hashes[i];
	}
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


long long int euclidean_hash_centroid(centroid p, int no_dimensions, int window,  int no_of_functions, int table_size, int table_id, double* t, int *r, double** random_vectors_temp){
	long long int index = 0;
	int i, j;
	long long int * first_hashes;
	first_hashes = (long long int *)malloc(no_of_functions * sizeof(long long int ));
	long double first_hash = 0.0;
	for(i = 0; i < no_of_functions; i++){
		for(j = 0; j < no_dimensions; j++){
			first_hash += ((((p) -> coordinates[j] * random_vectors_temp[i][j]) + t[i])/window);
		}
		
		first_hashes[i] = (long long int)floor(first_hash);
	
	}
	for(i = 0 ; i< no_of_functions; i++){
		(p) -> g_functions[table_id][i] = first_hashes[i];
	}
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
