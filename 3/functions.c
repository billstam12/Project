#include "functions.h"

/* NEW FUNCTIONS */

unique get_unique(int *all_user_ids, int *all_tweet_ids, int no_of_samples){
	unique unique_array;
	int curr_number = all_user_ids[0];
	int i, j;

	unique_array.count = 1;
	for(i = 1; i < no_of_samples; i++){
		if(all_user_ids[i]!=curr_number){
			unique_array.count ++;
		}
		curr_number = all_user_ids[i];
	}

	unique_array.tweets_per_user = (int*)malloc(sizeof(int)*unique_array.count);
	for(i = 0; i < unique_array.count; i++){
		unique_array.tweets_per_user[i] = 1;
	}

	curr_number = all_user_ids[0];
	j = 0;
	for(i = 0; i < no_of_samples; i++){
		curr_number = all_user_ids[i];
		while(all_user_ids[i+1]==curr_number){
			unique_array.tweets_per_user[j] ++;
			i++;
			//curr_number = all_user_ids[i];
		}
		j++;
	}

	unique_array.array = (int*)malloc(sizeof(int)*unique_array.count);
	unique_array.array[0] = all_user_ids[0];
	curr_number = all_user_ids[0];
	j = 1;
	for(i = 1; i < no_of_samples; i++){
		if(all_user_ids[i]!=curr_number){
			unique_array.array[j] = all_user_ids[i];
			j++;
		}
		curr_number = all_user_ids[i];
	}
	return unique_array;
}

cluster * parse_data_2(tweet* tweets, int no_of_tweets, int no_of_clusters){
	FILE * f = fopen("dataset/clustered_data.csv","r");
	int no_of_entries;
	int i,j;
	cluster * clusters;
	clusters = (cluster*)malloc(sizeof(cluster)*no_of_clusters);
	
	//The file created from project 2 is of the following type
	//no_of_entries score1    score2    score3
	// # entries followed by space and then the tab seperated scores for each cluster
	char * line = NULL;
	size_t len = 0;
	/* First count the number of lines of the document == number of tweets */
	ssize_t read;
	j = 0;
	int l;
	while((read = getline(&line, &len, f)) !=  -1 ){
		no_of_entries = atoi(strtok(line," "));
		clusters[j] = (cluster)malloc(sizeof(struct cluster));
		clusters[j]->tweets = (tweet*)malloc(sizeof(tweet)*no_of_entries);
		clusters[j]->size = no_of_entries;
		l = 0;
		while((line = strtok(NULL,"\t")) != NULL){	
			for(i = 0; i < no_of_tweets; i++){
				if(tweets[i]->id == atoi(line)){
					clusters[j]->tweets[l] = (tweet)malloc(sizeof(struct tweet));
					clusters[j]->tweets[l] = tweets[i];	
					tweets[i]->cluster_id = j;
					break;
				}
			}
			l++;
		}
		j++;
	}
	return clusters;
}
list_of_coins * parse_coins_query(FILE* f, int *no_of_coins){
	int i, j;
	rewind(f);
	char delimit[]="\t";

	char * line = NULL;
	size_t len = 0;
	/* First count the number of lines of the document == number of tweets */
	ssize_t read;
	int samps = 0;
	while((read = getline(&line, &len, f)) !=  -1){
		samps ++;
	}

	*no_of_coins = samps;
	list_of_coins * coins;
	coins = (list_of_coins*)malloc(sizeof(list_of_coins)*samps);
	//count words of each entry
	rewind(f);
	i = 0;
	while((read = getline(&line, &len, f)) !=  -1 ){
		coins[i]= (list_of_coins)malloc(sizeof(struct list_of_coins));
		line = strtok(line,delimit);
		int words = 1;
		while((line = strtok(NULL,delimit)) != NULL){		
			words ++;
		}
		coins[i]->array = (char**)malloc(sizeof(char*)* words);
		coins[i]->no_of_elements = words;
		i++;
	}
	i = 0;
	rewind(f);
	while((read = getline(&line, &len, f)) !=  -1 ){
		j = 0;
		line = strtok(line,delimit);
		coins[i]->array[j] = (char*)malloc(sizeof(char)*(strlen(line)+1));
		strcpy(coins[i]->array[j], line);	
		j++;
		while((line = strtok(NULL,delimit)) != NULL){		
			coins[i]->array[j] = (char*)malloc(sizeof(char)*(strlen(line)+1));
			strcpy(coins[i]->array[j], line);
			j++;
		}
		i++;
	}
	
	return coins;
}

users_and_tweets  parse_data(FILE* f){
	tweet * tweets;
	users_and_tweets data = (users_and_tweets)malloc(sizeof(struct users_and_tweets));


	rewind(f);
	char delimit[]="\t";

	char * line = NULL;
	size_t len = 0;
	ssize_t read;
	
	/* First count the number of lines of the document == number of tweets */
	printf("Parsing tweets..\n");
	int samps = 0;
	while((read = getline(&line, &len, f)) !=  -1){
		samps ++;
	}

	rewind(f);
	tweets = (tweet*)malloc(samps*sizeof(tweet));

	//Get no of dimensions of each twitter == number of words
	/* Each tweet has to ids, its users and its own and after that it has a list 
	of words that belong in a dictionary, that make up its sentiment score.
	For that we will first calculate the number of words in each sample.
	Meanwhile count the unique number of users, for later potential use. */

	int   * all_user_ids, *all_tweet_ids;
	all_user_ids = (int*)malloc(sizeof(int)*samps);
	all_tweet_ids = (int*)malloc(sizeof(int)*samps);

	int i = 0;
	int j;
	while((read = getline(&line, &len, f)) !=  -1 ){
		tweet new_tweet;
		new_tweet = (tweet) malloc(sizeof(struct tweet));
		new_tweet->coin_list = NULL;
		new_tweet->score = 0.0;
		line = strtok(line,delimit);
		new_tweet-> user_id = atof(line);
		all_user_ids[i] = atof(line); //user_id
		line = strtok(NULL,delimit);
		all_tweet_ids[i] = atof(line);
		new_tweet-> id = atof(line);
		new_tweet->no_of_words = 0;
		while((line = strtok(NULL,delimit)) != NULL){		
			new_tweet->no_of_words ++;
		}
		//printf("%d\n",new_tweet->no_of_words);
		tweets[i] = new_tweet;
		i++;
	}
	//Now load the list of words of the tweet

	rewind(f);
	i = 0;
	while((read = getline(&line, &len, f)) !=  -1 ){
		line = strtok(line,delimit);
		line = strtok(NULL,delimit);
		j = 0;
		tweets[i]->words = (char**)malloc(sizeof(char*)*tweets[i]->no_of_words);
		while((line = strtok(NULL,delimit)) != NULL){
			tweets[i]->words[j] = (char*)malloc(sizeof(char)*(strlen(line)+1));
			strcpy(tweets[i]->words[j], line);
			//printf("%s\n",tweets[i]->words[j]);
			j++;
		}
		//printf("%d\n",new_tweet->no_of_words);
		i++;
	}
	/* Moreover for each tweet we need to find the cryptocurrency it corresponds too
	and for that we will take advantage of the coins_queries */
	printf("Parsing tweets, Done!\n");
		
	char *coins_dic = "dataset/coins_queries.csv";
	FILE* c = fopen(coins_dic, "r");
	list_of_coins * coins;
	int no_of_coins;
	int z;
	int k,l;
	coins = parse_coins_query(c, &no_of_coins); //coins contains lists of each Coin Data
	
	for(i = 0; i < samps; i++){
		for(j = 0; j < no_of_coins; j++){
			for(z = 0; z < coins[j]->no_of_elements; z++){
				for(k = 0; k < tweets[i]->no_of_words; k ++){
					if(strcmp(tweets[i]->words[k],coins[j]->array[z])==0){
						if(tweets[i]->coin_list==NULL){
							tweets[i]->coin_list = (node_t*)malloc(sizeof(node_t));
							tweets[i]->coin_list->next = NULL;
							tweets[i]->coin_list->val =j;
						}
						else if(not_in_list(tweets[i]->coin_list, j)){
							//printf("ID=%d\n",tweets[i]->id);
							//print_list(tweets[i]->coin_list);
							push(tweets[i]->coin_list, j);
						}
					}
				}
			}
		}
	}

	printf("Creating coins dictionary..\n");
	char **dict = init_coin_dict();
	printf("Creating coins dictionary, Done!\n");

	// Here we will read the tweets from Project 2's dataset
	FILE * f2;
	f2 = fopen("dataset/twitter_dataset_small.csv", "r");
	char delimit2[]="\t,";
	line = NULL;
	len = 0;
	int samps2 = 0;
	while((read = getline(&line, &len, f2)) !=  -1){
		samps2 ++;
	}
	rewind(f2);
	//Get no of dimensions
	int flag = 1;
	int dims = 0;
	i = 0;	
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
	rewind(f2);

	/*  First create an array of points so we can store 
		each of the points in the memory. Then, we will begin
		hashing each point from our table of point-pointers 
		to the hashtable.
		So initialize the array of point-pointers.
	*/
	tweet2 * data2;
	data2 = (tweet2 *)malloc(samps * sizeof(tweet2));

	i = 0;
	while((read = getline(&line, &len, f2)) !=  -1 ){
		tweet2 new_point;
		new_point = malloc(sizeof(struct tweet2));
		new_point->coordinates = (double*) malloc(dims*sizeof(double));

		j = 0;		
		line = strtok(line,delimit2);
		new_point-> id = atof(line);
		while((line = strtok(NULL,delimit2)) != NULL){		
			new_point -> coordinates[j] = atof(line);
			data2[i] = new_point;
			j++;
		}
		i++;
	}
	
	// Now give the coordinates of these tweets to our real tweets
	for(i = 0; i < samps; i++){
		tweets[i]->coordinates = (double*)malloc(sizeof(double)*dims);
		for(j = 0; j < samps; j++){
			if(tweets[i]->id == data2[j]->id){
				for(z = 0; z < dims; z++){
					tweets[i]->coordinates[z] = data2[j]->coordinates[z];
				}
				break;
			}
		}
	}
	//Free the temporary tweet2 data
	for(i = 0; i < samps2; i++){
		free(data2[i]->coordinates);
	}
	free(data2);
	
	
	// Now that we have the database we will calculate the score of each tweet
	calculate_scores(tweets, samps);

	printf("Creating Users Database..\n");
	// Now for each user we will get the array of his tweets and the number of them 
	unique unique_users;
	unique_users = get_unique(all_user_ids, all_tweet_ids, samps);
	free(all_user_ids);
	free(all_tweet_ids);

	user * users;
	users = (user*)malloc(sizeof(user)*unique_users.count);
	j = 0;
	for(i = 0; i < unique_users.count; i++){
		users[i] = (user)malloc(sizeof(struct user));
		users[i]->id = unique_users.array[i];
		users[i]->no_of_tweets = unique_users.tweets_per_user[i];
		users[i]->tweets = (tweet*)malloc(sizeof(tweet)*unique_users.tweets_per_user[i]);
		users[i]->score_vector = (score_array *)malloc(sizeof(score_array)*100);
		for(z = 0; z < 100; z ++){
			users[i]->score_vector[z] = (score_array )malloc(sizeof(struct score_array));
			users[i]->score_vector[z]->value = 0;
			users[i]->score_vector[z]->id = z;
		}
		int no_of_tweet = 0;
		int threshold = j + unique_users.tweets_per_user[i];
		while(j < threshold){
			users[i]->tweets[no_of_tweet] = tweets[j];
			node_t * current = tweets[j]->coin_list;
			// ADD THE SCORES OF EACH TYPE OF COIN TWEET
			while(current!=NULL){
				users[i]->score_vector[current->val]->value += tweets[j]->score;
				current = current->next;
			}
			no_of_tweet ++;
			j++;
		}
	}
	printf("Creating Users Database, Done!\n");
	int no_of_clusters = 50;
	printf("Reading Clustering data from Project 2...\n");
	cluster* clusters = parse_data_2(tweets, samps, no_of_clusters);
	printf("Reading Clustering data from Project 2, Done!\n");

	data->no_of_users = unique_users.count;
	data->no_of_tweets = samps;
	data->no_of_clusters = no_of_clusters;
	data->no_of_coins = no_of_coins;
	
	rewind(f);
	fclose(f);

	for(i = 0; i < no_of_coins; i++){
		for(j = 0; j < coins[i]->no_of_elements; j++){
			free(coins[i]->array[j]);
		}
		free(coins[i]->array);
		free(coins[i]);
	}
	free(coins);
	fclose(c);

	data->users = users;
	data->tweets = tweets;
	data->clusters = clusters;
	data->dict = dict;
	return data;
}


void calculate_scores(tweet * tweets, int no_of_samples){
	int i, j, z;
	// Read Dictionary
	FILE * f;
	f = fopen("dataset/vader_lexicon.csv", "r");
	rewind(f);
	char delimit[]="\t";

	char * line = NULL;
	size_t len = 0;
	/* First count the number of lines of the document == number of words */
	ssize_t read;
	int samps = 0;
	while((read = getline(&line, &len, f)) !=  -1){
		samps ++;
	}

	rewind(f);
	printf("Calculating Scores..\n");
	i = 0;
	while((read = getline(&line, &len, f)) !=  -1){
		line = strtok(line,delimit);
		double score = atof(strtok(NULL,"\n"));
		for(i = 0; i < no_of_samples; i++){
			for(j = 0; j < tweets[i]->no_of_words; j++){
				if(strcmp(tweets[i]->words[j],line)==0){
					tweets[i]->score += score;
				}
			}			
		}
		i++;
	}
	/*Normalize Scores based on given formula */
	for(i = 0; i < no_of_samples; i++){
		tweets[i]->score = (tweets[i]->score)/sqrt((pow(tweets[i]->score ,2) + 15));
		//printf("%f\n",tweets[i]->score);
	}
	
	printf("Calculating Scores, Done!\n");

	/*
	for(i = 0; i < samps; i++){
		free(dictionary[i]->word);
		free(dictionary[i]);
	}
	*/	

	//free(dictionary);
	fclose(f);
}


/* OLD FUNCTIONS */

double cosine_similarity(score_array * c1, score_array * c2, int no_dimensions){
	int i;
	double similarity = 0.0;
	double c1_dist = 0.0;
	double c2_dist = 0.0;
	for(i = 0; i < no_dimensions; i++){
		c1_dist += pow(c1[i]->value, 2);
		c2_dist += pow(c2[i]->value, 2);
	}
	c1_dist = sqrt(c1_dist); c2_dist = sqrt(c2_dist);

	for(i = 0; i < no_dimensions; i++){
		similarity += (c1[i]->value * c2[i]->value);
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

double euclidean_distance(score_array * c1, score_array * c2, int no_dimensions){
	int i;
	double dist = 0.0;
	for(i = 0; i < no_dimensions; i++){
		dist += pow((c1[i]->value - c2[i]->value), 2);
	}	
	return sqrt(dist);
}

// QUESTION A STARTS HERE

int cosine_hash(user u, int no_dimensions, int no_of_functions, int table_size, double** random_vectors_temp){
	int index = 0;
	int i, j;
	char * hashes;
	hashes = (char *)malloc(no_of_functions * (sizeof( char )+1));
	long double first_hash = 0.0;
	for(i = 0; i < no_of_functions; i++){
		for(j = 0; j < no_dimensions; j++){
			first_hash += ((u) -> score_vector[j]->value) * random_vectors_temp[i][j];
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

double** init(hashtable* ht, user *data , int no_samples, int no_dimensions, int window, int no_of_functions, int table_size,int L, int table_id, double ** t,  int ** r){
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
	
	int index;
	for(i = 0; i < no_samples; i++){
		index = cosine_hash(data[i], no_dimensions, no_of_functions, table_size,  random_vectors_temp);
		/* initialize g functions here not in the hash */
		hashtable_insert(ht, data[i], index);
		
	}
	return random_vectors_temp;

}

int** cosine_lsh_implementation(user* users, int no_of_users, int p, int type, char ** dict, char* output){
	//Initialize the hashtable variables
	int L = 1;
	int window = 30;
	int no_of_dimensions = 100;
	int no_of_functions = 4;
	int i, j, z;
	FILE * o = fopen(output, "a");
	fprintf(o,"Cosine LSH\n");
	fclose(o);
	clock_t start, end;
	double cpu_time_used;
 	start = clock();
	double ** t; // t array
	t = (double**)malloc(L*sizeof(double*));
	int ** r; // r array
	r = (int**)malloc(L*sizeof(int*));
	double *** random_vectors; // random hyperplane vectors
	random_vectors = (double***)malloc(L*sizeof(double**));
	hashtable * hts;
	hts = (hashtable*)malloc(L * sizeof(hashtable));
	int table_size = pow(2,no_of_functions);


	// Initialize hashtables
	for(i = 0; i < L; i ++){
		hashtable_init(&hts[i], table_size);
	}
	for(i = 0; i < L; i++){
		random_vectors[i] = init(&hts[i], users , no_of_users,  no_of_dimensions, window,  no_of_functions, table_size, L,  i, &t[i], &r[i]);
		//hashtable_print(hts[i]);
		//Query the table again to find neighbors
		if(type == 0){
			cosine_lsh_query(users, hts[i],  no_of_users , no_of_dimensions,  L,  p,  window,  no_of_functions,  table_size,  t,  r,  random_vectors[i], type, dict, output);
		}
		else{
			return cosine_lsh_query(users, hts[i],  no_of_users , no_of_dimensions,  L,  p,  window,  no_of_functions,  table_size,  t,  r,  random_vectors[i], type, dict, output);
		}
	}
	for(i = 0; i < L; i++){
		hashtable_free(&hts[i], table_size);
		free(t[i]);
		free(r[i]);
		for(j = 0; j < no_of_functions; j++){
			free(random_vectors[i][j]);
		}
		free(random_vectors[i]);
	}
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	o = fopen(output, "a");
	fprintf(o, "Time:%f\n", cpu_time_used);
	fclose(o);

}


int** cosine_lsh_query(user* data, hashtable ht, int no_of_users,  int no_of_dimensions, int L, int no_of_neighbors, int window, int no_of_functions, int table_size, double** t, int** r, double** random_vectors, int type, char** dict, char* output){
	int i, j, z;	
	int ** coins = (int**)malloc(sizeof(int*)*no_of_users);
	FILE * o = fopen(output, "a");
	for ( i = 0; i < no_of_users; i++){
		long long int index = cosine_hash(data[i], no_of_dimensions, no_of_functions, table_size,  random_vectors);			
		user tmp=NULL;
		tmp = ht->buckets[index]->first;
		j = 0;
		user * nearest_neighbors = (user*)malloc(sizeof(user)*no_of_neighbors); //initialize nn list	
		while(tmp != NULL){
			double distance = cosine_similarity(tmp->score_vector, data[i]->score_vector, no_of_dimensions);

			if((distance!=0)){
				if(j < no_of_neighbors){
					nearest_neighbors[j] = tmp;
				}
				else{
					//Find maximum distance and replace
					int max_id = 0;
					double max_distance = 0.0;
					user max = nearest_neighbors[0];
					for(z = 1; z < no_of_neighbors; z++){
						double d1 = cosine_similarity(nearest_neighbors[z]->score_vector, data[i]->score_vector, no_of_dimensions);
						double d2 = cosine_similarity(max->score_vector, data[i]->score_vector, no_of_dimensions);
						
						if(d1 > d2){
							max_distance = d2;
							max_id = z;
						}
					}
					//If the maximum distance is maximum than the one you are checking, then change it
					if(distance <= max_distance){
						nearest_neighbors[max_id] = tmp;
					}
				}
				j++;

			}
			tmp = tmp -> next;
		}
		int real_no = no_of_neighbors;
		if(j < no_of_neighbors){
			real_no = j;
		}
		/* Now make predictions for the other coins based on the neighbors */
	
		int l = 0, k;
		for(k = 0; k < 100; k++){
			if(data[i]->score_vector[k]->value!=0){
				data[i]->score_vector[k]->value = -10;//set it to something real big so we won't have problems
			}
		}

		for(z = 0; z < real_no; z++){
			for(k = 0; k < 100; k++){
				if(data[i]->score_vector[k]->value==0){
					data[i]->score_vector[k]->value += (nearest_neighbors[z]->score_vector[k]->value)/no_of_neighbors;
				}
			}
		}
		// Now from the new coins find the 5 best
		qsort(data[i]->score_vector, 100, sizeof(data[i]->score_vector[0]), cmp);
	
		//printf("\n");
		if(type == 1){
			coins[i] = (int*)malloc(sizeof(int)*8);  //return more than 2, just in case we don't return stuff that the user has already seen
			for(j = 0; j < 8; j ++){
				coins[i][j] = data[i]->score_vector[j]->id;
			}
		}
		else{
			fprintf(o ,"<%d> ", data[i]->id);
			for(k = 0; k < 5; k++){
				fprintf(o, "%s ",dict[data[i]->score_vector[k]->id]);	
			}
			fprintf(o, "\n");
		}
	}
	fclose(o);
	return coins;
}

int cmp(const void *a, const void *b){
	score_array *a1 = (score_array*)a;
	score_array *a2 = (score_array*)b;
	if((*a1)->value > (*a2)->value) return -1;
	else if((*a1)->value < (*a2)-> value) return 1;
	else return 0;
}


/* OLD FUNCTIONS, CLUSTERING */
void print_coordinates_cent(centroid p, int no_dimensions){
	int i;
	for( i = 0; i < no_dimensions; i++)
		printf("%f ", p->score_vector[i]->value);	
	printf("\n");
}

void copy_centroid(centroid* c1, centroid c2, int no_of_dimensions){
	int i;
	(*c1)->id = c2->id;
	(*c1)->count = c2->count;
	(*c1)->silhouette_of_cluster = c2->silhouette_of_cluster;
	(*c1)->score_vector = (score_array*)malloc(no_of_dimensions*sizeof(score_array));
	for(i = 0; i < no_of_dimensions; i++){
		(*c1)->score_vector[i] = (score_array)malloc(sizeof(struct score_array));
		(*c1)->score_vector[i]->value = c2->score_vector[i]->value;
		(*c1)->score_vector[i]->id = c2->score_vector[i]->id;
	}
	(*c1)->assigned_users = (user*)malloc((*c1)->count*sizeof(user));
	for(i = 0; i < (*c1)->count; i++){
		copy_user(&((*c1)->assigned_users[i]), c2->assigned_users[i], no_of_dimensions);
	}
}

void copy_user(user* p1, user p2, int no_of_dimensions){
	int i;
	(*p1)->id = p2->id;
	(*p1)->centroid_id = p2->centroid_id;
	(*p1)->centroid2_id = p2->centroid2_id;
	(*p1)->score_vector = (score_array*)malloc(no_of_dimensions*sizeof(score_array));
	(*p1)->no_of_tweets = p2->no_of_tweets;
	(*p1)->tweets = (tweet*)malloc(sizeof(tweet)*(*p1)->no_of_tweets);
	for(i = 0; i < no_of_dimensions; i++){
		(*p1)->score_vector[i] = (score_array)malloc(sizeof(struct score_array));
		(*p1)->score_vector[i]->value = p2->score_vector[i]->value;
		(*p1)->score_vector[i]->id = p2->score_vector[i]->id;		
	}
	for(i = 0; i < (*p1)->no_of_tweets; i++){
		(*p1)->tweets[i] = (p2)->tweets[i];
	}
}

void cosine_cluster_implementation(user *users2, int no_of_users, int no_of_dimensions, int k, int p, int type, char** dict, char* output){
	FILE* o = fopen(output, "a");
	fprintf(o,"Cosine Cluster\n");
	clock_t start, end;
	double cpu_time_used;
 	start = clock();
	centroid * centroids;
	centroids = init_centroids(k, users2, no_of_users, no_of_dimensions, type);
	kmeans(centroids, users2, no_of_users, no_of_dimensions, k, p );
	k_means_recommend(centroids, users2, no_of_users, no_of_dimensions, k, p, 0, dict , output);
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	fprintf(o, "Time:%f\n", cpu_time_used);
	fclose(o);

}

centroid * init_centroids(int k, user* data, int no_of_samples, int no_of_dimensions, int type){
	centroid *centroids = (centroid*)malloc(k*sizeof(centroid));
	int i,j,z;
	if(type==0){
		for(i = 0; i< k ; i++){
			centroid new_centroid;
			new_centroid = (centroid)malloc(sizeof(struct centroid));
			centroids[i] = new_centroid;
			centroids[i]->id = i;
			centroids[i]->count = 0;
			centroids[i]->silhouette_of_cluster = 0;
			centroids[i]->prev_count = -1;
			centroids[i]->assigned_users = (user*)malloc(sizeof(user));
			//centroids[i]->coordinates = data[rand()%no_of_samples]->coordinates;
			centroids[i]->score_vector = (score_array*)malloc(no_of_dimensions*sizeof(score_array));
			int random_index = rand()%no_of_samples;
			for(j = 0; j< no_of_dimensions;j++){
				centroids[i]->score_vector[j] = (score_array )malloc(sizeof(struct score_array));
				centroids[i]->score_vector[j]->value = data[random_index]->score_vector[j]->value;
				centroids[i]->score_vector[j]->id = data[random_index]->score_vector[j]->id;
			}
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
		centroids[0]->assigned_users = (user*)malloc(sizeof(user));
		centroids[0]->score_vector = (score_array*)malloc(no_of_dimensions*sizeof(score_array));
		for(i = 0; i< no_of_dimensions;i++){
			centroids[0]->score_vector[i] = (score_array )malloc(sizeof(struct score_array));
			centroids[0]->score_vector[i]->value = data[id-1]->score_vector[i]->value;
			centroids[0]->score_vector[i]->id = data[id-1]->score_vector[i]->id;
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
						data[j]-> dist += euclidean_distance(centroids[i]->score_vector, data[j]->score_vector, no_of_dimensions);
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
			centroids[z+1]->prev_count = -1;
			centroids[z+1]->assigned_users = (user*)malloc(sizeof(user));
			centroids[z+1]->score_vector = (score_array*)malloc(no_of_dimensions*sizeof(score_array));
			for(i = 0; i< no_of_dimensions;i++){
				centroids[z+1]->score_vector[i] = (score_array )malloc(sizeof(struct score_array));
				centroids[z+1]->score_vector[i]->value = data[id-1]->score_vector[i]->value;
				centroids[z+1]->score_vector[i]->id = data[id-1]->score_vector[i]->id;
			}
			found_centers++;
		}
	}
	centroids[0]->id = 0;
	return centroids;
}


void kmeans(centroid* centroids, user* data, int no_of_samples, int no_of_dimensions, int k, int no_of_neighbors){
	int i,j;
	clock_t start, end;
	double cpu_time_used;
 	start = clock();
 	int max_loops = 0;
 	int converged = 0;
	while(converged == 0 && max_loops < 20){
		max_loops ++;
		// Save the centroids for later
		
		for( i = 0; i< k;i++){
			free(centroids[i]->assigned_users);
			centroids[i]->assigned_users = NULL;
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
		compute_cluster(centroids, data, no_of_samples, no_of_dimensions, k);
		/*
		for(i = 0; i < k; i++){
			print_coordinates_cent(old_cent[i], no_of_dimensions);
			print_coordinates_cent(centroids[i], no_of_dimensions);		
			printf("\n");
		}*/
		//Check for convergence
		converged = check_convergence(old_cent, centroids, no_of_dimensions, k);

		free(old_cent);
		//converged=-1;
	}
	for(i = 0; i < k; i++){
		//printf("CLUSTER: %d SIZE: %d ",i, centroids[i]->count );
		//print_coordinates_cent(centroids[i], no_of_dimensions);
	}
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	//printf("Clustering Time:%f\n", cpu_time_used);
	shilouette_evaluation(centroids, data, no_of_samples, no_of_dimensions, k);

	/* Now we will recommend based on the clusters */

}


int check_convergence(centroid * c1, centroid* c2, int no_of_dimensions, int k){
	int j = 0;
	int b = 1;
	int i = 0;
	
	for(i = 0;  i < k; i++){
		for( j = 0; j < no_of_dimensions; j++){
			if(c1[i]->score_vector[j]->value != c2[i]->score_vector[j]->value){
				return 0;
			}
		}
	}
	return 1;
}

void compute_cluster(centroid * centroids, user * data, int no_of_samples, int no_of_dimensions, int k){
	int i,j,z;
	// Check where each point falls
	for(i = 0; i < no_of_samples; i++){
		data[i]->centroid_id = -1;
		data[i]->centroid2_id = -1;
	}
	/* Type of assignment, lloyds, lsh and hypercube*/
	lloyds_assignment(centroids, data,  no_of_samples,  no_of_dimensions, k);
	basic_update(centroids, data, no_of_dimensions, k);
}

void lloyds_assignment(centroid * centroids, user* data, int no_of_samples, int no_of_dimensions, int k){
	int i,j,z;
	double min;
	int id, min2_id;

	for(i = 0; i < no_of_samples; i++){
		if(data[i]->centroid_id == -1){
			min = euclidean_distance((centroids[0])->score_vector, (data[i])->score_vector, no_of_dimensions);
			id = (centroids[0])->id;
			min2_id = id;
			double dist;
			for(j = 1; j < k; j++){
				dist = euclidean_distance((centroids[j])->score_vector, (data[i])->score_vector, no_of_dimensions);
				if(dist < min){
					min2_id = id;
					min = dist;
					id = (centroids[j])-> id;
				}
			}
			(data[i])->centroid_id = id;
			(data[i])->centroid2_id = min2_id;

			(centroids[id])->count +=1;		
			//printf("Point with ID:%ld was classified in %ld with distance %f\n",data[i]->id, id, dist);
		}

	}
	for(i=0; i< k ;i++){
		if(centroids[i]->count!=0){
			centroids[i]->assigned_users = (user*)malloc((centroids[i]->count* sizeof(user)));
			int jj = 0; //position in table
			for(j = 0; j < no_of_samples; j++){
				if(data[j]->centroid_id == i){
					centroids[i]->assigned_users[jj] = data[j];
					jj++;
				}
			}
		}
	}

}


void basic_update(centroid *centroids, user* data, int no_of_dimensions, int k){
	// Compute cluster mean and find new centroid
	int i,j,z ;

	for(i = 0; i < k; i++){
		// If for some reason one cluster is empty, then remove it 
		if(centroids[i]->count==0){
			for(j = i; j < k; j++){
				centroids[j] = centroids[j+1];
			}
			i -= 1;
			k -= 1;
			fprintf(stderr, "ERROR! Update resulted in an empty centroid, please rerun with different parameters\n");
			exit(-1);
		}
		else{
			for(j = 0; j < no_of_dimensions; j++){
				centroids[i]->score_vector[j]->value = 0;
				for(z = 0; z < centroids[i]->count; z++){
					centroids[i]->score_vector[j]->value += centroids[i]->assigned_users[z]->score_vector[j]->value;
				}
				centroids[i]->score_vector[j]->value /= centroids[i]->count;
			}
		}
	}
}

void shilouette_evaluation(centroid* centroids, user* data, int no_of_samples, int no_of_dimensions, int k){
	int i,j;
	for(i = 0 ; i < no_of_samples; i++){
		int a_id = data[i]->centroid_id;
		int b_id = data[i]->centroid2_id;
		long double a = 0;
		long double b = 0;
		for(j = 0; j < centroids[a_id]->count; j++){
			a += euclidean_distance(data[i]->score_vector, centroids[a_id]->assigned_users[j]->score_vector, no_of_dimensions);
		}
		for(j = 0; j < centroids[b_id]->count; j++){
			b += euclidean_distance(data[i]->score_vector, centroids[b_id]->assigned_users[j]->score_vector, no_of_dimensions);
		}
		double max;
		if(a > b){
			max = a;
		}
		else{ max = b;}
		data[i]->silhouette= ((b-a)/max);
	}
	//centroids[i]->silhouette_of_cluster = 0.0;
	for(i = 0 ; i < k; i++){
		for(j = 0; j < centroids[i]->count; j++){
			centroids[i]->silhouette_of_cluster+=centroids[i]->assigned_users[j]->silhouette;
		}
		centroids[i]->silhouette_of_cluster /= centroids[i]->count;
	}
}

int ** k_means_recommend(centroid* centroids, user* data, int no_of_samples, int no_of_dimensions, int k, int no_of_neighbors, int type, char** dict, char * output){
	int i, j, z;
	int ** coins = (int**)malloc(sizeof(int*)*no_of_samples);
	FILE * o;
	o = fopen(output, "a");
	//2d array containing the base mean scores of each cluster
	score_array** clusters_mean = (score_array**)malloc(sizeof(score_array*)*k);
	for(i = 0; i < k; i++){
		clusters_mean[i] = (score_array*)malloc(sizeof(score_array)*100);
		for(j = 0; j < 100; j++){
			clusters_mean[i][j] = (score_array)malloc(sizeof(struct score_array));
			clusters_mean[i][j]->value = 0;
		}
		for(j = 0; j < centroids[i]->count; j++){
			for(z = 0; z < 100; z++){
				clusters_mean[i][z]->value += (centroids[i]->assigned_users[j]->score_vector[z]->value/centroids[i]->count);
				clusters_mean[i][z]->id = z;
			}
		}
		// Now calculate for each user the mean of his unspoken cryptos
		for(j = 0; j < centroids[i]->count; j++){
			int real_no = no_of_neighbors;
			if(j < no_of_neighbors){
				real_no = j;
			}
			/* Now make predictions for the other coins based on the neighbors */
			for(z = 0; z < 100; z++){
				if(centroids[i]->assigned_users[j]->score_vector[z]->value!=0){
					centroids[i]->assigned_users[j]->score_vector[z]->value = -10;//set it to something real big so we won't have problems
				}
				else{
					centroids[i]->assigned_users[j]->score_vector[z]->value = clusters_mean[i][z]->value;
				}
			}
			
			// Now from the new coins find the 5 best
			qsort(centroids[i]->assigned_users[j]->score_vector, 100, sizeof(centroids[i]->assigned_users[j]->score_vector), cmp);
			
			//printf("\n");
			if(type == 1){
				coins[i] = (int*)malloc(sizeof(int)*8);  //return more than 2, just in case we don't return stuff that the user has already seen
				int z;
				for(z = 0; z < 8; z ++){
					coins[i][z] = centroids[i]->assigned_users[j]->score_vector[z]->id;
				}
			}
			else{
				fprintf(o ,"<%d> ", centroids[i]->assigned_users[j]->id);
				int z;			
				for(z = 0; z < 5; z++){
					fprintf(o ,"%s ",dict[centroids[i]->assigned_users[j]->score_vector[z]->id]);	
				}
				fprintf(o ,"\n");
			
			}
		}
	}
	//free centroids
	for(i = 0; i < k; i ++){
		for(j = 0; j < no_of_dimensions; j++){
			free(centroids[i]->score_vector[j]);
		}
		free(centroids[i]);
	}
	free(centroids);
	fclose(o);
	return coins;
}

//QUESTION A ENDS HERE

// QUESTION B STARTS HERE

void clustering_lsh_implementation(user* users, cluster* clusters, int no_of_users, int no_of_clusters, int p, char** dict, char *output){
	/* Create the c(j) vectors */
	int i,j,z;
	FILE *o = fopen(output, "a");
	fprintf(o,"Clustering LSH\n");
	clock_t start, end;
	double cpu_time_used;
 	start = clock();
	user * c = (user*)malloc(sizeof(user)*no_of_clusters);

	for(i = 0; i < no_of_clusters; i++){
		c[i] = (user)malloc(sizeof(struct user)*100);
		c[i]->score_vector = (score_array *)malloc(sizeof(score_array)*100);
		for(z = 0; z < 100; z ++){
			c[i]->score_vector[z] = (score_array )malloc(sizeof(struct score_array));
			c[i]->score_vector[z]->value = 0;
			c[i]->score_vector[z]->id = z;
		}
	}
	for(j = 0; j < 100; j++){
		for(i = 0; i < no_of_clusters; i++){
			for(z = 0; z < clusters[i]->size; z++){
				c[i]->score_vector[j]->value += clusters[i]->tweets[z]->coordinates[j];
			}

		}
	}

	// C is a no_of_clusters rows and 100 columns long array
	int** coins = cosine_lsh_implementation(c, no_of_clusters, p, 1, dict, output); // Now we will have to correspond the return list with the users
	for(i = 0; i < no_of_users;  i++){
		fprintf(o ,"<%d> ", users[i]->id);
		for(j = 0; j < no_of_clusters; j++){
			int count = 0; //count till 2
			for(z = 0; z < 8; z++){
				if(count < 2){
					int k;
					int flag = 0;
					for(k = 0; k < users[i]->no_of_tweets; k++){
						if(users[i]->score_vector[coins[j][z]]->value != 0){ // filter here
							flag = 1;
							break;
						}

					}
					if(flag == 0){
						fprintf(o, "%s ", dict[coins[j][z]]);
						count++;
					}
				}
				else break;
			}
			if(count >= 2){
				break;
			}
		}
		fprintf(o, "\n");
	}

	//free data
	for(i = 0; i < no_of_clusters; i++){
		free(coins[i]);
	}
	for(i = 0; i < no_of_clusters; i++){
		for(j = 0; j < 100; j++){
			free(c[i]->score_vector[j]);		
		}
		free(c[i]->score_vector);
		free(c[i]);
	}
	free(c);
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	fprintf(o, "Time:%f\n", cpu_time_used);
	fclose(o);

}



void clustering_cluster_implementation(user* users, cluster* clusters, int no_of_users, int no_of_clusters, int p, char** dict, char *output){
	int i,j,z;
	user * c = (user*)malloc(sizeof(user)*no_of_clusters);
	FILE *o = fopen(output, "a");
	fprintf(o,"Clustering LSH\n");
	clock_t start, end;
	double cpu_time_used;
 	start = clock();
	for(i = 0; i < no_of_clusters; i++){
		c[i] = (user)malloc(sizeof(struct user)*100);
		c[i]->score_vector = (score_array *)malloc(sizeof(score_array)*100);
		for(z = 0; z < 100; z ++){
			c[i]->score_vector[z] = (score_array )malloc(sizeof(struct score_array));
			c[i]->score_vector[z]->value = 0;
			c[i]->score_vector[z]->id = z;
		}
	}
	for(j = 0; j < 100; j++){
		for(i = 0; i < no_of_clusters; i++){
			for(z = 0; z < clusters[i]->size; z++){
				c[i]->score_vector[j]->value += clusters[i]->tweets[z]->coordinates[j];
			}
		}
	}

	// Looyds basic assignment with basic update, and k-means++ initialization
	int k = 5; 
	centroid * centroids;
	int type = 0;
	centroids = init_centroids(k, c, no_of_clusters, 100, type);
	kmeans(centroids, c, no_of_clusters, 100, k, p); 
	int** coins = k_means_recommend(centroids, c,  no_of_clusters,  100,  k,  p, 1, dict, output);
	for(i = 0; i < no_of_users;  i++){
		fprintf(o ,"<%d> ", users[i]->id);
		for(j = 0; j < no_of_clusters; j++){
			int count = 0; //count till 2
			for(z = 0; z < 8; z++){
				if(count < 2){
					int k;
					int flag = 0;
					for(k = 0; k < users[i]->no_of_tweets; k++){
						if(users[i]->score_vector[coins[j][z]]->value != 0){ // filter here
							flag = 1;
							break;
						}

					}
					if(flag == 0){
						fprintf(o ,"%s ", dict[coins[j][z]]);
						count++;
					}
				}
				else break;
			}
			if(count >= 2){
				break;
			}
		}
		fprintf(o, "\n");
	}
	//free data
	for(i = 0; i < no_of_clusters; i++){
		free(coins[i]);
	}
	for(i = 0; i < no_of_clusters; i++){
		for(j = 0; j < 100; j++){
			free(c[i]->score_vector[j]);
		}
		free(c[i]->score_vector);
		free(c[i]);
	}
	free(c);
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	fprintf(o, "Time:%f\n", cpu_time_used);
	fclose(o);
}

// Function for converting ids to coins

char ** init_coin_dict(){
	int i, j;
	FILE * f = fopen("dataset/coin_clean_dict.txt", "r");
	char * line = NULL;
	size_t len = 0;
	/* First count the number of lines of the document == number of tweets */
	ssize_t read;
	int samps = 0;
	while((read = getline(&line, &len, f)) !=  -1){
		samps ++;
	}
	rewind(f);
	char ** dict = (char**)malloc(sizeof(char*)*samps);
	i = 0;
	while((read = getline(&line, &len, f)) !=  -1){
		len = strlen(line);
		if (len > 0 && line[len-1] == '\n')
    		line[len-1] = 0;
		dict[i] = (char*)malloc(sizeof(char)*strlen(line));
		strcpy(dict[i], line);
		i ++;
	}
	return dict;
}

void free_data(tweet * tweets, int no_of_samples, user* users, int no_of_users,  cluster* clusters, char ** dict, int no_of_coins){
	int i,j;
	for(i = 0; i < no_of_samples; i++){
	    free_list(tweets[i]->coin_list);
		for(j = 0; j < tweets[i]->no_of_words; j++){
			free(tweets[i]->words[j]);
		}
		free(tweets[i]->words);
		free(tweets[i]);
	}
	free(tweets);
	for(i = 0; i < no_of_users; i++){
		for(j = 0; j < users[i]->no_of_tweets; j++){
			free(users[i]->score_vector[j]);
		}
		free(users[i]->score_vector);
		free(users[i]);
	}
	free(users);
	free(clusters);
	for(i = 0; i < no_of_coins; i++){
		free(dict[i]);
	}
	free(dict);
}

user * copy_users(user* users, int no_of_users){
	int i;
	user * users2;
	users2 = (user*)malloc(sizeof(user)*no_of_users);
	for(i = 0; i < no_of_users; i++){
		users2[i] = (user)malloc(sizeof(struct user));
		copy_user(&users2[i], users[i], 100);
	}
	return users2;
}

void free_users(user * users, int no_of_users){
	int i,j;
	for(i = 0; i < no_of_users; i++){
		for(j = 0; j < users[i]->no_of_tweets; j++){
			free(users[i]->score_vector[j]);
		}
		free(users[i]->score_vector);
		free(users[i]);
	}
	free(users);
}