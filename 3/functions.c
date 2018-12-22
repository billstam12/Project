#include "functions.h"

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


tweet * parse_data(FILE* f, int * no_of_samples){
	tweet * tweets;

	rewind(f);
	char delimit[]="\t";

	char * line = NULL;
	size_t len = 0;
	ssize_t read = getline(&line, &len, f);
	
	/* First count the number of lines of the document == number of tweets */

	int samps = 0;
	while((read = getline(&line, &len, f)) !=  -1){
		samps ++;
	}

	*no_of_samples = samps;
	rewind(f);
	tweets = (tweet*)malloc(samps*sizeof(tweet));

	//Get no of dimensions of each twitter == number of words
	/* Each tweet has to ids, its users and its own and after that it has a list 
	of words that belong in a dictionary, that make up its sentiment score.
	For that we will first calculate the number of words in each sample.
	Meanwhile count the unique number of users, for later potential use. */
	
	int * no_of_words = (int*)malloc(sizeof(int)*samps);
	int * all_user_ids = (int*)malloc(sizeof(int)*samps);
	int * all_tweet_ids = (int*)malloc(sizeof(int)*samps);

	int i = 0;
	int j;
	while((read = getline(&line, &len, f)) !=  -1 ){
		tweet new_tweet;
		new_tweet = malloc(sizeof(struct tweet));
				
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

	unique unique_users;
	unique_users = get_unique(all_user_ids, all_tweet_ids, samps);
	// Now for each user we will get the array of his tweets and the number of them 
			

	rewind(f);
	
	return tweets;
}
