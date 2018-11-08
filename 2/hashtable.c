#include "hashtable.h"

void hashtable_init(hashtable *ht, int size){
	(*ht) = (struct hashtable *)malloc(sizeof(hashtable));
	(*ht) -> size = size;
	(*ht) -> cnt = 0;
	(*ht) -> buckets = (struct bucket *)malloc(size * sizeof(struct bucket));
}

int hashtable_size(hashtable ht){
	return ht->cnt;
}

void bucket_insert(bucket b, point p){
	if((b)->first == NULL){
		(b)->first = p;
		(b)->first->next = NULL;
		(b)->last = (b)->first;
	}
	else{
		(b)->last->next = p;
		(b)->last = p;
		(b)->last->next =  NULL;
	}
}

void hashtable_print(hashtable ht){
	int i;
	for(i = 0; i< ht->size; i++){
		printf("BUCKET NO: %d\n", i);
		point tmp;
		tmp = ht->buckets[i].first;
		while(tmp != NULL){
			printf("%ld ", tmp->id);
			tmp = tmp -> next;
		}
		printf("\n");
	}
}

void hashtable_insert(hashtable *ht, point p, long long int index){
	bucket_insert(&((*ht)->buckets[index]), p);
}

void hashtable_free(hashtable *ht){
	free((*ht)->buckets);
	(*ht)->buckets = NULL;
	(*ht)->size = 0;
	free(*ht);
	(*ht) = NULL;
}

void bucket_print(hashtable ht, int index){
	printf("BUCKET NO: %d\n", index);
	point tmp;
	tmp = ht->buckets[index].first;
	while(tmp != NULL){
		printf("%ld ", tmp->id);
		tmp = tmp -> next;
	}
	printf("\n");
}

void nn_list_insert(nn_list * nn, long int id, double dist, int table_found){
	nearest_neighbor n = (struct nearest_neighbor*)malloc(sizeof(struct nearest_neighbor));
	n->id = id;
	n->distance = dist;
	n->table = table_found;
	if((*nn)->first == NULL){
		(*nn)->first = n;
		(*nn)->first->next = NULL;
		(*nn)->last = (*nn)->first;
	}
	else{
		(*nn)->last->next = n;
		(*nn)->last = n;
		(*nn)->last->next =  NULL;
	}
}