#include "hashtable.h"

/* LINKED LIST */
void print_list(node_t * head) {
    node_t * current = head;

    while (current != NULL) {
        printf("%d ", current->val);
        current = current->next;
    }
    printf("\n");
}

int not_in_list(node_t * head, int val) {
    node_t * current = head;

    while (current != NULL) {
        if(current->val==val){
        	return 0;
        }
        current = current->next;
    }
    return 1;
}

void push(node_t * head, int val) {
    node_t * current = head;
    while (current->next != NULL) {
        current = current->next;
    }

    /* now we can add a new variable */
    current->next = (node_t*)malloc(sizeof(node_t));
    current->next->val = val;
    current->next->next = NULL;
}


void free_list(node_t* head)
{
   node_t* tmp;

   while (head != NULL)
    {
       tmp = head;
       head = head->next;
       free(tmp);
    }

}

/* HASHTABLE */

void hashtable_init(hashtable *ht, int size){
    (*ht) = (struct hashtable *)malloc(sizeof(hashtable));
    (*ht) -> size = size;
    (*ht) -> cnt = 0;
    (*ht) -> buckets = (bucket *)malloc(size * sizeof(bucket));
    for(int i =0; i < size; i++){
        (*ht) -> buckets[i] = (bucket)malloc(sizeof(struct bucket));
        (*ht) -> buckets[i] -> first = NULL;
    }
}

int hashtable_size(hashtable ht){
    return ht->cnt;
}

void bucket_insert(bucket *b, user u){
    if((*b)->first == NULL){
        (*b)->first = u;
        (*b)->first->next = NULL;
        (*b)->last = u;
    }
    else{
        (*b)->last->next = u;
        (*b)->last = u;
        (*b)->last->next =  NULL;
    }
}

void hashtable_print(hashtable ht){
    int i;
    for(i = 0; i< ht->size; i++){
        printf("BUCKET NO: %d\n", i);
        user tmp;
        tmp = ht->buckets[i]->first;
        while(tmp != NULL){
            printf("%d ", tmp->id);
            tmp = tmp -> next;
        }
        printf("\n");
    }
}

void hashtable_insert(hashtable *ht, user u, long long int index){
    bucket_insert(&((*ht)->buckets[index]), u);
}

void hashtable_free(hashtable *ht, int size){
    for(int i = 0; i < size; i++){
        free((*ht)->buckets[i]);
        (*ht)->buckets[i]->first=NULL;
        (*ht)->buckets[i]->last=NULL;
    }
    free((*ht)->buckets);
    (*ht)->buckets = NULL;
    (*ht)->size = 0;
    (*ht)->cnt = 0;
    free(*ht);
    (*ht) = NULL;
}

void bucket_print(hashtable ht, int index){
    printf("BUCKET NO: %d\n", index);
    user tmp;
    tmp = ht->buckets[index]->first;
    while(tmp != NULL){
        printf("%d ", tmp->id);
        tmp = tmp -> next;
    }
    printf("\n");
}