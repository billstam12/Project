
all:cluster 

cluster: cluster.o functions.o hashtable.o
	gcc  -o  cluster  cluster.o functions.o hashtable.o  -lm

cluster.o: cluster.c 
	gcc -c cluster.c 

functions.o: functions.c 
	gcc -c functions.c

hashtable.o: hashtable.c 
	gcc -c hashtable.c 

clean:
	-rm cluster.o functions.o hashtable.o  cluster 
	-rm  -r output/*