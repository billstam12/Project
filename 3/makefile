
all:recommendation 

recommendation: recommendation.o functions.o hashtable.o
	gcc  -o  recommendation  recommendation.o functions.o hashtable.o  -lm

recommendation.o: recommendation.c 
	gcc -c recommendation.c 

functions.o: functions.c 
	gcc -c functions.c

hashtable.o: hashtable.c 
	gcc -c hashtable.c 

clean:
	-rm recommendation.o functions.o hashtable.o  recommendation 
	-rm  -r output/*