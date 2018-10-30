
all:lsh cube

lsh: lsh.o functions.o hashtable.o
	gcc  -o  lsh  lsh.o functions.o hashtable.o -lm

cube: cube.o functions.o hashtable.o
	gcc -o cube cube.o functions.o hashtable.o -lm

lsh.o: lsh.c 
	gcc -c lsh.c 

cube.o: cube.c 
	gcc -c cube.c 

functions.o: functions.c 
	gcc -c functions.c 

hashtable.o: hashtable.c 
	gcc -c hashtable.c 

clean:
	-rm lsh.o functions.o hashtable.o cube.o lsh cube
	-rm  -r output/*