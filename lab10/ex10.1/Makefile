CC = g++
#CFLAGS = -Wall -O3 --std=c++11
CFLAGS = -Wall -Wextra -pedantic
RR = ../../Random

ex10.1.exe : main.o $(RR)/Random.o TSP.o
	$(CC) $(RR)/Random.o TSP.o main.o -o ex10.1.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
TSP.o : TSP.cpp TSP.h
	$(CC) -c TSP.cpp -o TSP.o $(CFLAGS)
$(RR)/Random.o : $(RR)/Random.cpp $(RR)/Random.h
	$(CC) -c $(RR)/Random.cpp -o $(RR)/Random.o $(CFLAGS)
clean :
	rm *.o ex10.1.exe $(RR)/seed.out $(RR)/Random.o
