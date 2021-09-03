CC = g++
CFLAGS = -Wall -Wextra -pedantic -O3 --std=c++11
RR = ../../Random

Ising.exe : main.o Ising.o $(RR)/Random.o
	$(CC) -o Ising.exe main.o $(RR)/Random.o Ising.o 
Ising.o : Ising.cpp Ising.h $(RR)/Random.h
	$(CC) -c -o Ising.o Ising.cpp $(CFLAGS)
$(RR)/Random.o : $(RR)/Random.cpp $(RR)/Random.h
	$(CC) -c $(RR)/Random.cpp -o $(RR)/Random.o $(CFLAGS)
main.o : main.cpp Ising.h
	$(CC) -c -o main.o main.cpp $(CFLAGS)
clean :
	rm *.o Ising.exe $(RR)/seed.out
