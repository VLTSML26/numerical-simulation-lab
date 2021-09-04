CC = g++
CFLAGS = -Wall -Wextra -pedantic -O3 --std=c++11
RR = ../../Random

ex01.1.exe : main.o $(RR)/Random.o
	$(CC) $(RR)/Random.o main.o -o ex01.1.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
$(RR)/Random.o : $(RR)/Random.cpp $(RR)/Random.h
	$(CC) -c $(RR)/Random.cpp -o $(RR)/Random.o $(CFLAGS)
clean :
	rm *.o ex01.1.exe $(RR)/seed.out
