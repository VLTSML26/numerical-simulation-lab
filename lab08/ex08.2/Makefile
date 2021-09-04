CC = g++
CFLAGS = -Wall -O3 --std=c++11
RR = ../../Random

main.exe : main.o $(RR)/Random.o
	$(CC) $(RR)/Random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
$(RR)/Random.o : $(RR)/Random.cpp $(RR)/Random.h
	$(CC) -c $(RR)/Random.cpp -o $(RR)/Random.o $(CFLAGS)
clean :
	rm *.o main.exe $(RR)/seed.out $(RR)/Random.o
