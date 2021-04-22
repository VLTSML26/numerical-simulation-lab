CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o Random.o
	$(CC) Random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
Random.o : Random.cpp Random.h
	$(CC) -c Random.cpp -o Random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
