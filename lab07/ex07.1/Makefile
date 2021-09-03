CC = g++
CFLAGS = -Wall -Wextra -pedantic -O3 --std=c++11
RR = ../../Random
MD = ../../MolDyn

NVT.exe : main.o NVT.o $(RR)/Random.o $(MD)/MolDyn.o
	$(CC) -o NVT.exe main.o $(RR)/Random.o NVT.o $(MD)/MolDyn.o
$(MD)/MolDyn.o : $(MD)/MolDyn.cpp $(MD)/MolDyn.h $(RR)/Random.h
	$(CC) -c -o $(MD)/MolDyn.o $(MD)/MolDyn.cpp $(CFLAGS)
NVT.o : NVT.cpp NVT.h $(RR)/Random.h $(MD)/MolDyn.h
	$(CC) -c -o NVT.o NVT.cpp $(CFLAGS)
$(RR)/Random.o : $(RR)/Random.cpp $(RR)/Random.h
	$(CC) -c -o $(RR)/Random.o $(RR)/Random.cpp $(CFLAGS)
main.o : main.cpp NVT.h
	$(CC) -c -o main.o main.cpp $(CFLAGS)
clean :
	rm *.o NVT.exe $(RR)/seed.out $(MD)/MolDyn.o
