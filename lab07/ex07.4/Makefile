CC = g++
CFLAGS = -Wall -Wextra -pedantic -O3 --std=c++11
RR = ../../Random
MD = ../../MolDyn

.PHONY : all clean
all : NVE.exe NVT.exe
NVE.exe : main_NVE.o NVE.o $(RR)/Random.o $(MD)/MolDyn.o
	$(CC) -o NVE.exe main_NVE.o $(RR)/Random.o NVE.o $(MD)/MolDyn.o
NVT.exe : main_NVT.o NVT.o $(RR)/Random.o $(MD)/MolDyn.o
	$(CC) -o NVT.exe main_NVT.o $(RR)/Random.o NVT.o $(MD)/MolDyn.o
$(MD)/MolDyn.o : $(MD)/MolDyn.cpp $(MD)/MolDyn.h $(RR)/Random.h
	$(CC) -c -o $(MD)/MolDyn.o $(MD)/MolDyn.cpp $(CFLAGS)
NVE.o : NVE.cpp NVE.h $(RR)/Random.h $(MD)/MolDyn.h
	$(CC) -c -o NVE.o NVE.cpp $(CFLAGS)
NVT.o : NVT.cpp NVT.h $(RR)/Random.h $(MD)/MolDyn.h
	$(CC) -c -o NVT.o NVT.cpp $(CFLAGS)
$(RR)/Random.o : $(RR)/Random.cpp $(RR)/Random.h
	$(CC) -c -o $(RR)/Random.o $(RR)/Random.cpp $(CFLAGS)
main_NVT.o : main_NVT.cpp NVT.h
	$(CC) -c -o main_NVT.o main_NVT.cpp $(CFLAGS)
main_NVE.o : main_NVE.cpp NVE.h
	$(CC) -c -o main_NVE.o main_NVE.cpp $(CFLAGS)
clean :
	rm *.o NVT.exe NVE.exe $(RR)/seed.out $(MD)/MolDyn.o
