CC = g++
CFLAGS = -Wall -Wextra -O3 --std=c++11 #-pedantic
RR = ../../Random
MD = ../../MolDyn

ex04.1.exe : main.o NVE.o $(MD)/MolDyn.o $(RR)/Random.o
	$(CC) -o ex04.1.exe main.o NVE.o $(MD)/MolDyn.o $(RR)/Random.o
$(MD)/MolDyn.o : $(MD)/MolDyn.cpp $(MD)/MolDyn.h $(RR)/Random.h
	$(CC) -c -o $(MD)/MolDyn.o $(MD)/MolDyn.cpp $(CFLAGS)
NVE.o : NVE.cpp NVE.h $(MD)/MolDyn.h $(RR)/Random.h
	$(CC) -c -o NVE.o NVE.cpp $(CFLAGS)
$(RR)/Random.o : $(RR)/Random.cpp $(RR)/Random.h
	$(CC) -c -o $(RR)/Random.o $(RR)/Random.cpp $(CFLAGS)
main.o : main.cpp NVE.h
	$(CC) -c -o main.o main.cpp $(CFLAGS)
clean :
	rm *.o *.exe $(RR)/seed.out $(RR)/Random.o $(MD)/MolDyn.o
