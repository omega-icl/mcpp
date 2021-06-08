CC = g++
CFLAGS = -I..

.cpp.o:
	$(CC) $(CFLAGS) -c $<

all: remez

remez : remez.o
	$(CC) -o $@ $@.o

clean:
	-rm *.o remez
	
