LIBS = -lfftw3 -lm -lsndfile -lportaudio

all: main

main: main.o
	gcc $(DEBUG) $^ $(LIBS) -o $@

main.o: main.c
	gcc $(DEBUG) -Wall -c $<

clean:
	rm -f *.o main
