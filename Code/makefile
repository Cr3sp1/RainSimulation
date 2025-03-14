CC = clang++
CFLAGS = -Wall -O3 --std=c++11 
DEBUG_CFLAGS = -Wall -O3 --std=c++11 -pedantic -Winline -g
MYLIBS = body.o ray.o RainFunctions.o

main.exe: main.cpp ${MYLIBS}
	$(CC) -o $@ $^ $(CFLAGS) 

%.o : %.cpp
	$(CC) -c $< $(CFLAGS)

run:
	./main.exe

clean :
	rm *.o main.exe
 

debug: CFLAGS = $(DEBUG_CFLAGS)
debug: main.exe