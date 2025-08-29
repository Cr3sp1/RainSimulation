CC = clang++
CFLAGS = -Wall -O3 --std=c++17 
DEBUG_CFLAGS = -Wall -O3 --std=c++17 -pedantic -Winline -g
MYLIBS = body.o ray.o RainFunctions.o tinyxml2.o

main.exe: main.cpp ${MYLIBS}
	$(CC) -o $@ $^ $(CFLAGS) 

%.o: %.cpp
	$(CC) -c $< $(CFLAGS)

run:
	./main.exe

clean:
	rm *.o main.exe

format :
	clang-format -i *.cpp *.h
 

debug: CFLAGS = $(DEBUG_CFLAGS)
debug: main.exe