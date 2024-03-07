OS := $(shell uname)
ifeq ($(OS),Darwin)
        #CC      = /usr/local/opt/llvm/bin/clang++
        #CFLAGS  = -O3 -mavx -std=c++14 -w -march=native -I/usr/local/opt/llvm/include -fopenmp 
        #LDFLAGS = -L/usr/local/opt/llvm/lib
		CC      = /usr/local/bin/g++-9
        CFLAGS  = -O3 -mavx -std=c++14 -w -march=native -I/usr/local/opt/libomp/include -fopenmp -I/usr/local/include/
        LDFLAGS = -L/usr/local/opt/libomp/lib -L/usr/local/lib/
else
        CC      = g++
        CFLAGS  = -O3 -std=c++14 -w -fopenmp
        DEBUGCFLAGS  = -std=c++14 -w -fopenmp
		VALGRINDCFLAGS  = -O0 -std=c++14 -w -fopenmp
        LDFLAGS = -L/usr/local/lib/ -lgeos_c -I/usr/local/include/
endif

HEADERS= $(shell ls libvbyte-master/*h)

#common
COMPRESSION = libvbyte-master/vbyte.o libvbyte-master/varintdecode.o
SOURCES = containers/relation.cpp $(shell ls topological/*cpp)
OBJECTS = $(SOURCES:.cpp=.o)

all: main

main: $(HEADERS) main.cpp $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(COMPRESSION) main.cpp -o sj $(LDFLAGS) 

debug:  $(OBJECTS)
	$(CC) -g $(DEBUGCFLAGS) $(OBJECTS) $(COMPRESSION) main.cpp -o sj $(LDFLAGS)  

valgrind: $(OBJECTS)
	$(CC) -g $(VALGRINDCFLAGS) $(LDFLAGS)  $(OBJECTS) $(COMPRESSION) main.cpp -o sj -Iinclude

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

.cc.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf containers/*.o
	rm -rf algorithms/*.o
	rm -rf partitioning/*.o
	rm -rf scheduling/*.o
	rm -rf topological/*.o
	rm -rf opengl_rasterizer/*.o
	rm -rf sj
