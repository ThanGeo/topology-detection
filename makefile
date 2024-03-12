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
        CFLAGS  = -std=c++14 -w -fopenmp
		DEBUG_FLAGS =
		VALGRINDCFLAGS  = -O0 -std=c++14 -w -fopenmp
        LDFLAGS = 
endif

HEADERS= $(shell ls libvbyte-master/*h)

#common
COMPRESSION = libvbyte-master/vbyte.o libvbyte-master/varintdecode.o
SOURCES = containers/relation.cpp $(shell ls APRIL/*cpp)
OBJECTS = $(SOURCES:.cpp=.o)

all: CFLAGS += -O3
all: main

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) $(COMPRESSION) main.cpp -o sj -Iinclude

debug:	DEBUG_FLAGS += -g
debug: $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) $(COMPRESSION) main.cpp -o sj -Iinclude

valgrind: $(OBJECTS)
	$(CC) -g $(VALGRINDCFLAGS) $(LDFLAGS) $(OBJECTS) $(COMPRESSION) main.cpp -o sj -Iinclude

.cpp.o:
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) -c $< -o $@

.cc.o:
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) -c $< -o $@

clean:
	rm -rf containers/*.o
	rm -rf algorithms/*.o
	rm -rf partitioning/*.o
	rm -rf scheduling/*.o
	rm -rf APRIL/*.o
	rm -rf opengl_rasterizer/*.o
	rm -rf sj
	rm -rf *.o