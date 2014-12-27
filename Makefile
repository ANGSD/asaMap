FLAGS=-O3 -D_USE_KNETFILE


CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)

PRG=line


all: $(PRG)

.PHONY: clean

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

-include $(OBJ:.o=.d)

%.o: %.c
	$(CC) -c  $(CFLAGS) $*.c
	$(CC) -MM $(CFLAGS) $*.c >$*.d
%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp >$*.d


line: $(OBJ)
	$(CXX) $(FLAGS)  -o line *.o -lz -lpthread

clean:
	rm  -f *.o *.d line *~
