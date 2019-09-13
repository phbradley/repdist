# the compiler
CC = g++

## compiler flags:

# debugging
#CCFLAGS  = -g -std=c++11 -Wall

# production
CCFLAGS  = -O3 -std=c++11 -Wall

# 
INCLUDES = -I ../include/

# recompile if any .hh files changed
HHS = *.hh

all: ../bin/repdist

../bin/repdist:  repdist.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o ../bin/repdist repdist.cc

clean:
	-rm bin/*
