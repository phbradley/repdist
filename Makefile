# the compiler
CC = g++

## compiler flags:

# debugging
CCFLAGS  = -g -std=c++11 -Wall

# production
#CCFLAGS  = -O3 -std=c++11 -Wall

# EDIT THIS to point to the location where you installed boost
# this should be the main directory containing boost/ doc/ libs/ etc as subdirectories
BOOSTDIR = /home/pbradley/include/boost_1_67_0

# 
INCLUDES = -I ./include/  -I $(BOOSTDIR)

# recompile if any .hh files changed
HHS = src/*.hh

all: bin/repdist

bin/repdist:  src/repdist.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o bin/repdist src/repdist.cc

clean:
	-rm bin/*
