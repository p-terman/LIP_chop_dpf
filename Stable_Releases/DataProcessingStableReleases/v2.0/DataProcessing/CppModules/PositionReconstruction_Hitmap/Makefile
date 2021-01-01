#
# Makefile for many Unix compilers using the
# "standard" command name CC
#
CC=g++
DEBUG=-g
CFLAGS=-Wall -c $(DEBUG) $(BOOST)
LFLAGS=-Wall $(DEBUG)
OBJS=RQFile_IO.o 

ifeq ($(shell uname), Darwin)
CC=clang++
CFLAGS=-std=c++0x $(BOOST) -c $(DEBUG)
LFLAGS=-std=c++0x $(DEBUG)
endif

all: ../bin/PositionReconstruction_Hitmap

clean:
	rm -f *.o ../bin/PositionReconstruction_Hitmap

../bin/PositionReconstruction_Hitmap: ../Utilities/RQFile_IO.cxx PBomb.cpp Pulse.cpp PositionReconstruction_Hitmap.cpp Pulse.h PBomb.h
	g++ -g -Wall $(BOOST) ../Utilities/RQFile_IO.cxx PBomb.cpp Pulse.cpp PositionReconstruction_Hitmap.cpp -o ../bin/PositionReconstruction_Hitmap
