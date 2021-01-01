#
# Makefile for many Unix compilers using the
# "standard" command name CC
#
CC=g++
DEBUG=-g
CFLAGS=-Wall -std=c++0x -c $(DEBUG) $(BOOST)
LFLAGS=-Wall -std=c++0x $(DEBUG)
OBJS=RQFile_IO.o

ifeq ($(shell uname), Darwin)
CC=clang++
#OSFLAGS+= -gstabs 
endif

PROG=PositionReconstruction_CorrCentroid

all: $(OBJS) $(PROG)

$(PROG): $(OBJS) $(PROG).o
	$(CC) $(LFLAGS) $^  -o ../bin/$@ 

$(PROG).o:  $(PROG).cpp CorCentroid.h 
	$(CC) $(CFLAGS) $< -c -o $@ 

clean:
	\rm -f *.o ../bin/$(PROG)


# do not change these
RQFile_IO.o: ../Utilities/RQFile_IO.cxx ../Utilities/RQFile_IO.hh 
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@ 
