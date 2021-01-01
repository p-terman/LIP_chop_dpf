#
# Makefile for many Unix compilers using the
# "standard" command name CC
#
CC=g++
DEBUG=-g
CFLAGS=-Wall -std=c++0x -c $(DEBUG) $(BOOST)
LFLAGS=-Wall -std=c++0x $(DEBUG)
OBJS=LEvtUtilities.o LEvtChannel.o LEvtEvent.o LEvtFile.o RQFile_IO.o 

ifeq ($(shell uname), Darwin)
CC=clang++
#OSFLAGS+= -gstabs 
endif

EXEC=PulseQualityCheck_CppTemplate
all: $(OBJS) $(EXEC)

$(EXEC): $(OBJS) $(EXEC).o
	$(CC) $(OSFLAGS) $(LFLAGS) $^ -o ../bin/$@

$(EXEC).o: $(EXEC).cpp
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@ 

clean:
	\rm -f *.o ../bin/PulseQualityCheck_CppTemplate

# do not change these
LEvtUtilities.o: ../Utilities/LEvtUtilities.cpp ../Utilities/LEvtUtilities.h
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@ 

LEvtChannel.o: ../Utilities/LEvtChannel.cpp ../Utilities/LEvtChannel.h 
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@

LEvtEvent.o: ../Utilities/LEvtEvent.cpp ../Utilities/LEvtEvent.h 
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@ 

LEvtFile.o: ../Utilities/LEvtFile.cpp ../Utilities/LEvtFile.h ../Utilities/LEvtEvent.h
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@

RQFile_IO.o: ../Utilities/RQFile_IO.cxx ../Utilities/RQFile_IO.hh 
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@ 

