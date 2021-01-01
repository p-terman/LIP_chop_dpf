#
# Makefile for many Unix compilers using the
# "standard" command name CC
#
CC=g++
DEBUG=-g
CFLAGS=-Wall -std=c++0x -c $(DEBUG) $(BOOST)
LFLAGS=-Wall -std=c++0x $(DEBUG) 
OBJS=LEvtUtilities.o LEvtChannel.o LEvtEvent.o LEvtFile.o LCvtChannel.o LCvtEvent.o LCvtFile.o RQFile_IO.o

ifeq ($(shell uname), Darwin)
CC=clang++
#OSFLAGS+= -gstabs 
endif


PulseQuantities_PseudoPromptFraction: $(OBJS) PulseQuantities_PseudoPromptFraction.o
	$(CC) $(OSFLAGS) $(LFLAGS) $^ -o ../bin/$@

PulseQuantities_PseudoPromptFraction.o: PulseQuantities_PseudoPromptFraction.cpp PseudoPromptFraction.h
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@ 

clean:
	\rm -f *.o ../bin/PulseQuantities_PseudoPromptFraction

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

LCvtChannel.o: ../Utilities/LCvtChannel.cpp ../Utilities/LCvtChannel.h  ../Utilities/LEvtChannel.h ../Utilities/LEvtChannel.cpp
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@ 

LCvtEvent.o: ../Utilities/LCvtEvent.cpp ../Utilities/LCvtEvent.h ../Utilities/LEvtEvent.h ../Utilities/LEvtEvent.cpp
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@ 

LCvtFile.o: ../Utilities/LCvtFile.cpp ../Utilities/LCvtFile.h ../Utilities/LCvtEvent.h ../Utilities/LEvtFile.h ../Utilities/LEvtFile.cpp ../Utilities/LEvtEvent.h
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@ 

