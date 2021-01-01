#
# Makefile for many Unix compilers using the
# "standard" command name CC
#

include ../build/Makefile.config

OBJS  = ../build/LEvtUtilities.o
OBJS += ../build/LEvtChannel.o
OBJS += ../build/LEvtEvent.o
OBJS += ../build/LEvtFile.o
OBJS += ../build/LCvtChannel.o
OBJS += ../build/LCvtEvent.o
OBJS += ../build/LCvtFile.o
OBJS += ../build/RQFile_IO.o

all: ../bin/PulseTiming_CompensationForTRC

../bin/PulseTiming_CompensationForTRC: $(OBJS) ../build/PulseTiming_CompensationForTRC.o
	$(CC) $(OSFLAGS) $(LFLAGS) $^ -o $@

../build/PulseTiming_CompensationForTRC.o: PulseTiming_CompensationForTRC.cpp TRCAlgorithms.h
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@

include ../build/Makefile.objs

clean:
	rm -f ../bin/PulseTiming_CompensationForTRC
	rm -f ../build/PulseTiming_CompensationForTRC.o

distclean:
	rm -f ../bin/PulseTiming_CompensationForTRC
	rm -f ../build/PulseTiming_CompensationForTRC.o
	@echo Removing all object files from build directory.
	@rm -f $(OBJS)
