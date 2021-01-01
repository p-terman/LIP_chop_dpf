#
# Makefile for many Unix compilers using the
# "standard" command name CC
#

include ../build/Makefile.config

OBJS  = ../build/LEvtUtilities.o
OBJS += ../build/LEvtChannel.o
OBJS += ../build/LEvtEvent.o
OBJS += ../build/LEvtFile.o
OBJS += ../build/RQFile_IO.o

all: ../bin/PulseFinder_PODoverlap

../bin/PulseFinder_PODoverlap: $(OBJS) ../build/PulseFinder_PODoverlap.o
	$(CC) $(OSFLAGS) $(LFLAGS) $^ -o $@

../build/PulseFinder_PODoverlap.o: PulseFinder_PODoverlap.cpp PODoverlap.h
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@

include ../build/Makefile.objs

clean:
	rm -f ../bin/PulseFinder_PODoverlap
	rm -f ../build/PulseFinder_PODoverlap.o

distclean:
	rm -f ../bin/PulseFinder_PODoverlap
	rm -f ../build/PulseFinder_PODoverlap.o
	@echo Removing all object files from build directory.
	@rm -f $(OBJS)
