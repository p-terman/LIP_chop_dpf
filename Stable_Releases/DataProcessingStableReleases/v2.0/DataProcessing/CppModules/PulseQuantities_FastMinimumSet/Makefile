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
OBJS += ../build/Fits.o

all: ../bin/PulseQuantities_FastMinimumSet

../bin/PulseQuantities_FastMinimumSet: $(OBJS) ../build/PulseQuantities_FastMinimumSet.o
	$(CC) $(OSFLAGS) $(RCFLAGS) $(LFLAGS) $(RLFLAGS) $^ $(LIBS) -lMinuit -o $@

../build/PulseQuantities_FastMinimumSet.o: PulseQuantities_FastMinimumSet.cpp Quantities.h Fits.h
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $< -o $@ 

../build/Fits.o: Fits.cxx Fits.h
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $< -o $@ 

include ../build/Makefile.objs

clean:
	rm -f ../bin/PulseQuantities_FastMinimumSet
	rm -f ../build/PulseQuantities_FastMinimumSet.o

distclean:
	rm -f ../bin/PulseQuantities_FastMinimumSet
	rm -f ../build/PulseQuantities_FastMinimumSet.o
	@echo Removing all object files from build directory.
	@rm -f $(OBJS)
