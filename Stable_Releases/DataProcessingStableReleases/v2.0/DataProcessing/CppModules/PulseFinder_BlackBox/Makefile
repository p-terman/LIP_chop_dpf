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

all: ../bin/PulseFinder_BlackBox

../bin/PulseFinder_BlackBox: $(OBJS) ../build/SplitGaussian.o ../build/PulseFinder.o ../build/PulseFinder_BlackBox.o
	$(CC) $(OSFLAGS) $(BOOST) $(RCFLAGS) $(LFLAGS) $(RLFLAGS) $^ $(LIBS) -o $@

../build/PulseFinder_BlackBox.o: PulseFinder_BlackBox.cpp 
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $< -o $@ 

../build/PulseFinder.o: PulseFinder.cxx PulseFinder.h
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $< -o $@ 

../build/SplitGaussian.o: SplitGaussian.cxx SplitGaussian.h  
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $< -o $@ 


include ../build/Makefile.objs

clean:
	rm -f ../bin/PulseFinder_BlackBox
	rm -f ../build/PulseFinder_BlackBox.o
	rm -f ../build/PulseFinder.o
	rm -f ../build/SplitGaussian.o

distclean:
	rm -f ../bin/PulseFinder_BlackBox
	rm -f ../build/PulseFinder_BlackBox.o
	rm -f ../build/PulseFinder.o
	rm -f ../build/SplitGaussian.o
	@echo Removing all object files from build directory.
	@rm -f $(OBJS)
