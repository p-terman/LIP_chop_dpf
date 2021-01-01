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

all: ../bin/PulseCalibration_BaselineSanity 

../bin/PulseCalibration_BaselineSanity: $(OBJS) ../build/PulseCalibration_BaselineSanity.o ../build/BaselineSanity.o
	$(CC) $(OSFLAGS) $(LFLAGS) $^ -o $@

../build/PulseCalibration_BaselineSanity.o: PulseCalibration_BaselineSanity.cpp
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@

../build/BaselineSanity.o: BaselineSanity.cxx BaselineSanity.h
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@

include ../build/Makefile.objs

clean:
	rm -f ../bin/PulseCalibration_BaselineSanity
	rm -f ../build/PulseCalibration_BaselineSanity.o
	rm -f ../build/BaselineSanity.o

distclean:
	rm -f ../bin/PulseCalibration_BaselineSanity
	rm -f ../build/PulseCalibration_BaselineSanity.o
	rm -f ../build/BaselineSanity.o
	@echo Removing all object files from build directory.
	@rm -f $(OBJS)
