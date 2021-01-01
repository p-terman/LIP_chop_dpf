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

all: ../bin/PulseCalibration_MakeCvtFile 

../bin/PulseCalibration_MakeCvtFile: $(OBJS) ../build/PulseCalibration_MakeCvtFile.o
	$(CC) $(OSFLAGS) $(LFLAGS) $^ -o $@

../build/PulseCalibration_MakeCvtFile.o: PulseCalibration_MakeCvtFile.cpp
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@

include ../build/Makefile.objs

clean:
	rm -f ../bin/PulseCalibration_MakeCvtFile 
	rm -f ../build/PulseCalibration_MakeCvtFile .o

distclean:
	rm -f ../bin/PulseCalibration_MakeCvtFile 
	rm -f ../build/PulseCalibration_MakeCvtFile .o
	@echo Removing all object files from build directory.
	@rm -f $(OBJS)
