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

all: ../bin/PulseQuantities_WaterPmtRQs

../bin/PulseQuantities_WaterPmtRQs: $(OBJS) ../build/PulseQuantities_WaterPmtRQs.o
	$(CC) $(OSFLAGS) $(LFLAGS) $^ -o $@

../build/PulseQuantities_WaterPmtRQs.o: PulseQuantities_WaterPmtRQs.cpp Timing.h Quantities.h
	$(CC) $(OSFLAGS) $(CFLAGS) $< -o $@

include ../build/Makefile.objs

clean:
	rm -f ../bin/PulseQuantities_WaterPmtRQs
	rm -f ../build/PulseQuantities_WaterPmtRQs.o

distclean:
	rm -f ../bin/PulseQuantities_WaterPmtRQs
	rm -f ../build/PulseQuantities_WaterPmtRQs.o
	@echo Removing all object files from build directory.
	@rm -f $(OBJS)
