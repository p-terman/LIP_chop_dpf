#
# Makefile for many Unix compilers using the
# "standard" command name CC
#

include ../build/Makefile.config

OBJS += ../build/RQFile_IO.o ../build/PHET.o

all: ../bin/PulseClassifier_PHET

../bin/PulseClassifier_PHET: $(OBJS) ../build/PulseClassifier_PHET.o
	$(CC) $(OSFLAGS) $(BOOST) $(RCFLAGS) $(LFLAGS) $(RLFLAGS) $^ $(LIBS) -o $@

../build/PulseClassifier_PHET.o: PulseClassifier_PHET.cpp PHET.h
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $< -o $@ 

../build/PHET.o: PHET.cxx PHET.h
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $< -o $@ 



include ../build/Makefile.objs

clean:
	rm -f ../bin/PulseClassifier_PHET
	rm -f ../build/PulseClassifier_PHET.o
	rm -f ../build/PHET.o

distclean:
	rm -f ../bin/PulseClassifier_PHET
	rm -f ../build/PulseClassifier_PHET.o
	rm -f ../build/PHET.o
	@echo Removing all object files from build directory.
	@rm -f $(OBJS)
