##############################################################################
#
# A very simple Makefile
#
# Change log:
#   2012 Dec 04     Initial Submission. (Michael Woods)
##############################################################################

# Avoid the .dSYM produced on a Mac when source code is provided to the linker
# instead of an object file.
ifeq ($(shell uname), Darwin)
OSFLAGS     +=  -gstabs
endif

all: Livetime

Livetime: Livetime.cc
	g++ -g -Wall $(OSFLAGS) -o ../bin/Livetime  Livetime.cc

clean: 
	rm -f ../bin/Livetime

again:
	make clean
	make 


