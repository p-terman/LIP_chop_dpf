##############################################################################
#
# A very simple Makefile
#
# Change log:
#   2012 Dec 04     Initial Submission. (Michael Woods)
##############################################################################

# Avoid the .dSYM produced on a Mac when source code is provided to the linker
# instead of an object file.
CC=g++

ifeq ($(shell uname), Darwin)
CC=clang++
#OSFLAGS+= -gstabs 
endif

all: Livetime

Livetime: Livetime.cc
	$(CC) -g -Wall $(OSFLAGS) -o ../bin/Livetime  Livetime.cc

clean: 
	rm -f ../bin/Livetime

again:
	make clean
	make 


