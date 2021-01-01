
include Makefile.arch

# Flags ------------------------------------------------------------------------

SVNV      := `svnversion -n`
CXXFLAGS  += -g 
#LIBS      += -lNew		# Don't know why but this broke comp'n on Mac.
#GLIBS     += -lMinuit 


# ROOT->RQ1 Converter ----------------------------------------------------------

RQ12ROOT_MAIN_E := ../bin/rqroot_to_rq 
RQ12ROOT_MAIN_S := rqroot_to_rq.cxx 

# RQ1->ROOT Converter ----------------------------------------------------------

ROOT2RQ1_MAIN_E := ../bin/rq_to_rqroot
ROOT2RQ1_MAIN_S := rq_to_rqroot.cxx 

# RQ1->ASCII Converter ---------------------------------------------------------

RQ12ASCII_MAIN_E := ../bin/rq_to_ascii
RQ12ASCII_MAIN_S := rq_to_ascii.cxx 

# Overall ----------------------------------------------------------------------

PROGRAMS  := $(RQ12ROOT_MAIN_E) $(ROOT2RQ1_MAIN_E) $(RQ12ASCII_MAIN_E)

#------------------------------------------------------------------------------

ifeq ($(shell uname), Darwin)
#CXXFLAGS     +=  -gstabs
endif

all: $(PROGRAMS) 

$(RQ12ROOT_MAIN_E): $(RQ12ROOT_MAIN_S) 
	@bash -c 'echo -n "Compiling converter $@... "'
	@$(CXX) $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ 
	@echo " done"

$(ROOT2RQ1_MAIN_E): $(ROOT2RQ1_MAIN_S) 
	@bash -c 'echo -n "Compiling converter $@... "'
	@$(CXX) $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ 
	@echo " done"

$(RQ12ASCII_MAIN_E): $(RQ12ASCII_MAIN_S) 
	@bash -c 'echo -n "Compiling converter $@... "'
	@$(CXX) $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ 
	@echo " done"

$(ROOT2RQ1_MAIN_S): 

$(RQ12ROOT_MAIN_S): 

$(RQ12ASCII_MAIN_S):

clean:
		@rm -f $(PROGRAMS)

distclean:  clean
