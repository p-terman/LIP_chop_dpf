# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

include Makefile.arch
-include ../MyConfig.mk

#------------------------------------------------------------------------------
	
EVENTO        = LUX_rq1_channel.$(ObjSuf) LUX_rq1_pulse.$(ObjSuf) LUX_rq1_event.$(ObjSuf)  LUX_rq1_header.$(ObjSuf) eventdict.$(ObjSuf) 
EVENTS        = LUX_rq1_channel.$(SrcSuf) LUX_rq1_pulse.$(SrcSuf)  LUX_rq1_event.$(SrcSuf) LUX_rq1_header.$(SrcSuf) eventdict.$(SrcSuf) 
EVENTSO       = libTEvent.$(DllSuf)
EVENT         = read_rq1_built$(ExeSuf)
ifeq ($(PLATFORM),win32)
EVENTLIB      = libTEvent.lib
else
EVENTLIB      = $(shell pwd)/$(EVENTSO)
endif

MAINEVENTO    = read_rq1_built.$(ObjSuf)
MAINEVENTS    = read_rq1_built.$(SrcSuf)

READO       = rq1_plot_ex.$(ObjSuf)
READS       = rq1_plot_ex.$(SrcSuf)
READ        = rq1_plot_ex$(ExeSuf) 

OBJS          = $(EVENTO) $(MAINEVENTO) $(READO)


PROGRAMS      = $(EVENT) $(READ)

ifeq ($(ARCH),macosx)
#CXXFLAGS += -m32
SOFLAGS += -m32
#$(error, "It is working")
endif

#------------------------------------------------------------------------------
		
.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)

$(EVENTSO):     $(EVENTO)
	
ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
ifeq ($(PLATFORM),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def \
		   $(OutPutOpt)$(EVENTLIB)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) \
		   $(OutPutOpt)$@
		$(MT_DLL)
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
endif
endif
endif
		@echo "$@ done"
	
$(EVENT):       $(EVENTSO) $(MAINEVENTO)
		$(LD) $(LDFLAGS) $(MAINEVENTO) $(EVENTLIB) $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(READ):        $(READO)
		$(LD) $(LDFLAGS) $^ $(EVENTLIB) $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(PROGRAMS) $(EVENTSO) $(EVENTLIB) *dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d *.log .def so_locations
		@rm -rf cxx_repository files

.SUFFIXES: .$(SrcSuf)

###

Event.$(ObjSuf): LUX_rq1_event.h
read_rq1_built.$(ObjSuf): LUX_rq1_event.h

EventDict.$(SrcSuf): LUX_rq1_event.h TEventLinkDef.h
	@echo "Generating dictionary $@..."
	$(ROOTCINT) -f $@ -c $^
	

.$(SrcSuf).$(ObjSuf):
	$(CXX)  $(CXXFLAGS) -c $<

#$(EVENTO): %.$(ObjSuf): %.$(SrcSuf)
#	@echo the CXX Flags are: $(CXXFLAGS)
#	$(CXX) $(CXXFLAGS) -c $< -o $@


ifeq ($(GCC_MAJOR),4)
ifeq ($(GCC_MINOR),1)
TBenchDict.o: CXXFLAGS += -Wno-strict-aliasing
endif
endif
