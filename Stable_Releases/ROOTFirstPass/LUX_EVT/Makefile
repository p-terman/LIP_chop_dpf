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
	
EVENTO        = Lux_EVT_Event.$(ObjSuf) Lux_EVT_Channel.$(ObjSuf) Lux_EVT_Pulse.$(ObjSuf) Lux_EVT_Header.$(ObjSuf) XmlParser.$(ObjSuf) eventdict.$(ObjSuf) 
EVENTS        = Lux_EVT_Event.$(SrcSuf) Lux_EVT_Channel.$(SrcSuf) Lux_EVT_Pulse.$(SrcSuf) Lux_EVT_Header.$(SrcSuf) XmlParser.$(SrcSuf) eventdict.$(SrcSuf)
EVENTSO       = libTEvent.$(DllSuf)
EVENT         = Make_EVT_Tree$(ExeSuf)
ifeq ($(PLATFORM),win32)
EVENTLIB      = libTEvent.lib
else
EVENTLIB      = $(shell pwd)/$(EVENTSO)
endif

MAINEVENTO    = Make_EVT_Tree.$(ObjSuf)
MAINEVENTS    = Make_EVT_Tree.$(SrcSuf)

READO       = read_tree.$(ObjSuf)
READS       = read_tree.$(SrcSuf)
READ        = read_tree$(ExeSuf)

WRITEO      = writeRQ1file.$(ObjSuf)
WRITES      = writeRQ1file.$(SrcSuf)
WRITE       = writeRQ1file$(ExeSuf)

OBJS          = $(EVENTO) $(MAINEVENTO) $(READO) $(WRITEO) 

PROGRAMS      = $(EVENT) $(READ) $(WRITE)

#Fix to get this to work in snow leopard
ifeq ($(ARCH),macosx)
SOFLAGS += -m32
endif

XMLFLAGS = -lXMLIO

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
		$(LD) $(LDFLAGS) $(MAINEVENTO) $(EVENTLIB) $(LIBS) $(XMLFLAGS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(READ):        $(READO)
		$(LD) $(LDFLAGS) $^ $(EVENTLIB) $(LIBS) $(XMLFLAGS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(WRITE):        $(WRITEO)
		$(LD) $(LDFLAGS) $^ $(EVENTLIB) $(LIBS) $(XMLFLAGS) $(OutPutOpt)$@
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

Event.$(ObjSuf): Lux_EVT_Event.h
Make_EVT_Tree.$(ObjSuf): Lux_EVT_Event.h

EventDict.$(SrcSuf): Lux_EVT_Event.h TEventLinkDef.h
	@echo "Generating dictionary $@..."
	$(ROOTCINT) -f $@ -c $^
	

.$(SrcSuf).$(ObjSuf):
	$(CXX)  $(CXXFLAGS) -c $<

ifeq ($(GCC_MAJOR),4)
ifeq ($(GCC_MINOR),1)
TBenchDict.o: CXXFLAGS += -Wno-strict-aliasing
endif
endif
