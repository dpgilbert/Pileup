#
# stuff to make
#
ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else

all: 	build
help:
	@echo "Available Targets:";\
	cat Makefile | perl -ne 'printf("\t%-15s %s\n",$$1,$$2) if(/^(\S+):[^#]+(#.*)$$/)'

ifndef VERBOSE
  QUIET := @
endif

#ROOFITINCLUDE = 
#ifdef CMSSW_VERSION
#	ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
#endif

CC = g++
CMSROOT = ./
INCLUDE = $(shell root-config --cflags) -I$(CMSROOT) -I$(CMSROOT)/CORE
CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(INCLUDE) $(EXTRACFLAGS)
ROOTLIBS = $(shell root-config --ldflags --cflags --libs) -lTMVA #-lEG -lGenVector
COREDIR = CORE

DICTINCLUDE = $(ROOTSYS)/include/Math/QuantFuncMathCore.h $(ROOTSYS)/include/TLorentzVector.h $(ROOTSYS)/include/Math/Vector4D.h

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags --libs) -lEG -lGenVector -lTMVA -lRooFit

DIR = ./

SOURCES = $(wildcard *.cc)
OBJECTS = $(SOURCES:.cc=.o)
LIB = libPileupCorrection.so

CORESOURCES = $(wildcard ../../CORE/*.cc) $(wildcard ../../CORE/Tools/*.cc) 
COREOBJECTS = $(CORESOURCES:.cc=.o) 
CORELIB = libPileupCorrectionCORE.so


DICT = LinkDef_out.o

LIBS = $(LIB) $(CORELIB)

EXE = PileupCorrection.exe

#
# how to make it
#

$(CORELIB): $(DICT) $(COREOBJECTS)
	$(QUIET) echo "Linking $@"; \
	echo "$(LINKER) -shared -o $@ $(COREOBJECTS) $(DICT) $(LINKERFLAGS)"; \
	$(LINKER) -shared -o $@ $(COREOBJECTS) $(DICT) $(LINKERFLAGS)

$(LIB):	$(DICT) $(OBJECTS)
	$(QUIET) echo "Linking $@"; \
	echo "$(LINKER) -shared -o $@ $(OBJECTS) $(DICT) $(LINKERFLAGS)"; \
	$(LINKER) -shared -o $@ $(OBJECTS) $(DICT) $(LINKERFLAGS)

#LinkDef_out.cxx: LinkDef.h
#	$(QUIET) echo "Making CINT dictionaries"; \
	rootcint -f LinkDef_out.cc -c -p $(DICTINCLUDE)  LinkDef.h; \
	cat LinkDef_out.cc > LinkDef_out.cxx; rm LinkDef_out.cc

# the option "-Wl,-rpath,./" adds ./ to the runtime search path for libraries
$(EXE):  $(LIBS)
	$(QUIET) echo "Building $@"; \
	echo "$(CC) -o $@ $(LIBS) $(ROOTLIBS) -Wl,-rpath,./"; \
	$(CC) -o $@ $(LIBS) $(ROOTLIBS) -Wl,-rpath,./ 

%.o: 	%.cc %.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.C %.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cc
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.C
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o:    %.cxx 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

libs:	$(LIBS)

build:  $(EXE)

b: build

clean:
	rm -f \
	LinkDef_out* \
	*.o \
	$(LIBS) \
	$(EXE) \
        ../../CORE/*.o \
        ../../CORE/*.d \
        ../../CORE/*.so \
        ../../CORE/Tools/*.o \
	../../CORE/Tools/*.d \
        ../../CORE/Tools/*.so \

localclean:
	rm -f \
	LinkDef_out* \
	*.o \
	*.so 

endif
