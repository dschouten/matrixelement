## Makefile

PACKNAME=$(shell basename $(PWD))

#############################################################
ARCH         := $(shell root-config --arch)
CXX           =
ROOTLDFLAGS  := $(shell root-config --ldflags) 
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
#############################################################

ifndef FF
FF = gfortran
endif

# gcc compiler (i686)
ifeq ($(ARCH),linux)
CXX           = g++
CXXFLAGS      = -g -O3 -Wall -fPIC -m32 -I../external/include/
FFLAGS        = -g -O3 -c -x f77-cpp-input -ffixed-line-length-132 -ffast-math -fstrength-reduce -fexpensive-optimizations
LD            = g++
LDFLAGS       = -g -O3 -shared -m32
SOFLAGS       = -shared
endif

# gcc compiler (x86_64)
ifeq ($(ARCH),linuxx8664gcc)
CXX           = g++
CXXFLAGS      = -g -O3 -Wall -fPIC -m64 -I../external/include/
FFLAGS        = -g -O3 -c -fPIC -ffixed-line-length-132 -ffast-math -fstrength-reduce -fexpensive-optimizations
LD            = g++
LDFLAGS       = -g -O3 -shared
SOFLAGS       = -shared
endif

ifeq ($(CXX),)
$(error $(ARCH) invalid architecture)
endif

SHARED=$(shell dirname $(PWD))/shlib/lib$(PACKNAME).so
DICT=$(shell dirname $(PWD))/shlib/lib$(PACKNAME)_dict.so

INCFLAGS = $(shell find ../ -maxdepth 1 -type d | sed -e "s:\(.*\):-I\1:g" | grep -v "../shlib" )
CXXFLAGS += $(ROOTCFLAGS) $(INCFLAGS) -fnon-call-exceptions

SRCS =  $(wildcard src/*.cc) 
FSRCS = $(wildcard fsrc/*.f)

DICTOBJS =  $(patsubst %_linkdef.h, %.o, \
		      $(patsubst dict/%, objects/dict_%, \
			  $(wildcard dict/*_linkdef.h) ) )

OBJS =   $(patsubst src/%, objects/%, \
		    $(patsubst %.cc, %.o, $(SRCS)) )
FOBJS =  $(patsubst fsrc/%, objects/%, \
		    $(patsubst %.f, %.o, $(FSRCS)) )

EXPLLINKLIBS = -L../shlib -L../external/lib -lcuba -lgsl # -lcubature

all:	shared

shared:	$(SHARED) $(DEPS)

dict: $(DICTOBJS) $(SHARED)
	@$(LD) $(LDFLAGS) $(DICTOBJS) $(EXPLLINKLIBS) $(SHARED) -o $(DICT)

$(SHARED): $(OBJS) $(FOBJS)
	@echo "Creating library $(SHARED)"
	@$(LD) $(LDFLAGS) $(OBJS) $(FOBJS) $(EXPLLINKLIBS) -o $(SHARED) 
	@echo "$(SHARED) successfully compiled!"


objects/dict_%.o: $(PACKNAME)/%.hh dict/%_linkdef.h
	@echo "Generating dictionary for $<"
	@$(ROOTSYS)/bin/rootcint -f $(patsubst %.o, %.C, $@) -c -Idict $(INCFLAGS) $^
	@$(CXX) -c $(CXXFLAGS) -I../$(PACKNAME)/objects/ -o $@ $(patsubst %.o, %.C, $@)

objects/%.o: src/%.cc $(PACKNAME)/%.hh
	@echo "Compiling $<"
	@echo $(CXX) $(CXXFLAGS) -c -o $@ $<
	@$(CXX) $(CXXFLAGS) -c -o $@ $<

objects/%.o: fsrc/%.f
	@echo "Compiling $<"
	@$(FF) $(FFLAGS) -c -o $@ $<

objects/%.d: src/%.cc
	@echo "Generating dependencies for $<"
	@set -e; $(CXX) -M $(CXXFLAGS) $< \
	| sed 's%\($*\)\.o[ :]*%\1.o $@ : %g' > $@; \
	[ -s $@ ] || rm -f $@

objects/%.d: fsrc/%.f
	@echo "Generating dependencies for $<"
	@set -e; $(FF) -M $(FFLAGS) $< \
	| sed 's%\($*\)\.o[ :]*%\1.o $@ : %g' > $@; \
	[ -s $@ ] || rm -f $@

clean:
	@echo "Cleaning everything..."
	@rm -f objects/*.d objects/*.o core* objects/dict_* $(SHARED) $(DICT)

