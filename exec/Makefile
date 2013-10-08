## Makefile

PACKNAME=$(shell basename $(PWD))
PACKDIR=$(shell dirname $(PWD))

#################################################################
ARCH         := $(shell root-config --arch)
CXX           =
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
#################################################################

LHAPDFCFLAGS := # $(shell ../lhapdf/bin/lhapdf-config --cppflags)
LHAPDFLIBS   := # $(shell ../lhapdf/bin/lhapdf-config --ldflags)

PACKLIB = -L../shlib -lmatrix -lintegrator -ldhelas -lcuba -lnr -lttanalysis -lMinuit
PACKOBJ = 

ifndef FF
FF = gfortran
endif

# gcc compiler (i686)
ifeq ($(ARCH),linux)
CXX           = g++
CXXFLAGS      = -g -O3 -Wall -fPIC -m32
FFLAGS        = -O3 -c -x f77-cpp-input
LD            = g++
LDFLAGS       = -g -O2 -shared -m32 
SOFLAGS       = -shared
endif

# gcc compiler (x86_64)
ifeq ($(ARCH),linuxx8664gcc)
CXX           = g++
CXXFLAGS      = -g -O3 -Wall -fPIC -m64
FFLAGS        = -O3 -c -fPIC
LD            = g++
LDFLAGS       = -g -O3 --no_exceptions -shared -m64
SOFLAGS       = -shared
PACKLIB      += ## $(shell locate lib64/libg2c)
endif

ifeq ($(CXX),)
$(error $(ARCH) invalid architecture)
endif

INCFLAGS = $(shell find ../ -maxdepth 1 -type d | sed -e "s:\(.*\):-I\1:g" | grep -v "../shlib" )
CXXFLAGS += $(ROOTCFLAGS) $(LHAPDFCFLAGS) $(INCFLAGS) 

#

SRCS =  $(wildcard src/*.cc) 

INCS = 	$(wildcard $(PACKNAME)/*.hh)

DEPS =  $(patsubst src/%, objects/%, \
		  $(patsubst %.cc, %.d, $(wildcard src/*.cc)) )

# We want dictionaries only for classes that have _linkdef.h files
DICTOBS =  $(patsubst %_linkdef.h, %.o, \
		      $(patsubst dict/%, objects/dict_%, \
			  $(wildcard dict/*_linkdef.h) ) )

OBJS =   $(patsubst src/%, objects/%, \
		    $(patsubst %.cc, %.o, $(SRCS)) )

OBJS += $(DICTOBS)

all: python

bin: 

objects/dict_%.o: $(PACKNAME)/%.hh dict/%_linkdef.h
	@echo "Generating dictionary for $<"
	@$(ROOTSYS)/bin/rootcint -f $(patsubst %.o, %.C, $@) -c -Idict $(INCFLAGS) $(notdir $^) 
	@$(CXX) -c $(CXXFLAGS) -o $@ $(patsubst %.o, %.C, $@)

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
	@rm -f objects/*.d objects/*.o core* objects/dict_* 

python: $(wildcard python/*.py)
	@echo "Copying python sources"
	@cp -rvf python/*.py ../python/.
