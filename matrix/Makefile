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

# Remove LHAPDF libs and just use CTEQ standalone

LHAPDFCFLAGS := # $(shell ../lhapdf/bin/lhapdf-config --cppflags)
LHAPDFLIBS   := # $(shell ../lhapdf/bin/lhapdf-config --ldflags)

LIBS = $(ROOTLIBS) $(LHAPDFLIBS) -L../shlib -lttanalysis

ifndef FF
  FF = gfortran
endif

F90 = gfortran # f95

# gcc optimization: see http://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html
# g77 optimization: see http://gcc.gnu.org/onlinedocs/g77/Optimize-Options.html

# gcc compiler (i686)
ifeq ($(ARCH),linux)
CXX           = g++
CXXFLAGS      = -g -O3 -Wall -fPIC -m32 -Wpacked -malign-double -mpreferred-stack-boundary=8 -I../external/include/
FFLAGS        = -O3 -m32 -c -fPIC -x f77-cpp-input -ffixed-line-length-132 -fno-second-underscore -ext-names -malign-double
F90FLAGS      = -O3 -std=gnu -m32 -fPIC -ffree-line-length-none -I objects/
LD            = g++
LDFLAGS       = -g -O3 -m32 -shared -fPIC $(LIBS) $(ROOTLDFLAGS)
SOFLAGS       = -shared
LIBS         += ## -L/usr/lib -lg2c 
endif

# gcc compiler (x86_64)
ifeq ($(ARCH),linuxx8664gcc)
CXX           = g++
CXXFLAGS      = -g -O3 -Wall -fPIC -m64 -Wpacked -I../external/include/
FFLAGS        = -O3 -m64 -c -fPIC -ffixed-line-length-132 -ffast-math -fstrength-reduce -fexpensive-optimizations -fno-second-underscore -ext-names 
F90FLAGS      = -O3 -std=gnu -m64 -fPIC -ffree-line-length-none -I objects/
LD            = g++
LDFLAGS       = -g -O3 -shared -fPIC $(LIBS) $(ROOTLDFLAGS)
SOFLAGS       = -shared
LIBS         += ## /usr/lib64/$(shell /bin/ls /usr/lib64 | grep libg2c | head -1 )
endif

ifeq ($(CXX),)
$(error $(ARCH) invalid architecture)
endif

SHARED=$(shell dirname $(PWD))/shlib/lib$(PACKNAME).so
DICT=$(shell dirname $(PWD))/shlib/lib$(PACKNAME)_dict.so

INCFLAGS = $(shell find ../ -maxdepth 1 -type d | sed -e "s:\(.*\):-I\1:g" | grep -v "../shlib" )
CXXFLAGS += $(ROOTCFLAGS) $(LHAPDFCFLAGS) $(INCFLAGS) 

CXXSRCS =  $(wildcard src/*.cc) 
FSRCS = $(wildcard fsrc/*.f)
F90SRCS = $(wildcard fsrc/*.f90)

INCS = 	$(wildcard $(PACKNAME)/*.hh)

CXXDEPS =  $(patsubst src/%, objects/%, \
		$(patsubst %.cc, %.d, $(wildcard src/*.cc)) )

FDEPS = $(patsubst fsrc/%, objects/%, \
		$(patsubst %.f, %.d, $(wildcard fsrc/*.f)) )
F90DEPS = $(patsubst fsrc/%, objects/%, \
		$(patsubst %.f90, %.d, $(wildcard fsrc/*.f90)) )

DICTOBJS =  $(patsubst %_linkdef.h, %.o, \
		$(patsubst dict/%, objects/dict_%, \
		$(wildcard dict/*_linkdef.h) ) )

CXXOBJS =   $(patsubst src/%, objects/%, \
		$(patsubst %.cc, %.o, $(CXXSRCS)) )
FOBJS =  $(patsubst fsrc/%, objects/%, \
		$(patsubst %.f, %.o, $(FSRCS)) )
F90OBJS =  $(patsubst fsrc/%, objects/%, \
		$(patsubst %.f90, %.o, $(F90SRCS)) )

all:	shared

shared:	$(SHARED) $(CXXDEPS)

dict: $(DICTOBJS) $(SHARED)
	@$(LD) $(LDFLAGS) $(DICTOBJS) $(EXPLLINKLIBS) $(SHARED) -o $(DICT)

$(SHARED): $(CXXOBJS) $(FOBJS) $(F90OBJS)
	@echo "Creating library $(SHARED)"
	@echo $(LD) $(LDFLAGS) $(CXXOBJS) $(FOBJS) $(F90OBJS) $(EXPLLINKLIBS) -o $(SHARED)
	@$(LD) $(LDFLAGS) $(CXXOBJS) $(FOBJS) $(F90OBJS) $(EXPLLINKLIBS) -o $(SHARED)
	@echo "$(SHARED) successfully compiled!"

objects/dict_%.o: $(PACKNAME)/%.hh dict/%_linkdef.h
	@echo "Generating dictionary for $<"
	@$(ROOTSYS)/bin/rootcint -f $(patsubst %.o, %.C, $@) -c -Idict $(INCFLAGS) $^
	@$(CXX) -c $(CXXFLAGS) -o $@ $(patsubst %.o, %.C, $@)

objects/%.d: src/%.cc
	@echo "Generating dependencies for $<"
	@$(CXX) $(CXXFLAGS) -MF $@ -MM -MT $@ -MT $(basename $@).o $<

objects/%.o: src/%.cc $(PACKNAME)/%.hh 
	@echo "Compiling $<"
	@$(CXX) $(CXXFLAGS) -c -o $@ $<

objects/%.o: fsrc/%.f
	@echo "Compiling $<"
	@$(FF) $(FFLAGS) -c -o $@ $<

objects/%.o: fsrc/%.f90
	@echo "Compiling $<"
	@$(F90) $(F90FLAGS) -c -o $@ $<

clean:
	@echo "Cleaning everything..."
	@rm -f objects/*.d objects/*.o core* objects/dict_* $(SHARED) $(DICT)

##
## kludge to get the module dependencies right in Fortran-90
##

objects/mod_Graviton.o: objects/mod_Parameters.o
objects/mod_Graviton.o: objects/mod_Misc.o
objects/mod_Higgs.o: objects/mod_Parameters.o
objects/mod_Higgs.o: objects/mod_Misc.o

include $(CXXDEPS)
