PACKNAME=$(shell basename $(PWD))
PACKDIR=$(shell dirname $(PWD))

ARCH         := $(shell root-config --arch)
CXX           =

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

PACKLIB = 
PACKOBJ = 

# gcc compiler (i686)
ifeq ($(ARCH),linux)
CXX           = g++
FXX           = g77
CXXFLAGS      = -g -O3 -Wall -fPIC -m32
FXXFLAGS      = -O3 -c -x f77-cpp-input
LD            = g++
LDFLAGS       = -g -O2 -shared -m32 
SOFLAGS       = -shared
endif

# gcc compiler (x86_64)
ifeq ($(ARCH),linuxx8664gcc)
CXX           = g++
FXX           = g77
CXXFLAGS      = -g -O3 -Wall -fPIC -m64
FXXFLAGS      = -O3 -c -fPIC
LD            = g++
LDFLAGS       = -g -O3 --no_exceptions -shared -m64
SOFLAGS       = -shared
endif

ifeq ($(CXX),)
$(error $(ARCH) invalid architecture)
endif

INCFLAGS = $(shell find ../ -maxdepth 1 -type d | sed -e "s:\(.*\):-I\1:g" | grep -v "../shlib" )
CXXFLAGS += $(ROOTCFLAGS) $(INCFLAGS) 

#

SRCS =  $(wildcard src/*.cpp) 

INCS = 	$(wildcard $(PACKNAME)/*.h)

DEPS =  $(patsubst src/%, objects/%, \
		  $(patsubst %.cpp, %.d, $(wildcard src/*.cpp)) )

# We want dictionaries only for classes that have _linkdef.h files
DICTOBS =  $(patsubst %_linkdef.h, %.o, \
		      $(patsubst dict/%, objects/dict_%, \
			  $(wildcard dict/*_linkdef.h) ) )

OBJS =   $(patsubst src/%, objects/%, \
		    $(patsubst %.cpp, %.o, $(SRCS)) )

OBJS += $(DICTOBS)

all: libmtt.so

libmtt.so: $(OBJS) $(DICTOBJS)
	@echo "Creating library libmtt.so"
	@echo "@$(LD) $(LDFLAGS) $(OBJS) $(EXPLLINKLIBS) -o ../shlib/libmtt.so"
	@$(LD) $(LDFLAGS) $(OBJS) $(EXPLLINKLIBS) -o ../shlib/libmtt.so
	@echo "libmtt.so successfully compiled!"

objects/dict_%.o: $(PACKNAME)/%.h dict/%_linkdef.h
	@echo "Generating dictionary for $<"
	@$(ROOTSYS)/bin/rootcint -f $(patsubst %.o, %.C, $@) -c -Idict $(INCFLAGS) $(notdir $^) 
	@$(CXX) -c $(CXXFLAGS) -o $@ $(patsubst %.o, %.C, $@)

objects/%.o: src/%.cpp $(PACKNAME)/%.h
	@echo "Compiling $<"
	@echo $(CXX) $(CXXFLAGS) -c -o $@ $<
	@$(CXX) $(CXXFLAGS) -c -o $@ $<

objects/%.o: fsrc/%.f
	@echo "Compiling $<"
	@$(FXX) $(FXXFLAGS) -c -o $@ $<

objects/%.d: src/%.cpp
	@echo "Generating dependencies for $<"
	@set -e; $(CXX) -M $(CXXFLAGS) $< \
	| sed 's%\($*\)\.o[ :]*%\1.o $@ : %g' > $@; \
	[ -s $@ ] || rm -f $@

objects/%.d: fsrc/%.f
	@echo "Generating dependencies for $<"
	@set -e; $(FXX) -M $(FXXFLAGS) $< \
	| sed 's%\($*\)\.o[ :]*%\1.o $@ : %g' > $@; \
	[ -s $@ ] || rm -f $@

clean:
	@echo "Cleaning everything..."
	@rm -f objects/*.d objects/*.o core* objects/dict_* $(SHARED)

python: $(wildcard python/*.py)
	@echo "Copying python sources"
	@cp -rf $< ../python/.
