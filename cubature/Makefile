PACKNAME=$(shell basename $(PWD))

ARCH         := $(shell root-config --arch)
CXX           =

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

LIBS = 

ifndef FF
FF = gfortran
endif

# gcc optimization: see http://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html
# g77 optimization: see http://gcc.gnu.org/onlinedocs/g77/Optimize-Options.html

# gcc compiler (i686)
ifeq ($(ARCH),linux)
CXX           = g++
CXXFLAGS      = -g -O3 -Wall -fPIC -m32
FFLAGS        = -O3 -c -x f77-cpp-input -ffixed-line-length-132 -fno-second-underscore -ext-names
LD            = g++
LDFLAGS       = -g -O3 -shared $(LIBS) 
SOFLAGS       = -shared
endif

# gcc compiler (x86_64)
ifeq ($(ARCH),linuxx8664gcc)
CXX           = g++
CXXFLAGS      = -g -O3 -Wall -fPIC -m64
FFLAGS        = -O3 -c -fPIC -ffixed-line-length-132 -ffast-math -fstrength-reduce -fexpensive-optimizations -fno-second-underscore -ext-names
LD            = g++
LDFLAGS       = -g -O3 --no_exceptions -shared $(LIBS) 
SOFLAGS       = -shared
endif

ifeq ($(CXX),)
$(error $(ARCH) invalid architecture)
endif

SHARED = $(shell dirname $(PWD))/shlib/lib$(PACKNAME).so

INCFLAGS = $(shell find ../ -maxdepth 1 -type d | sed -e "s:\(.*\):-I\1:g" | grep -v "../shlib" )
CXXFLAGS += $(ROOTCFLAGS) $(INCFLAGS) 

SRCS =  $(wildcard *.c) 

INCS = 	$(wildcard $(PACKNAME)/*.h)

OBJS =   $(patsubst %, objects/%, \
		    $(patsubst %.c, %.o, $(SRCS)) )

all:	shared

shared:	$(SHARED) $(DEPS)

$(SHARED): $(OBJS)
	@echo "Creating library $(SHARED)"
	@echo "@$(LD) $(LDFLAGS) $(OBJS) $(EXPLLINKLIBS) -o $(SHARED)"
	@$(LD) $(LDFLAGS) $(OBJS) $(EXPLLINKLIBS) -o $(SHARED)
	@echo "$(SHARED) successfully compiled!"

objects/%.o: %.c %.h 
	@echo "Compiling $<"
	@echo $(CXX) $(CXXFLAGS) -c -o $@ $<
	@$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	@echo "Cleaning everything..."
	@rm -f objects/*.o core* $(SHARED)


