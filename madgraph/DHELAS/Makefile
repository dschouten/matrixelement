#################################################################
ARCH         := $(shell root-config --arch)
#################################################################

ifndef FF
  FF = gfortran
endif

ifeq ($(origin FFLAGS),undefined)
  FFLAGS= -O -w
endif

ifeq ($(ARCH),linuxx8664gcc)
  FFLAGS += -m64
  LDFLAGS += -m64
endif

ifeq ($(ARCH),linux)
  FFLAGS += -m32
  LDFLAGS += -m32
endif

FFLAGS += -ffixed-line-length-132 -fPIC -fno-second-underscore

UNAME := $(shell uname)

ifdef static
  LIBEXT = a
  CREATELIB = $(AR) cru
else
  LIBEXT = so
  FFLAGS+= -fPIC
  LDFLAGS += -shared
  CREATELIB = $(FF) $(FFLAGS) $(LDFLAGS) -o
endif

LIBRARY	= libdhelas.$(LIBEXT)
LIBDIR  = ../../shlib/

BASIC_OBJS = AlohaFunctions.o

include AlohaMakefile.inc

all: $(LIBDIR)$(LIBRARY)

%.o: %.f
	@echo "Compiling $<"
	@$(FF) $(FFLAGS) -c -o $@ $<

$(LIBDIR)$(LIBRARY): $(BASIC_OBJS) $(ALOHARoutine)
	$(CREATELIB) $@ $^

clean:
	$(RM) *.o $(LIBDIR)$(LIBRARY)
