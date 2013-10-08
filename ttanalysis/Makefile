#
#            M A I N
#           =========

default: lib

info:
	@echo Please type: 'gmake main.exe'
	@echo or:          'gmake lib'
	@echo or:          'gmake clean'

#################################################################
ARCH         := $(shell root-config --arch)

#################################################################
SRC_cc = $(wildcard src/*.cc)
SRC_C += $(wildcard src/*.C)
SRC_cpp += $(wildcard src/*.cpp)


OBJ_cc = $(patsubst src/%.cc, obj/%.o, $(SRC_cc))
OBJ_C = $(patsubst src/%.C, obj/%.o, $(SRC_C))
OBJ_cpp = $(patsubst src/%.cpp, obj/%.o, $(SRC_cpp))
ANAOBJ = $(OBJ_cc) $(OBJ_C) $(OBJ_cpp)

#################################################################

CFLAGS = -g -O3 -fPIC

ifeq ($(ARCH),linuxx8664gcc)
CFLAGS += -m64
endif

ifeq ($(ARCH),linux)
CFLAGS += -m32
endif

INC = -Ittanalysis/

SYSTEM:= $(shell uname)

##################################################################
ifeq ($(SYSTEM),Linux)
  CPPCOMP = gcc -Wno-deprecated -pthread
  LINKER = gcc -lstdc++
endif
##################################################################

obj/%.o :: src/%.cc ttanalysis/*.h*
	$(CPPCOMP) -c $(CFLAGS) $(INC) $< -o $@

obj/%.o :: src/%.C ttanalysis/*.h*
	$(CPPCOMP) -c $(CFLAGS) $(INC) $< -o $@

obj/%.o :: src/%.cpp ttanalysis/*.h*
	$(CPPCOMP) -c $(CFLAGS) $(INC) $< -o $@

##############################################################

main.exe:$(ANAOBJ)
	$(LINKER) $(ANAOBJ) \
	-o main.exe

lib: obj/tTSolver.o
	$(LINKER) $(CFLAGS) -shared obj/tTSolver.o -o ../shlib/libttanalysis.so

clean:
	rm -rf obj/*.o
	rm -f *.x
	rm -rf ../shlib/libttanalysis.so




