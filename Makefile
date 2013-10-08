NCORES = 1 ## $(shell cat /proc/cpuinfo | wc -l)

all: compile

.PHONY : clean compile dict

compile:
	@mkdir -p lib/
	@mkdir -p shlib/
	@mkdir -p ttanalysis/obj
	@mkdir -p integrator/objects
	@mkdir -p matrix/objects
	@mkdir -p exec/objects

	@cd cuba && make -j$(NCORES) -f Makefile && cd ../
	@cd madgraph/DHELAS && make -j$(NCORES) -f Makefile && cd ../../
	@cd ttanalysis && make -j$(NCORES) -f Makefile && cd ../

	@cd integrator && make -j$(NCORES) -f Makefile && cd ../
	@cd integrator && make -j$(NCORES) -f Makefile dict && cd ../
	@cd matrix && make -j$(NCORES) -f Makefile && cd ../
	@cd matrix && make -j$(NCORES) -f Makefile dict && cd ../
	@cd exec && make -j$(NCORES) -f Makefile && cd ../

clean:
	@cd cuba && make -j$(NCORES) -f Makefile clean && cd ../
	@cd madgraph/DHELAS && make -j$(NCORES) -f Makefile clean && cd ../../
	@cd ttanalysis && make -j$(NCORES) -f Makefile clean && cd ../

	@cd integrator && make -f Makefile clean && cd ../
	@cd matrix && make -f Makefile clean && cd ../
	@cd exec && make -f Makefile clean && cd ../
