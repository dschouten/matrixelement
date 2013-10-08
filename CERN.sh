echo "using: /cvmfs/atlas.cern.ch"

echo "   GCC: /cvmfs/atlas.cern.ch/repo/sw/atlas-gcc/432/$(uname -m)"
source /cvmfs/atlas.cern.ch/repo/sw/atlas-gcc/432/$(uname -m)/setup.sh

echo "   PYTHON: /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/$(uname -m)/python/2.6.5-$(uname -m)-slc5-gcc43"
PYTHONBASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/$(uname -m)/python/2.6.5-$(uname -m)-slc5-gcc43
PYTHONBASE=${PYTHONBASE}/sw/lcg/external/Python/2.6.5/$(uname -m)-slc5-gcc43-opt
export LD_LIBRARY_PATH=${PYTHONBASE}/lib:$LD_LIBRARY_PATH
export PATH=${PYTHONBASE}/bin:$PATH    

echo "   ROOT: /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/$(uname -m)/root/5.28.00g-slc5-gcc4.3"
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/$(uname -m)/root/5.28.00g-slc5-gcc4.3/bin/thisroot.sh

unset HEPPY

source setup.sh
