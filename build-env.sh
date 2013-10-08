#Setup GCC
export XS_GCC_PATH=/afs/cern.ch/sw/lcg/external/gcc/4.6.2/x86_64-slc5-gcc46-opt/
source $XS_GCC_PATH/setup.sh
export LD_LIBRARY_PATH=$XS_GCC_PATH/lib64:$LD_LIBRARY_PATH

#Setup ROOT
export XS_ROOT_PATH=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.02/x86_64-slc5-gcc43-opt/
source $XS_ROOT_PATH/root/bin/thisroot.sh

#Setup xrootd
export XS_XROOTD_PATH=/afs/cern.ch/sw/lcg/external/xrootd/3.1.0/x86_64-slc5-gcc46-opt/
export PATH=$XS_XROOTD_PATH/bin:$PATH
export LD_LIBRARY_PATH=$XS_XROOTD_PATH/lib64:$LD_LIBRARY_PATH

#Setup CMake
export XS_CMAKE_PATH=/afs/cern.ch/sw/lcg/external/CMake/2.8.6/x86_64-slc5-gcc46-opt/
export PATH=$XS_CMAKE_PATH/bin:$PATH

#Setup Python
export XS_PYTHON_PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.2/x86_64-slc5-gcc46-opt/
export PATH=$XS_PYTHON_PATH/bin:$PATH
export LD_LIBRARY_PATH=$XS_PYTHON_PATH/lib:$LD_LIBRARY_PATH

#Add python analysis packages
export XS_PYANALYSIS_PATH=/afs/cern.ch/sw/lcg/external/pyanalysis/1.3_python2.7/x86_64-slc5-gcc46-opt/
export PATH=$XS_PYANALYSIS_PATH/bin/:$PATH
export PYTHONPATH=$XS_PYANALYSIS_PATH/lib/python2.7/site-packages/:$PYTHONPATH

#Add python tools
export XS_PYTOOLS_PATH=/afs/cern.ch/sw/lcg/external/pytools/1.6_python2.7/x86_64-slc5-gcc46-opt/
export PATH=$XS_PYTOOLS_PATH/bin/:$PATH
export PYTHONPATH=$XS_PYTOOLS_PATH/lib/python2.7/site-packages/:$PYTHONPATH

#Setup Boost
export XS_BOOST_PATH=/afs/cern.ch/sw/lcg/external/Boost/1.48.0_python2.7/x86_64-slc5-gcc46-opt/
export BOOST_INCLUDEDIR=$XS_BOOST_PATH/include/boost-1_48/
export BOOST_LIBRARYDIR=$XS_BOOST_PATH/lib/
export LD_LIBRARY_PATH=$BOOST_LIBRARYDIR:$LD_LIBRARY_PATH

#Setup Valgrind
export XS_VALGRIND_PATH=/afs/cern.ch/sw/lcg/external/valgrind/3.6.0/x86_64-slc5-gcc46-opt/
export PATH=$XS_VALGRIND_PATH/bin:$PATH
export LD_LIBRARY_PATH=$XS_VALGRIND_PATH/lib:$LD_LIBRARY_PATH
