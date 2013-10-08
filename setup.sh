MELOC=$( dirname ${BASH_SOURCE} )

if [ "x${MELOC}" == "x" ]; then
    if [ "x$HEPPY" != "x" ]; then
	MELOC=$HEPPY/mymeanalysis
    else
	MELOC=$PWD
    fi
fi
MELOC=$( readlink -nf ${MELOC} )

echo ${MELOC}
export MELOC

export PYTHONPATH=$MELOC/exec/python:$PYTHONPATH
export LD_LIBRARY_PATH=$MELOC/shlib:$LD_LIBRARY_PATH
export PATH=$MELOC/exec/share:$PATH

if [ -d $MELOC/matrix ]; then
    mkdir -p $MELOC/matrix/objects
    mkdir -p $MELOC/integrator/objects
    mkdir -p $MELOC/shlib
    mkdir -p $MELOC/lib
fi

export PATH=$PATH:$MELOC/exec/run/grid:$MELOC/exec/run/prepare
export PATH=$PATH:$MELOC/exec/run/pbsjobs
