if [ "x$HEPPY" != "x" ]; then
    MELOC=$HEPPY/mymeanalysis.git
else
    MELOC=$PWD
fi

if [ "x${BASE}" == "x" ]; then
    BASE="/tmp"
fi

cd ${BASE}

rm -rf mymeanalysis.sandbox
git clone $MELOC mymeanalysis.sandbox

cd mymeanalysis.sandbox

git log -n 1 --oneline | cut -d" " -f1 > version.git
cp ${MELOC}/exec/run/grid/*.sh .
cp ${MELOC}/exec/run/grid/Makefile.submit .
touch job.cfg

if [ "x${PWD}" == "x${BASE}/mymeanalysis.sandbox" ]; then 
    rm -rf .git/
fi

cd $( dirname $MELOC )
