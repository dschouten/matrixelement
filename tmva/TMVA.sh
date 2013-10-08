## rm -rf TMVA-*/
## wget http://sourceforge.net/projects/tmva/files/latest/download
## wget http://trshare.triumf.ca/~${BASE}/TMVA-v4.1.2patch.tgz
## tar xzf TMVA-*.tgz

cd TMVA-*/
make -j

rm -rf ../../shlib/libTMVA.so
ln -s $PWD/lib/libTMVA.so ../../shlib/libTMVA.so

ln -s test macros

cd ../

rm -rf TMVA
ln -s TMVA-*/TMVA TMVA

root -l -q -b 'compile.C'
mv TMVAWrapper_C.so ../shlib/libTMVAWrapper.so
rm TMVAWrapper_C.d
