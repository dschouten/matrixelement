{
  gSystem->Load( "libCore.so" );
  gSystem->Load( "libTreePlayer.so" );
  gSystem->Load( "libTMVA.so" );
  gROOT->ProcessLine( ".L TMVAWrapper.C++" );
}
