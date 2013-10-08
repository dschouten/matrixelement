//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/WW0jFeynDiagram.hh"
#include "matrix/Constants.hh"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>
#include <TLorentzVector.h>

// ------------------------- ======= ------------------------- ======= -------------------------
WW0jFeynDiagram::WW0jFeynDiagram()
{
  ida = hepstd::kunknown;
  idb = hepstd::kunknown;
  
  clear();
}

// ------------------------- ======= ------------------------- ======= -------------------------
WW0jFeynDiagram& WW0jFeynDiagram::operator=( const WW0jFeynDiagram& clone )
{
  if( this != &clone )
  {
    hepstd::tlvCopy( pa, clone.pa );
    hepstd::tlvCopy( pb, clone.pb );
    hepstd::tlvCopy( lp, clone.lp );
    hepstd::tlvCopy( lm, clone.lm );
    
    hepstd::tlvCopy( nul, clone.nul );
    hepstd::tlvCopy( nur, clone.nur );
    
    hepstd::tlvCopy( met, clone.met );
    
    ida = clone.ida;
    idb = clone.idb;
    
    hepstd::tlvCopy( recoil, clone.recoil );
  } 
  return (*this);
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WW0jFeynDiagram::clear()
{
  lp.SetPxPyPzE(0, 0, 0, 0);
  lm.SetPxPyPzE(0, 0, 0, 0);
  
  nul.SetPxPyPzE(0, 0, 0, 0);
  nur.SetPxPyPzE(0, 0, 0, 0);
  
  pa.SetPxPyPzE(0, 0, 0, 0);
  pb.SetPxPyPzE(0, 0, 0, 0);
  
  recoil.SetPxPyPzE(0, 0, 0, 0);
  
  met.SetPxPyPzE(0, 0, 0, 0);
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WW0jFeynDiagram::phaseSpace() const
{
  double f = 16.0 * nul.E() * nur.E() * lp.E() * lm.E();  
  if( f != 0 )
    return 1.0 / f;
  return 0.;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WW0jFeynDiagram::sHat() const
{
  return pow((pa + pb).M(),2);
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WW0jFeynDiagram::HT() const
{
  return (lp.Pt() + lm.Pt() + nul.Pt() + nur.Pt());
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WW0jFeynDiagram::total( TLorentzVector& total ) const
{
  total = lp + lm + nul + nur;
}

// ------------------------- ======= ------------------------- ======= -------------------------
TLorentzVector WW0jFeynDiagram::total( ) const
{
  return ( lp + lm + nul + nur );
}
