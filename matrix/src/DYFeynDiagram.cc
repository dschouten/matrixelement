//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/DYFeynDiagram.hh"
#include "matrix/Constants.hh"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>
#include <TLorentzVector.h>

// ------------------------- ======= ------------------------- ======= -------------------------
DYFeynDiagram::DYFeynDiagram()
{
   ida = hepstd::kunknown;
   idb = hepstd::kunknown;

   clear();
}

// ------------------------- ======= ------------------------- ======= -------------------------
DYFeynDiagram& DYFeynDiagram::operator=( const DYFeynDiagram& clone )
{
  if( this != &clone )
  {
    hepstd::tlvCopy( pa , clone.pa );
    hepstd::tlvCopy( pb , clone.pb );
    hepstd::tlvCopy( lp , clone.lp );
    hepstd::tlvCopy( lm , clone.lm );
    hepstd::tlvCopy( jet, clone.jet );

    ida = clone.ida;
    idb = clone.idb;
  } 
  return (*this);
}

// ------------------------- ======= ------------------------- ======= -------------------------
void DYFeynDiagram::clear()
{
   lp.SetPxPyPzE(0, 0, 0, 0);
   lm.SetPxPyPzE(0, 0, 0, 0);

   pa.SetPxPyPzE(0, 0, 0, 0);
   pb.SetPxPyPzE(0, 0, 0, 0);

   jet.SetPxPyPzE(0, 0, 0, 0);
}


// ------------------------- ======= ------------------------- ======= -------------------------
double DYFeynDiagram::phaseSpace() const
{
  double f = 8 * jet.E() * lp.E() * lm.E();  
   if( f != 0 )
    return 1.0 / f;
  return 0.;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double DYFeynDiagram::sHat() const
{
  return pow((pa + pb).M(),2);
}

// ------------------------- ======= ------------------------- ======= -------------------------
double DYFeynDiagram::HT() const
{
  return (lp.Pt() + lm.Pt() + jet.Pt());
}

// ------------------------- ======= ------------------------- ======= -------------------------
void DYFeynDiagram::total( TLorentzVector& total ) const
{
  total = lp + lm + jet;
}
