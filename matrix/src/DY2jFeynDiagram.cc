//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/DY2jFeynDiagram.hh"
#include "matrix/Constants.hh"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>
#include <TLorentzVector.h>

// ------------------------- ======= ------------------------- ======= -------------------------
DY2jFeynDiagram::DY2jFeynDiagram()
{
   ida = hepstd::kunknown;
   idb = hepstd::kunknown;

   clear();
}

// ------------------------- ======= ------------------------- ======= -------------------------
DY2jFeynDiagram& DY2jFeynDiagram::operator=( const DY2jFeynDiagram& clone )
{
  if( this != &clone )
  {
    hepstd::tlvCopy( pa , clone.pa );
    hepstd::tlvCopy( pb , clone.pb );
    hepstd::tlvCopy( lp , clone.lp );
    hepstd::tlvCopy( lm , clone.lm );
    hepstd::tlvCopy( jeta, clone.jeta );
    hepstd::tlvCopy( jetb, clone.jetb );

    ida = clone.ida;
    idb = clone.idb;
  } 
  return (*this);
}

// ------------------------- ======= ------------------------- ======= -------------------------
void DY2jFeynDiagram::clear()
{
   lp.SetPxPyPzE(0, 0, 0, 0);
   lm.SetPxPyPzE(0, 0, 0, 0);

   pa.SetPxPyPzE(0, 0, 0, 0);
   pb.SetPxPyPzE(0, 0, 0, 0);

   jeta.SetPxPyPzE(0, 0, 0, 0);
   jetb.SetPxPyPzE(0, 0, 0, 0);
}


// ------------------------- ======= ------------------------- ======= -------------------------
double DY2jFeynDiagram::phaseSpace() const
{
  double f = 16 * jeta.E() * jetb.E() * lp.E() * lm.E();  
   if( f != 0 )
    return 1.0 / f;
  return 0.;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double DY2jFeynDiagram::sHat() const
{
  return pow((pa + pb).M(),2);
}

// ------------------------- ======= ------------------------- ======= -------------------------
double DY2jFeynDiagram::HT() const
{
  return (lp.Pt() + lm.Pt() + jeta.Pt() + jetb.Pt());
}

// ------------------------- ======= ------------------------- ======= -------------------------
void DY2jFeynDiagram::total( TLorentzVector& total ) const
{
  total = lp + lm + jeta + jetb;
}

// ------------------------- ======= ------------------------- ======= -------------------------
TLorentzVector DY2jFeynDiagram::total( ) const
{
  return lp + lm + jeta + jetb;
}
