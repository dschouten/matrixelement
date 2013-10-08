//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/WW2jFeynDiagram.hh"
#include "matrix/Constants.hh"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>
#include <TLorentzVector.h>

// ------------------------- ======= ------------------------- ======= -------------------------
WW2jFeynDiagram::WW2jFeynDiagram()
{
  ida = hepstd::kunknown;
  idb = hepstd::kunknown;

  clear();
}

// ------------------------- ======= ------------------------- ======= -------------------------
WW2jFeynDiagram& WW2jFeynDiagram::operator=( const WW2jFeynDiagram& clone )
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

    hepstd::tlvCopy( jeta, clone.jeta );

    hepstd::tlvCopy( jetb, clone.jetb );

    hepstd::tlvCopy( recoil, clone.recoil );

  } 
  return (*this);
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WW2jFeynDiagram::clear()
{
  lp.SetPxPyPzE(0, 0, 0, 0);
  lm.SetPxPyPzE(0, 0, 0, 0);

  nul.SetPxPyPzE(0, 0, 0, 0);
  nur.SetPxPyPzE(0, 0, 0, 0);

  pa.SetPxPyPzE(0, 0, 0, 0);
  pb.SetPxPyPzE(0, 0, 0, 0);

  jeta.SetPxPyPzE(0, 0, 0, 0);
  jetb.SetPxPyPzE(0, 0, 0, 0);

  recoil.SetPxPyPzE(0, 0, 0, 0);

  met.SetPxPyPzE(0, 0, 0, 0);
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WW2jFeynDiagram::phaseSpace() const
{
  double f = 64.0 * nul.E() * nur.E() * lp.E() * lm.E() * jeta.E() * jetb.E();
  if( f != 0 )
    return 1.0 / f;
  return 1.;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WW2jFeynDiagram::sHat() const
{
  return pow((pa + pb).M(),2);
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WW2jFeynDiagram::HT() const
{
  return (jeta.Pt() + jetb.Pt() + lp.Pt() + lm.Pt() + nul.Pt() + nur.Pt());
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WW2jFeynDiagram::total( TLorentzVector& total ) const
{
  total = lp + lm + nul + nur + jeta + jetb;
}
