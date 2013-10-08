//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/WW1jFeynDiagram.hh"
#include "matrix/Constants.hh"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>
#include <TLorentzVector.h>

// ------------------------- ======= ------------------------- ======= -------------------------
WW1jFeynDiagram::WW1jFeynDiagram()
{

  ida = hepstd::kunknown;
  idb = hepstd::kunknown;

  clear();
}

// ------------------------- ======= ------------------------- ======= -------------------------
WW1jFeynDiagram& WW1jFeynDiagram::operator=( const WW1jFeynDiagram& clone )
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

    hepstd::tlvCopy( jet, clone.jet );

    hepstd::tlvCopy( recoil, clone.recoil );

  } 
  return (*this);
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WW1jFeynDiagram::clear()
{
  lp.SetPxPyPzE(0, 0, 0, 0);
  lm.SetPxPyPzE(0, 0, 0, 0);

  nul.SetPxPyPzE(0, 0, 0, 0);
  nur.SetPxPyPzE(0, 0, 0, 0);

  pa.SetPxPyPzE(0, 0, 0, 0);
  pb.SetPxPyPzE(0, 0, 0, 0);

  jet.SetPxPyPzE(0, 0, 0, 0);

  recoil.SetPxPyPzE(0, 0, 0, 0);

  met.SetPxPyPzE(0, 0, 0, 0);
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WW1jFeynDiagram::phaseSpace() const
{
  double f = 32.0 * nul.E() * nur.E() * lp.E() * lm.E() * jet.E();
  if( f != 0 )
    return 1.0 / f;
  return 1.;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WW1jFeynDiagram::sHat() const
{
  return pow((pa + pb).M(),2);
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WW1jFeynDiagram::HT() const
{
  return (jet.Pt() + lp.Pt() + lm.Pt() + nul.Pt() + nur.Pt());
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WW1jFeynDiagram::total( TLorentzVector& total ) const
{
  total = lp + lm + nul + nur + jet;
}

// ------------------------- ======= ------------------------- ======= -------------------------
TLorentzVector WW1jFeynDiagram::total( ) const
{
  return (lp + lm + nul + nur + jet);
}
