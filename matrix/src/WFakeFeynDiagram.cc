//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/WFakeFeynDiagram.hh"
#include "matrix/Constants.hh"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>
#include <TLorentzVector.h>

// ------------------------- ======= ------------------------- ======= -------------------------
WFakeFeynDiagram::WFakeFeynDiagram()
{
  ida = hepstd::kunknown;
  idb = hepstd::kunknown;
  
  clear();
}

// ------------------------- ======= ------------------------- ======= -------------------------
WFakeFeynDiagram& WFakeFeynDiagram::operator=( const WFakeFeynDiagram& clone )
{
  if( this != &clone )
  {
    hepstd::tlvCopy( pa, clone.pa );
    hepstd::tlvCopy( pb, clone.pb );
    hepstd::tlvCopy( fake, clone.fake );
    hepstd::tlvCopy( real, clone.real );
    
    hepstd::tlvCopy( nu, clone.nu );
    
    hepstd::tlvCopy( met, clone.met );
    
    ida = clone.ida;
    idb = clone.idb;
    
    fakecharge = clone.fakecharge;
    realcharge = clone.realcharge;

    hepstd::tlvCopy( recoil, clone.recoil );
  } 
  return (*this);
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WFakeFeynDiagram::clear()
{
  fake.SetPxPyPzE(0, 0, 0, 0);
  real.SetPxPyPzE(0, 0, 0, 0);
  
  nu.SetPxPyPzE(0, 0, 0, 0);
  
  pa.SetPxPyPzE(0, 0, 0, 0);
  pb.SetPxPyPzE(0, 0, 0, 0);
  
  recoil.SetPxPyPzE(0, 0, 0, 0);
  
  met.SetPxPyPzE(0, 0, 0, 0);
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WFakeFeynDiagram::phaseSpace() const
{
  double f = 8.0 * nu.E() * fake.E() * real.E();  
  // for( unsigned int i = 0; i < jets.size(); ++i ) 
  // {
  //   f *= 2.0 * jets[i].E();
  // }
  if( f != 0 )
    return 1.0 / f;
  return 0.;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WFakeFeynDiagram::sHat() const
{
  return pow((pa + pb).M(),2);
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WFakeFeynDiagram::HT() const
{
  return (fake.Pt() + real.Pt() + nu.Pt());
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WFakeFeynDiagram::total( TLorentzVector& total ) const
{
  total = fake + real + nu;
}
