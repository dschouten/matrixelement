//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef WFAKEFEYNDIAGRAM_HH
#define WFAKEFEYNDIAGRAM_HH

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <iomanip>
#include <TLorentzVector.h>

#include "matrix/Constants.hh"
#include "matrix/Utility.hh"

//////////////////////////////////////////////////////////////
//
// Class stores the parton information for calculating matrix
// elements with final state W's (W->l,nu)
//
//////////////////////////////////////////////////////////////

class WFakeFeynDiagram
{
public:

  // ------------------- event data ----------------------------

  double weight;

  TLorentzVector recoil; // << hadronic recoil

  std::vector<TLorentzVector> jets; 

  TLorentzVector fake; // << the fake lepton
  TLorentzVector real;
  TLorentzVector nu;
  TLorentzVector pa; // << incoming partons
  TLorentzVector pb;
  hepstd::PartonType ida;
  hepstd::PartonType idb;

  int fakecharge;
  int realcharge;

  TLorentzVector met; // << missing ET

  // ------------------- public methods ----------------------------
  
  WFakeFeynDiagram();

  WFakeFeynDiagram& operator=( const WFakeFeynDiagram& );

  void clear();

  // perform some basic parton-level calculations

  double sHat() const;
  double phaseSpace() const;
  double HT() const;
  void total( TLorentzVector& ) const;
  
  TLorentzVector wplus() const;
  TLorentzVector wminus() const;
  
  // dump the event to stdout

  void print( ) const;
  
};

inline TLorentzVector WFakeFeynDiagram::wminus() const 
{
  return realcharge > 0 ? TLorentzVector(0,0,0,0) : (real + nu);
}

inline TLorentzVector WFakeFeynDiagram::wplus() const 
{
  return realcharge > 0 ? (real + nu) : TLorentzVector(0,0,0,0);
}

inline void WFakeFeynDiagram::print( ) const 
{  
  std::cout << std::setprecision( 5 ) 	  
	<< "qa:\t"   << pa.Px()    << "\t" << pa.Py()   << "\t" << pa.Pz()   << "\t" << pa.E()   << std::endl            
	<< "qb:\t"   << pb.Px()    << "\t" << pb.Py()   << "\t" << pb.Pz()   << "\t" << pb.E()   << std::endl            
	<< "fake:\t" << fake.Px()  << "\t" << fake.Py() << "\t" << fake.Pz() << "\t" << fake.E() << std::endl          
	<< "l:\t"    << real.Px()  << "\t" << real.Py() << "\t" << real.Pz() << "\t" << real.E() << std::endl           
        << "v:\t"    << nu.Px()    << "\t" << nu.Py()   << "\t" << nu.Pz()   << "\t" << nu.E()   << std::endl;           
}
#endif
