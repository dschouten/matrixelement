//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef WW0JFEYNDIAGRAM_HH
#define WW0JFEYNDIAGRAM_HH

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

class WW0jFeynDiagram
{
public:

  // ------------------- event data ----------------------------

  double weight;

  TLorentzVector recoil; // << hadronic recoil

  TLorentzVector lp; // << outgoing leptons
  TLorentzVector nur;
  TLorentzVector lm;
  TLorentzVector nul;
  TLorentzVector pa; // << incoming partons
  TLorentzVector pb;
  hepstd::PartonType ida;
  hepstd::PartonType idb;
  TLorentzVector met; // << missing ET

  // ------------------- public methods ----------------------------
  
  WW0jFeynDiagram();

  WW0jFeynDiagram& operator=( const WW0jFeynDiagram& );

  void clear();

  // perform some basic parton-level calculations

  double sHat() const;
  double phaseSpace() const;
  double HT() const;
  TLorentzVector total(  ) const;
  void total( TLorentzVector& ) const;

  TLorentzVector wplus() const;
  TLorentzVector wminus() const;
  
  // dump the event to stdout

  void print( ) const;
  
};

inline TLorentzVector WW0jFeynDiagram::wminus() const 
{
  return lm + nur; 
}

inline TLorentzVector WW0jFeynDiagram::wplus() const 
{
  return lp + nul; 
}

inline void WW0jFeynDiagram::print( ) const 
{  
  std::cout << "qa:\t"  << pa.Px()  << "\t" << pa.Py()  << "\t" << pa.Pz()  << "\t" << pa.E()  << std::endl            
	    << "qb:\t"  << pb.Px()  << "\t" << pb.Py()  << "\t" << pb.Pz()  << "\t" << pb.E()  << std::endl            
	    << "l-:\t"  << lm.Px()  << "\t" << lm.Py()  << "\t" << lm.Pz()  << "\t" << lm.E()  << std::endl          
	    << "l+:\t"  << lp.Px()  << "\t" << lp.Py()  << "\t" << lp.Pz()  << "\t" << lp.E()  << std::endl           
	    << "vr:\t"  << nul.Px() << "\t" << nul.Py() << "\t" << nul.Pz() << "\t" << nul.E() << std::endl           
	    << "vl:\t"  << nur.Px() << "\t" << nur.Py() << "\t" << nur.Pz() << "\t" << nur.E() << std::endl; 
}
#endif
