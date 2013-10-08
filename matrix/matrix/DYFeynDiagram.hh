//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//     
//

#ifndef DYJFEYNDIAGRAM_HH
#define DYJFEYNDIAGRAM_HH

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
// elements with Z/gamma* -> leptons and one jet
//
//////////////////////////////////////////////////////////////

class DYFeynDiagram
{
public:

  // ------------------- event data ----------------------------

  TLorentzVector jet; // << hadronic recoil/jet
  TLorentzVector lp; // << outgoing leptons
  TLorentzVector lm;
  TLorentzVector pa; // << incoming partons
  TLorentzVector pb;
  hepstd::PartonType ida;
  hepstd::PartonType idb;

  // ------------------- public methods ----------------------------

  DYFeynDiagram();

  DYFeynDiagram& operator=( const DYFeynDiagram& );

  void clear();

  // perform some basic parton-level calculations

  double sHat() const;
  double phaseSpace() const;
  double HT() const;
  void total( TLorentzVector& ) const;

  TLorentzVector zgamma() const;
  
  // dump the event to stdout

  void print( ) const;
  
};

inline TLorentzVector DYFeynDiagram::zgamma() const 
{
  return lm + lp; 
}

inline void DYFeynDiagram::print( ) const 
{  
  std::cout << "qa:\t"  << pa.Px()  << "\t" << pa.Py()  << "\t" << pa.Pz()  << "\t" << pa.E()  << std::endl            
	    << "qb:\t"  << pb.Px()  << "\t" << pb.Py()  << "\t" << pb.Pz()  << "\t" << pb.E()  << std::endl            
	    << "l-:\t"  << lm.Px()  << "\t" << lm.Py()  << "\t" << lm.Pz()  << "\t" << lm.E()  << std::endl          
	    << "l+:\t"  << lp.Px()  << "\t" << lp.Py()  << "\t" << lp.Pz()  << "\t" << lp.E()  << std::endl           
	    << "j :\t"  << jet.Px() << "\t" << jet.Py() << "\t" << jet.Pz() << "\t" << jet.E() << std::endl; 
}
#endif
