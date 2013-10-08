//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//      
//

#ifndef DY2JFEYNDIAGRAM_HH
#define DY2JFEYNDIAGRAM_HH

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
// elements with Z/gamma* -> leptons and two jets
//
//////////////////////////////////////////////////////////////

class DY2jFeynDiagram
{
public:

  // ------------------- event data ----------------------------

  TLorentzVector jeta;
  TLorentzVector jetb;
  TLorentzVector lp; // << outgoing leptons
  TLorentzVector lm;
  TLorentzVector pa; // << incoming partons
  TLorentzVector pb;
  hepstd::PartonType ida;
  hepstd::PartonType idb;

  // ------------------- public methods ----------------------------

  DY2jFeynDiagram();

  DY2jFeynDiagram& operator=( const DY2jFeynDiagram& );

  void clear();

  // perform some basic parton-level calculations

  double sHat() const;
  double phaseSpace() const;
  double HT() const;
  void total( TLorentzVector& ) const;
  TLorentzVector total( ) const;

  TLorentzVector zgamma() const;
  
  // dump the event to stdout

  void print( ) const;
  
};

inline TLorentzVector DY2jFeynDiagram::zgamma() const 
{
  return lm + lp; 
}

inline void DY2jFeynDiagram::print( ) const 
{  
  std::cout << "qa:\t"  << pa.Px()  << "\t" << pa.Py()  << "\t" << pa.Pz()  << "\t" << pa.E()  << std::endl            
	    << "qb:\t"  << pb.Px()  << "\t" << pb.Py()  << "\t" << pb.Pz()  << "\t" << pb.E()  << std::endl            
	    << "l-:\t"  << lm.Px()  << "\t" << lm.Py()  << "\t" << lm.Pz()  << "\t" << lm.E()  << std::endl          
	    << "l+:\t"  << lp.Px()  << "\t" << lp.Py()  << "\t" << lp.Pz()  << "\t" << lp.E()  << std::endl           
	    << "j :\t"  << jeta.Px() << "\t" << jeta.Py() << "\t" << jeta.Pz() << "\t" << jeta.E() << std::endl
	    << "j :\t"  << jetb.Px() << "\t" << jetb.Py() << "\t" << jetb.Pz() << "\t" << jetb.E() << std::endl; 
}
#endif
