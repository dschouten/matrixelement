//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef WW1JFEYNDIAGRAM_HH
#define WW1JFEYNDIAGRAM_HH

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
// elements with final state W's (W->l,nu) + 2 jets (partons)
//
//////////////////////////////////////////////////////////////

class WW1jFeynDiagram
{
public:

  class TaggedJet
  {
  public:
    TaggedJet( );
    TaggedJet( const TLorentzVector& lv, bool tagged );
    TaggedJet( const TaggedJet& );
    TaggedJet operator=( const TaggedJet& );
    TLorentzVector lv;
    bool tagged;
    static bool sort(const TaggedJet& left, const TaggedJet& right);
  };

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

  TLorentzVector jet; // << outgoing partons

  TLorentzVector met; // << missing ET

  // ------------------- public methods ----------------------------

  WW1jFeynDiagram();

  WW1jFeynDiagram& operator=( const WW1jFeynDiagram& );

  void clear();

  // perform some basic parton-level calculations

  double sHat() const;
  double phaseSpace() const;
  double HT() const;
  void total( TLorentzVector& ) const;

  TLorentzVector total( ) const;
  TLorentzVector wplus() const;
  TLorentzVector wminus() const;
  
  // dump the event to stdout

  void print( ) const;
  
};

inline WW1jFeynDiagram::TaggedJet::TaggedJet( )
{ 
  tagged = false; 
  lv.SetPxPyPzE(0,0,0,0); 
}

inline WW1jFeynDiagram::TaggedJet::TaggedJet( const TLorentzVector& lv, bool tagged ) : 
  lv( lv ), tagged( tagged )
{ 

}

inline WW1jFeynDiagram::TaggedJet::TaggedJet( const TaggedJet& cl )
{
  tagged = cl.tagged;
  hepstd::tlvCopy( lv, cl.lv );
}

inline WW1jFeynDiagram::TaggedJet  WW1jFeynDiagram::TaggedJet::operator=( const TaggedJet& cl )
{
  tagged = cl.tagged;
  hepstd::tlvCopy( lv, cl.lv );
  return (*this);
}

inline TLorentzVector WW1jFeynDiagram::wminus() const 
{
  return lm + nur; 
}

inline TLorentzVector WW1jFeynDiagram::wplus() const 
{
  return lp + nul; 
}

inline void WW1jFeynDiagram::print( ) const 
{  
  std::cout << "qa:\t"  << pa.Px()  << "\t" << pa.Py()  << "\t" << pa.Pz()  << "\t" << pa.E()  << std::endl
            << "qb:\t"  << pb.Px()  << "\t" << pb.Py()  << "\t" << pb.Pz()  << "\t" << pb.E()  << std::endl
            << "l-:\t"  << lm.Px()  << "\t" << lm.Py()  << "\t" << lm.Pz()  << "\t" << lm.E()  << std::endl
            << "l+:\t"  << lp.Px()  << "\t" << lp.Py()  << "\t" << lp.Pz()  << "\t" << lp.E()  << std::endl
            << "vl:\t"  << nul.Px() << "\t" << nul.Py() << "\t" << nul.Pz() << "\t" << nul.E() << std::endl
            << "vr:\t"  << nur.Px() << "\t" << nur.Py() << "\t" << nur.Pz() << "\t" << nur.E() << std::endl
            << "j:\t"   << jet.Px() << "\t" << jet.Py() << "\t" << jet.Pz() << "\t" << jet.E() << std::endl; 
}

#endif
