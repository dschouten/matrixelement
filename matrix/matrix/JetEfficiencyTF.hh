//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef JETEFFICIENCY_TF
#define JETEFFICIENCY_TF

#include "matrix/TransferFunction.hh"

//////////////////////////////////////////////////////////
//
// Transfer function for jets 
//
// - NOTE provides probability for parton pT given no observed
//   jet (1 - efficiency), i.e. efficiency to NOT be measured
//
// Version 1.0
// - threshold function provided for single pseudorapidity bin 
//
// Version 2.0
// - binned in rapidity, using function
//    epsilon(pT) = b + (a - b) / (1 + exp( pT^{d} - x_{0} ) / c) + e
// 
//   parameters are listed for rapidity bins in order: x_{0},a,b,c,d,e
// 
//////////////////////////////////////////////////////////

class JetEfficiencyTF : public TransferFunction
{
public:
  JetEfficiencyTF( const std::string& defn );
  virtual ~JetEfficiencyTF() { }

  // return limits [a,b] such that 1 - epsilon < 1.0e-3
  // parameters
  //    char : 'eff' -> efficiency (error function)
  //    double : pT -> pT of the jet
  //    double [] : parameters
  //    double [] : limits [a,b]
  virtual bool limits( const char*, const double&, const double[], double[] );
  virtual double operator()( const char*, const double&, const double[] );

private:
  std::vector<double> m_par;
  std::vector<double> m_eta_bins;

  unsigned m_nbins;

  unsigned m_version; //<< version of the TF defined

  double m_lower_bound;
  double m_upper_bound;
};

#endif
