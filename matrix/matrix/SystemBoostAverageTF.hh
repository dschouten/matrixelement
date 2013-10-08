//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef SYSTEMBOOSTAVERAGE_TF
#define SYSTEMBOOSTAVERAGE_TF

#include "matrix/TransferFunction.hh"

class TH1;
class TF1;
class TSpline;

//////////////////////////////////////////////////////////
//
// Transfer function for system transverse boost
//
// - provides probability distribution for pT from theory input
//
// Version 1.0
//
//////////////////////////////////////////////////////////

class SystemBoostAverageTF : public TransferFunction
{
public:
  SystemBoostAverageTF( );
  SystemBoostAverageTF( const std::string& defn );
  virtual ~SystemBoostAverageTF();

  virtual bool limits( const char*, const double&, const double[], double[] );
  virtual double operator()( const char*, const double&, const double[] );

private:
  mutable TF1* m_theory_fit;
  mutable TSpline* m_theory_spl;

  mutable TH1* m_theory_dat; //< theory spectrum for system pT

  double m_a, m_b; //< limits of the system pT spectrum
};

#endif
