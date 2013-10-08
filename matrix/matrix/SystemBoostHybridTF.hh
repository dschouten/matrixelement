//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef SYSTEMBOOSTHYBRID_TF
#define SYSTEMBOOSTHYBRID_TF

#include "matrix/TransferFunction.hh"
#include "matrix/Utility.hh"

#include <cmath>

class TH1;
class TF1;
class TSpline;

//////////////////////////////////////////////////////////
//
// Transfer function for system transverse boost
//
// - provides probability distribution for pT from theory input
//   but uses phi estimate from recoil measurement
//
// phi ~ n(0,sigma) + offset
// pT is sampled from theory spectrum
//
// Version 1.0
//
//////////////////////////////////////////////////////////

class SystemBoostHybridTF : public TransferFunction
{
public:
  SystemBoostHybridTF( );
  SystemBoostHybridTF( const std::string& defn );
  virtual ~SystemBoostHybridTF();

  virtual bool limits( const char*, const double&, const double[], double[] );
  virtual double operator()( const char*, const double&, const double[] );

private:
  mutable TF1* m_theory_fit;
  mutable TSpline* m_theory_spl;
  mutable TH1* m_theory_dat; //< theory spectrum for system pT
  mutable TF1* m_sigma_phi_fun; ///< function for scaling phi resolution with pT

  double m_a, m_b; //< limits of the system pT spectrum
  
  mutable TF1* m_delta_phi_fun; ///< function for defining offset term for phi TF

  double phi_integral( double sigma, double delta ) const;
  double phi_kernel( double x, double sigma, double delta ) const;
};

//------------------------- ======= ------------------------- ======= -------------------------
inline double SystemBoostHybridTF::phi_integral( double sigma, double delta ) const
{
  return hepstd::erf( M_PI / (M_SQRT_2 * sigma) ) + delta; // << \int[-pi,pi] phi_kernel
}

//------------------------- ======= ------------------------- ======= -------------------------
inline double SystemBoostHybridTF::phi_kernel( double dx, double sigma, double delta ) const
{
  return 1. / ( M_SQRT_2PI * sigma ) * std::exp( -((dx)*(dx)) / (2.0 * sigma * sigma) ) + delta / (2*M_PI); // phi_kernel ~ n(0,sigma) + delta
}

#endif
