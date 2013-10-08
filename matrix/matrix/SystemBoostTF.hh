//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//         
//

#ifndef SYSTEMBOOST_TF
#define SYSTEMBOOST_TF

#include "matrix/TransferFunction.hh"
#include "matrix/Utility.hh"

class TF1;

//////////////////////////////////////////////////////////
//
// Transfer function for system transverse boost
//
// - provides probability distribution for pT, \phi given 
//   measured recoil r
//
// Version 2.0
// - momentum and azimuthal resolution are f( pT )
//   pT  ~ n(mu,sigma), mu = mu(rT), sigma = sigma(rT)
//   phi ~ n(0,sigma) + offset, sigma = sigma(rT), offset = delta(rT)
//
//////////////////////////////////////////////////////////

class SystemBoostTF : public TransferFunction
{
public:
  SystemBoostTF( const std::string& defn );
  virtual ~SystemBoostTF();

  virtual bool limits( const char*, const double&, const double[], double[] );
  virtual double operator()( const char*, const double&, const double[] );

private:
  mutable TF1* m_sigma_mom_fun; ///< function for scaling momentum resolution with rT
  mutable TF1* m_mu_mom_fun; ///< function for setting mean pT with rT
  mutable TF1* m_sigma_phi_fun; ///< ... and for angular resolution 
  mutable TF1* m_delta_phi_fun; ///< ... and for angular offset term

  bool m_use_phi_range; ///< flag to use phi resolution to set region of integration

  double phi_integral( double sigma, double delta ) const;
  double phi_kernel( double x, double sigma, double delta ) const;

  double mom_integral( double mu, double sigma ) const;
  double mom_kernel( double dx, double sigma ) const;
};

//------------------------- ======= ------------------------- ======= -------------------------
inline double SystemBoostTF::phi_integral( double sigma, double delta ) const
{
  return hepstd::erf( M_PI / (M_SQRT_2 * sigma) ) + delta; // << \int[-pi,pi] phi_kernel
}

//------------------------- ======= ------------------------- ======= -------------------------
inline double SystemBoostTF::phi_kernel( double dx, double sigma, double delta ) const
{
  return 1. / ( M_SQRT_2PI * sigma ) * std::exp( -((dx)*(dx)) / (2.0 * sigma * sigma) ) + delta / (2*M_PI); // phi_kernel ~ n(0,sigma) + delta
}

//------------------------- ======= ------------------------- ======= -------------------------
inline double SystemBoostTF::mom_integral( double mu, double sigma ) const
{
  double den = (M_SQRT_2 * sigma);
  return 0.5 * ( hepstd::erf( 4*sigma / den ) -
		 hepstd::erf( (std::max(0.0,mu-4*sigma) - mu) / den ) ); // << \int[max(0,mu-4sigma),mu+4sigma)] mom_kernel
}

//------------------------- ======= ------------------------- ======= -------------------------
inline double SystemBoostTF::mom_kernel( double dx, double sigma ) const
{
  return 1. / ( M_SQRT_2PI * sigma ) * std::exp( -((dx)*(dx)) / (2.0 * sigma * sigma) ); // mom_kernel ~ n(mu,sigma)
}

#endif

  
  
