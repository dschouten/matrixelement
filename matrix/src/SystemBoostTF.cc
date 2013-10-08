//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/SystemBoostTF.hh"

#include <iostream>
#include <string>
#include <math.h>
#include <stdexcept>

#include <TF1.h>

/**
 * 
 * Read in functions that parametrize the pT an phi resolution as a function of the
 * estimated recoil pT. The TF's for both pT and phi are assumed to be N(0,sigma)
 *
 * notes/caveats: 
 *    - the functions (of 'x') must be specified without spaces in the .tf file
 *    - the integration limits for phi are reported in terms of +/- dphi instead 
 *      of absolute coordinates to avoid unnecessarily complex code to deal
 *      with phi-wrapping
 *    - the TF's for pT and phi are accessed with names "pT" and "azimuth"
 *
 **/

// ------------------------- ======= ------------------------- ======= -------------------------
SystemBoostTF::SystemBoostTF( const std::string& defn ) : 
  TransferFunction( defn ),
  m_sigma_mom_fun( 0x0 ),
  m_sigma_phi_fun( 0x0 ),
  m_use_phi_range( false )
{
  m_sigma_mom_fun = new TF1( "SIGMA_MOM", str_parameter( "sigma_mom" ).c_str() );
  m_mu_mom_fun    = new TF1( "MU_MOM",    str_parameter( "mean_mom" ).c_str() );
  m_sigma_phi_fun = new TF1( "SIGMA_PHI", str_parameter( "sigma_phi" ).c_str() );
  m_delta_phi_fun = new TF1( "DELTA_PHI", str_parameter( "offset_phi" ).c_str() );

  str_parameter( "use_phi_range", m_use_phi_range );
}

// ------------------------- ======= ------------------------- ======= -------------------------
SystemBoostTF::~SystemBoostTF()
{
  if( m_sigma_mom_fun ) delete m_sigma_mom_fun;
  if( m_mu_mom_fun ) delete m_mu_mom_fun;
  if( m_sigma_phi_fun ) delete m_sigma_phi_fun;
  if( m_delta_phi_fun ) delete m_delta_phi_fun;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool SystemBoostTF::limits( const char* name, const double& val, const double par[], double vlim[] ) 
{
  double sigma, mu, delta;
  switch( name[0] )
  {
    case 'p':
      sigma = fabs(m_sigma_mom_fun->Eval( par[0] )); 
      mu = fabs(m_mu_mom_fun->Eval( par[0] ));
      vlim[0] = std::max( 0., mu - 4.0 * sigma );
      vlim[1] = std::max( 0., mu + 4.0 * sigma );
      break;
    case 'a':
      sigma = fabs(m_sigma_phi_fun->Eval( par[0] )); 
      delta = fabs(m_delta_phi_fun->Eval( par[0] )); 
      vlim[0] = !m_use_phi_range ? -M_PI : -std::min( M_PI, 4.0 * sigma )*(delta==0) - M_PI*(delta!=0); // if delta != 0, integrate over whole phase space
      vlim[1] = !m_use_phi_range ?  M_PI :  std::min( M_PI, 4.0 * sigma )*(delta==0) + M_PI*(delta!=0);
      break;
    default:
      std::cout << "WARNING variable [" << name << "] not known." << std::endl;
      return false;
  }
  return true;
}
  
// ------------------------- ======= ------------------------- ======= -------------------------
double SystemBoostTF::operator()( const char* name, const double& x, const double par[] ) 
{

  double sigma, delta, mu, dx;

  double retval = 1;

  switch( name[0] )
  {
    case 'p':
      sigma	= fabs(m_sigma_mom_fun->Eval( par[0] )); // @todo: cache the TF1::Eval's since par[0] changes only per event
      mu	= fabs(m_mu_mom_fun->Eval( par[0] ));
      dx        = fabs( x - mu );
      retval    = mom_kernel( dx, sigma )*((double)(x>0)) / mom_integral( mu, sigma ); // mom_integral != 1 for mu < 4sigma since pT is > 0
      break;
    case 'a':
      sigma	= fabs(m_sigma_phi_fun->Eval( par[0] ));
      delta	= fabs(m_delta_phi_fun->Eval( par[0] ));
      dx        = fabs(hepstd::phiMPiPi( x - par[1] ));
      retval    = phi_kernel( dx, sigma, delta ) / phi_integral( sigma, delta ); // phi_integral accounts for offset term, and for sigma > pi/4
      break;
    default:
      std::cout << "WARNING variable [" << name << "] not known." << std::endl;
      return 0;
  }

  return retval;
}
