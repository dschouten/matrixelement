//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/METResolutionTF.hh"

#include <iostream>
#include <string>
#include <math.h>
#include <stdexcept>

// ------------------------- ======= ------------------------- ======= -------------------------
METResolutionTF::METResolutionTF( const std::string& defn ) : 
  TransferFunction( defn ),
  m_sigma_x_scale( -1 ),
  m_sigma_y_scale( -1 )
{
  m_sigma_x_scale = parameter("sigma_x_scale");
  m_sigma_y_scale = parameter("sigma_y_scale");
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool METResolutionTF::limits( const char* name, const double& val, const double par[], double vlim[] ) 
{
  double sigma = -1;
  if( name[0] == 'x' ) 
  {
    sigma = m_sigma_x_scale * sqrt( par[0] );
  }
  else
  {
    if( name[0] == 'y' )
    {
      sigma = m_sigma_y_scale * sqrt( par[0] );
    }
    else
    {
      std::cout << "WARNING variable [" << name << "] not known." << std::endl;
      return false;
    }
  }
  vlim[0] = val - 3 * sigma;
  vlim[1] = val + 3 * sigma;
  return true;
}
  
// ------------------------- ======= ------------------------- ======= -------------------------
double METResolutionTF::operator()( const char* name, const double& x, const double par[] ) 
{
  static double sigma = -1;
  if( name[0] == 'x' ) 
  {
    sigma = m_sigma_x_scale * sqrt( par[0] );
  }
  else
  {
    if( name[0] == 'y' )
    {
      sigma = m_sigma_y_scale * sqrt( par[0] );
    }
  }  
  if( fabs(x - par[1]) > 4.0 * sigma )
  {
    SHOW_DEBUG( std::cout << "sigma: " << sigma << ", (x - mu) = " << x - par[1] << ", TF = " 
		<< 1.0 / ( M_SQRT_2PI * sigma ) * exp( -((x - par[1])*(x - par[1])) / (2 * sigma * sigma) )
		<< std::endl );
    throw std::runtime_error( "variable is outside the bounds of the TF!" );
  }
  double arg = -((x - par[1])*(x - par[1])) / (2.0 * sigma * sigma);
  return 1.0 / ( M_SQRT_2PI * sigma ) * hepstd::fexp( arg );
}
