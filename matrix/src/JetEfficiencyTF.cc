//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/JetEfficiencyTF.hh"

#include <iostream>
#include <string>
#include <math.h>
#include <stdexcept>

#include "TF1.h"

#define VERSION_1 1
#define VERSION_2 2

namespace
{
  class Functor
  {
  public:
    Functor( JetEfficiencyTF* tf_ptr ) : m_tf_ptr( tf_ptr ) { }
    
    double operator()( double* x, double* par )
    {
      if( m_tf_ptr == NULL )
      {
	return 0.;
      }
      return (*m_tf_ptr)( "eff", (*x), par );
    }
    
    JetEfficiencyTF* m_tf_ptr;
  };  
  
  class FunctionWrapper : public TF1
  {
  public:
    FunctionWrapper( Functor* f, double l, double h ) : TF1( "f", f, l, h, 1 ) { }
    double operator()( double x ) { return Eval( x ); }
  };
}

#define NUM_PAR 6
#define MAX_ETA 4.5

// ------------------------- ======= ------------------------- ======= -------------------------
JetEfficiencyTF::JetEfficiencyTF( const std::string& defn ) : 
  TransferFunction( defn )
{
  m_par = std::vector<double>( NUM_PAR, 0. );

  bool flag = false;
  m_par[0] = parameter( "erf_threshold", flag );
  m_par[1] = parameter( "erf_width", flag );
  
  m_version = VERSION_1;

  if( !flag ) // didn't find parameters for VERSION_1 implementation
  {
    m_version = VERSION_2;

    m_nbins = (unsigned) parameter( "num_eta_bins", flag );    
    if( !flag )
    {
      std::cout << "ERROR could not read parameter [num_eta_bins] for JetEfficiencyTF!" << std::endl;
      throw std::runtime_error("could not read parameter");
    }    

    char par[16]; 
    
    for( unsigned i = 0; i < m_nbins; ++i )
    {
      sprintf( par, "eta_bin_%d", i );
      m_eta_bins.push_back( parameter( par, flag ) );
      if( !flag )
      {
	std::cout << "ERROR could not read parameter [" << par << "] for JetEfficiencyTF!" << std::endl;
	throw std::runtime_error("could not read parameter");
      }
    }
    m_eta_bins.push_back( MAX_ETA ); // read in the rapidity bins

    m_par = std::vector<double>( NUM_PAR * m_nbins, 0. );
 
    for( unsigned i = 0; i < m_par.size(); ++i )
    {
      sprintf( par, "par_%d", i );
      m_par[i] = parameter( par, flag );
      if( !flag )
      {
	std::cout << "ERROR could not read parameter [" << par << "] for JetEfficiencyTF!" << std::endl;
	throw std::runtime_error("could not read parameter");
      }
    } // read the function parameters for each rapidity bin
  }

  m_lower_bound = parameter( "lower_bound" );
  m_lower_bound = m_lower_bound <= 0 ? 0 : m_lower_bound;
  
  m_upper_bound = parameter( "upper_bound" );
  m_upper_bound = m_upper_bound <= 0 ? 50 : m_upper_bound;

  SHOW_DEBUG( std::cout << "Initialized JetEfficiency TF v" << m_version << std::endl );  
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool JetEfficiencyTF::limits( const char* name, const double& val, const double par[], double vlim[] ) 
{
  if( name[0] == 'e' )
  {
    vlim[0] = m_lower_bound;
    vlim[1] = m_upper_bound;
    return true;
  }
  return false;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double JetEfficiencyTF::operator()( const char* name, const double& x, const double par[] ) 
{
  if( name[0] == 'e' )
  {

    if( x > m_upper_bound )
      return 0;

    if( x < m_lower_bound )
      return (*this)( name, m_lower_bound, par ); // cut-off at lower bound for b-quark E

    if( m_version == VERSION_1 )
    {
      double arg = (x - m_par[0]) / (M_SQRT_2 * m_par[1]);
      return 0.5 * (1.0 - hepstd::erf( arg )); // probability to not observe a parton with erf(x) jet turn-on 
    }

    if( m_version == VERSION_2 )
    {
      if( fabs( par[0] ) >= MAX_ETA )
      {
	return 1.0; // don't include the phase space outside of TF fiducial volume
      }
      
      unsigned ibin = 0;

      while( ibin < m_nbins - 1 )
      {
	if( m_eta_bins[ibin] < fabs( par[0] ) &&
	    m_eta_bins[ibin + 1] >= fabs( par[0] ) )
	{
	  break;
	}
	ibin += 1;
      }
	    
      unsigned ipar = ibin * NUM_PAR;

      unsigned ipar_0 = ipar;
      unsigned ipar_1 = ipar + 1;
      unsigned ipar_2 = ipar + 2;
      unsigned ipar_3 = ipar + 3;
      unsigned ipar_4 = ipar + 4;
      unsigned ipar_5 = ipar + 5;

      if( ipar + NUM_PAR > m_par.size() )
      {
	return 1.0;
      }      
      return fabs( std::max( 1.0 - ( m_par[ipar_2] + 
			       ( (m_par[ipar_1]-m_par[ipar_2]) / 
				 (1.0 + hepstd::fexp((pow(x,m_par[ipar_4])-m_par[ipar_0])/m_par[ipar_3])+m_par[ipar_5]) ) ), 1. ) );
    }
  }
  return 1.0;
}
