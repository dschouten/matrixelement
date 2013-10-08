//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/JetEnergyResolutionTF.hh"

#include <TF1.h>

#include <iostream>
#include <string>
#include <math.h>
#include <stdexcept>

#define NPAR 10

#define VERSION_1 1
#define VERSION_2 2

namespace
{
  class Functor
  {
  public:
    Functor( JetEnergyResolutionTF* tf_ptr ) : m_tf_ptr( tf_ptr ) { }
    
    double operator()( double* x, double* par )
    {
      if( m_tf_ptr == NULL )
      {
	return 0.;
      }
      return (*m_tf_ptr)( "e", (*x), par );
    }
    
    JetEnergyResolutionTF* m_tf_ptr;
  };  
  
  class FunctionWrapper : public TF1
  {
  public:
    FunctionWrapper( Functor* f, double l, double h ) : TF1( "f", f, l, h, 1 ) { }
    double operator()( double x ) { return Eval( x ); }
  };
}

// ------------------------- ======= ------------------------- ======= -------------------------
JetEnergyResolutionTF::JetEnergyResolutionTF( const std::string& defn ) : 
  TransferFunction( defn ),
  m_sigma_e( -1 ),
  m_par( NPAR, -1 ),
  m_lower_bound( 20.0 ),
  m_upper_bound( 500.0 ),
  m_parameters_override( false ),
  m_norm_flag( true )
{
  m_version = VERSION_1;
  
  bool flag = false;
  m_sigma_e = parameter( "sigma", flag );
  
  if( !flag )
  {
    m_version = VERSION_2;  
    char par[16];  
    for( unsigned i = 0; i < NPAR; ++i )
    {
      sprintf( par, "par_%d", i );
      m_par[i] = parameter( par, flag );
      if( !flag )
      {
	std::cout << "ERROR could not read parameter [" << par << "] for JetEnergyResolutionTF!" << std::endl;
	throw std::runtime_error("could not read parameter");
      }
    }
    m_lower_bound = parameter( "lower_bound", flag, 20.0  );
    m_upper_bound = parameter( "upper_bound", flag, 500.0 );
  }
  
  m_pt_cutoff = parameter( "pt_cutoff",            flag, -1   );
  m_ps_frac   = parameter( "phase_space_fraction", flag, 0.90 );
  
  SHOW_DEBUG( std::cout << "Initialized JetEnergyResolution TF v" << m_version << std::endl );
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool JetEnergyResolutionTF::limits( const char* name, const double& val, const double par[], double vlim[] ) 
{
  if( m_version == VERSION_1 )
  {
    if( name[0] == 'p' || name[0] == 'e' )
    {
      double sigma = val * (m_sigma_e / sqrt( val ));
      vlim[0] = val - 4 * sigma;
      vlim[1] = val + 4 * sigma;
      return true;
    }
  }
  
  if( m_version == VERSION_2 )
  {
    if( name[0] == 'p' || name[0] == 'e' )
    {
      bool lflag = false;
      bool hflag = false;
      
      double dE  = 0.0;
      double ETA = m_parameters_override ? m_par_eta : par[0];
      double E   = m_parameters_override ? m_par_emeas : val;

      const double delta = 0.5;
      
      vlim[0] = m_pt_cutoff > 0 ? std::max( m_lower_bound, m_pt_cutoff * cosh(ETA) ) : m_lower_bound;
      vlim[1] = m_upper_bound;
            
      m_norm_flag = false;

      Functor* fobj         = new Functor( this );
      FunctionWrapper* fptr = new FunctionWrapper( fobj, m_lower_bound, m_upper_bound );

      fptr->SetParameter( 0, E );  
      fptr->SetParameter( 1, ETA );

      double peak  = fptr->GetMaximum(); 
      double start = fptr->GetMaximumX();
      
      // std::cout << "peak X: " << start << std::endl;
      // std::cout << "peak Y: " << peak << std::endl;
      
      // now iterate in steps of size delta to find approximate bounds
      // such that range contains ~ 95% 

      while( dE < (m_upper_bound - m_lower_bound) )
      {	
	if( !hflag && ( (start + dE) > m_upper_bound || 
			fptr->Integral( start + dE, m_upper_bound ) / fptr->Integral( m_lower_bound, m_upper_bound ) < 0.01 ) )
	{
	  hflag = true;
	  vlim[1] = std::min( m_upper_bound, start - dE ); 
	}

	if( !lflag && ( (start - dE) < m_lower_bound || 
			fptr->Integral( m_lower_bound, start - dE ) / fptr->Integral( m_lower_bound, m_upper_bound ) < 0.01 ) )
	{
	  lflag = true;
	  vlim[0] = std::max( m_lower_bound, start - dE );
	}
	
	// std::cout << std::max( m_lower_bound, start - dE ) << " " 
	// 	  << std::min( m_upper_bound, start + dE ) << " " 
	// 	  << fptr->Integral( std::max( m_lower_bound, start - dE ), 
	// 			     std::min( m_upper_bound, start + dE ) ) << std::endl;
	
	if( lflag && hflag )
	  break;

	if( fptr->Integral( std::max( m_lower_bound, start - dE ), 
			    std::min( m_upper_bound, start + dE ) ) / fptr->Integral( m_lower_bound, m_upper_bound ) > m_ps_frac )
	  break;
	
	dE += delta;
      }

      m_norm_flag = true;

      if( !lflag )
      {
	vlim[0] = std::max( m_lower_bound, start - dE );
      }
      if( !hflag )
      {
	vlim[1] = std::min( m_upper_bound, start + dE );
      }
      
      delete fobj;
      delete fptr;

      return true;
    }
  }
  return false;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double JetEnergyResolutionTF::operator()( const char* name, const double& x, const double par[] ) 
{
  double emeas   = m_parameters_override ? m_par_emeas : par[0];
  double eta     = m_parameters_override ? m_par_eta   : par[1];
  double eparton = x;
  
  if( m_pt_cutoff > 0 && eparton / cosh(eta) < m_pt_cutoff )
  {
    return 0;
  }

  if( eparton < m_lower_bound || eparton > m_upper_bound )
  {
    return 0;
  }

  double norm = m_norm_flag ? integral( emeas, eta ) : 1.0;

  if( m_version == VERSION_1 )
  {
    if( name[0] == 'p' || name[0] == 'e' )
    {
      double sigma = emeas * (m_sigma_e / sqrt( emeas ));
      double arg = -((eparton - emeas)*(eparton - emeas)) / (2.0 * sigma * sigma);
      return ( 1.0 / ( M_SQRT_2PI * sigma ) * hepstd::fexp( arg ) ) / norm; // probability for parton E given measured E
    }
  }
  if( m_version == VERSION_2 )
  {
    if( name[0] == 'p' || name[0] == 'e' )
    {
      double difference = emeas - eparton;
      
      double p01 = m_par[0] + m_par[1] * eparton;
      double p23 = m_par[2] + m_par[3] * eparton;
      double p45 = m_par[4] + m_par[5] * eparton;
      double p67 = m_par[6] + m_par[7] * eparton;
      double p89 = m_par[8] + m_par[9] * eparton;
      
      double fxy = hepstd::fexp(-.5 * pow((difference - p01) / p23,2));
      fxy += p45 * hepstd::fexp(-.5 * pow((difference - p67) / p89,2));
      fxy /= p23 + p45 * p89;
      fxy /= M_SQRT_2PI; 
      return fxy / norm;
    }
  }
  return 1.0;
}

// ------------------------- ======= ------------------------- ======= -------------------------
inline double JetEnergyResolutionTF::integral( double emeas, double eta )
{
  static double emeas_cache = 0;
  static double eta_cache   = 0;
  static double integ_cache = 1;

  if( fabs( emeas_cache - emeas ) < 1.0e-8 && fabs( eta_cache - eta ) < 1.0e-8 )
  {
    return integ_cache;
  }

  Functor fobj( this );
  FunctionWrapper fptr( &fobj, m_lower_bound, m_upper_bound );

  fptr.SetParameter( 0, emeas );
  fptr.SetParameter( 1, eta );

  m_norm_flag = false;

  emeas_cache = emeas;
  eta_cache   = eta;
  integ_cache = fptr.Integral( m_lower_bound, m_upper_bound );

  m_norm_flag = true;

  return integ_cache;
}

// FIXME THIS IS BULLSHIT ...
