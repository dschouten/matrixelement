//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "integrator/GSLIntegratorVegas.hh"

#include <iostream>
#include <stdexcept>
#include <vector>

#include <cmath>

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

using std::vector;

namespace 
{
  integrand_t f_ptr; //< function pointer for integrand
  double ret_val; //< store y for y = f(x)
  long ncalls; //< number of function evaluations recorded
  
  double integrand( double x[], size_t dim, void* params )
  {
    ncalls++;
    (f_ptr)( NULL, x, NULL, &ret_val, params);
    return ret_val;
  }
};

// ------------------------- ======= ------------------------- ======= -------------------------
GSLIntegratorVegas::GSLIntegratorVegas() :
  Integrator("GSL VEGAS"), m_nsteps( 10 ), m_max_chisqr( 1. ), m_allow_null( true )
{
  
}

// ------------------------- ======= ------------------------- ======= -------------------------
void GSLIntegratorVegas::doIntegral( double returnval[], double relerr[],
				     int* ifail, int* nfnevl, double prob[] ) const
{
  ncalls = 0;
  
  double params = 0;
  double result, error;
  
  unsigned ndim = getNDimensions();
  
  f_ptr = getIntegrand();
  
  gsl_monte_function F;
  
  F.f = integrand;
  F.dim = ndim;
  F.params = &params;
  
  double* xl = new double[ndim];
  double* xh = new double[ndim];
  for( unsigned i = 0; i < ndim; ++i )
  {
    xl[i]=0;
    xh[i]=1;
  }
  
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc( ndim );
  
  gsl_rng_env_setup();
  
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *r = gsl_rng_alloc (T);
  
  unsigned npts_startup = static_cast<unsigned>( getMinEval() );
  unsigned npts = getMaxEval() - getMinEval();
  
  gsl_monte_vegas_integrate(&F, xl, xh, ndim, npts_startup, r, s, &result, &error);
  
  prob[0]	= 0;
  (*ifail)	= 0;      
  
  unsigned ipart = 0;
  do
  {
    gsl_monte_vegas_integrate(&F, xl, xh, ndim, static_cast<unsigned>( ((double)npts) / m_nsteps ), r, s, &result, &error);
    ipart += 1;
    if( ipart > m_nsteps + 1 )
    {
      (*ifail) = 1;
      break;
    }
    prob[0] = 0.0; // gsl_monte_vegas_chisq(s); /* GSL 1.10 doesn't have this function defined - dschoute */
    std::cout << "\tVEGAS iteration #" << ipart 
	      << " (" << static_cast<unsigned>( ((double)npts) / m_nsteps ) 
	      << "," << prob[0]
	      << "," << (result != 0 ? fabs( error / result ) : 1.)
	      << ")" << std::endl;
  }
  while( fabs(prob[0]) > m_max_chisqr || 
	 ( result != 0 && fabs( error / result ) > getEpsilonRel() ) || ( result == 0 and !m_allow_null ) );
  
  returnval[0]	= result;
  relerr[0]	= error;
  (*nfnevl)	= (ipart) * static_cast<unsigned>( ((double)npts) / m_nsteps ) + npts_startup; 
  
  gsl_rng_free (r);         
  gsl_monte_vegas_free (s);
}
