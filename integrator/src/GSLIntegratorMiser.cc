//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "integrator/GSLIntegratorMiser.hh"

#include <iostream>
#include <stdexcept>
#include <vector>

#include <cmath>

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>

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

  const unsigned NSTEPS = 5;
};

// ------------------------- ======= ------------------------- ======= -------------------------
GSLIntegratorMiser::GSLIntegratorMiser() :
  Integrator("GSL MISER")
{

}

// ------------------------- ======= ------------------------- ======= -------------------------
void GSLIntegratorMiser::doIntegral( double returnval[], double relerr[],
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

  gsl_monte_miser_state* s = gsl_monte_miser_alloc( ndim );
  gsl_monte_miser_init( s );      

  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *r = gsl_rng_alloc (T);

  unsigned npts = getMinEval();
  unsigned npts_delta = getMaxEval() - getMinEval();

  prob[0]	= 0;
  (*ifail)	= 0;      

  while( npts < static_cast<unsigned>( getMaxEval() ) )
  {
    gsl_monte_miser_integrate( &F, xl, xh, ndim, npts, r, s, &result, &error ); 
    if( result != 0 && fabs( error / result ) > getEpsilonRel() ) 
    {
      npts += static_cast<unsigned>( npts_delta / ((float)NSTEPS) );
      continue;
    }
    else
    {
      break;
    }
  }

  if( result != 0 && fabs( error / result ) > getEpsilonRel() )
  {
    (*ifail) = 1;
  }
  
  returnval[0]	= result;
  relerr[0]	= error;
  (*nfnevl)	= npts; 

  gsl_rng_free (r);         
  gsl_monte_miser_free (s);
}
