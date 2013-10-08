//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "integrator/GSLIntegrator1D.hh"

#include <iostream>
#include <stdexcept>
#include <vector>

#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

using std::vector;

namespace 
{
  integrand_t f_ptr; //< function pointer for integrand
  double ret_val; //< store y for y = f(x)
  long ncalls; //< number of function evaluations recorded

  double integrand( double x, void* params )
  {
    ncalls++;
    (f_ptr)( NULL, &x, NULL, &ret_val, params);
    return ret_val;
  }

  void err_handler( const char * reason, 
		    const char * file, 
		    int line, 
		    int gsl_errno )
  {
    throw std::runtime_error( reason );
  }

};

// ------------------------- ======= ------------------------- ======= -------------------------
GSLIntegrator1D::GSLIntegrator1D() :
  Integrator("GSL 1D")
{

}

// ------------------------- ======= ------------------------- ======= -------------------------
void GSLIntegrator1D::doIntegral( double returnval[], double err[],
				   int* ifail, int* nfnevl, double prob[] ) const
{
  ncalls = 0;

  double params = 0;
  double result(1), error(1);

  f_ptr = getIntegrand();

  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc( getMaxEval() );
  
  gsl_function F;

  F.function = integrand;
  F.params   = &params;

  gsl_error_handler_t* handler = gsl_set_error_handler( &err_handler );

  try
  {
    gsl_integration_qags (&F, 0, 1, 0,
    			  getEpsilonRel(), getMaxEval(), w, &result, &error); 
    
    // long neval = getMinEval();

    // while( ncalls == 0 || fabs(error / result) > getEpsilonRel() )
    // {
    //   result  = myQuad( getIntegrand(), error, 0, 1, neval );
    //   ncalls += neval;
    //   if( ncalls >= getMaxEval() )
    // 	break;
    //   neval  += (getMaxEval() - getMinEval()) / 3;
    // }
  }
  catch( std::runtime_error err )
  {
    std::cout << "WARNING std::runtime exception in GSLIntegrator1D: " << err.what() << std::endl;
    result = -2;
    (*ifail) = 1;
    gsl_integration_workspace_free (w);
    return;
  }

  gsl_integration_workspace_free (w);

  gsl_set_error_handler( handler );
       
  returnval[0]	= result;
  err[0]	= error;
  (*nfnevl)	= ncalls;
  prob[0]	= 0;
  (*ifail)	= !( fabs(error / result) <= getEpsilonRel() );

}
