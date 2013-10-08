//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#ifndef GSLINTEGRATOR1D_HH
#define GSLINTEGRATOR1D_HH

#include "integrator/Integrator.hh"

#include <cmath>
#include <cstdlib>

class GSLIntegrator1D : public Integrator
{
public:
  GSLIntegrator1D();
  virtual ~GSLIntegrator1D() {}
  
  virtual void doIntegral(double returnVal[], double error[], int* fail,
			  int* neval, double prob[]) const;
  
private:
  double myQuad( integrand_t f_ptr, double& errest, double a, double b, int n ) const;
  
};


/* ==============================================================
    Numerical integration of f(x) on [a,b]
    method: Monte-Carlo method
    written by: Alex Godunov (February 2007)
----------------------------------------------------------------
input:
    f   - a single argument real function (supplied by the user)
    a,b - the two end-points of the interval of integration
    n   - number of random points for xi
output:
    r      - result of integration
    errest - estimated error
Comments:
    be sure that following headers are included
    #include <cstdlib>
    #include <ctime>     
================================================================ */

double GSLIntegrator1D::myQuad( integrand_t f_ptr, double& errest, double a, double b, int n ) const
{
  using namespace std;
  
  double r, x, u, gs, fs;
  /* variables fs and gs are used to estimate an error of integration */    
  
  srand(time(NULL));               /* initial seed value (use system time) */
  
  fs = 0.0;
  gs = 0.0;

  double buff;

  for (int i = 1; i <= n; i=i+1)
  {
    u = 1.0*rand()/(RAND_MAX+1); /* random number between 0.0 and 1.0 */
    x = a + (b-a)*u;
    (f_ptr)( NULL, &x, NULL, &buff, NULL );
    fs = fs + buff;  // fs + f(x);
    gs = gs + buff*buff; // gs + f(x)*f(x);
  }
  r = fs*(b-a)/n;

  /* evaluate integration error */
  fs = fs/n;
  gs = gs/n;
  errest = (b-a)*sqrt((gs - fs*fs)/n);
  return r;
}

#endif
