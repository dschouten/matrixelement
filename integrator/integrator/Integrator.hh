//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#ifndef INTEGRATOR_HH
#define INTEGRATOR_HH

#include <string>
#include <vector>
#include <iostream>


#ifndef __CINT__
namespace cuba {
#include "cuba/cuba.h"
}

namespace cubature {
#include "cubature/cubature.h"
}

#endif

typedef int (*integrand_t)(const int *ndim, const double x[],
			   const int *ncomp, double f[], void *userdata);

//////////////////////////////////////////////////
//
// Base integrator class
//
// Styled after the functionality in Cuba library
//
//////////////////////////////////////////////////

class Integrator
{
public:
  Integrator( std::string name );
  virtual ~Integrator() {}
  
  // perform the integration
  virtual void doIntegral( double returnVal[], double error[], int* fail,
			   int* neval, double prob[] ) const = 0;

  ///////////////////
  // setters
  ///////////////////

  virtual void setNDimensions( int input )  { m_ndim = input; }
  virtual void setNComp( int input )  { m_ncomp = input; }
  virtual void setIntegrand( integrand_t input )  { m_integrand = input; }
  virtual void setUserdata( void *userdata ) { m_userdata = userdata; }
  virtual void setEpsilon( double relative, double absolute ) { m_epsrel = relative; m_epsabs = absolute; }
  virtual void setNEval( int min, int max )  { m_mineval = min; m_maxeval = max; }
  virtual void setRandomSeed( int seed )  { m_ranseed = seed; }

  inline void setSampleSet( bool input );
  inline void setVerbose( int input );
  inline void setPseudoRandom( bool input );

  inline void setTimeout( unsigned timeout );

  ///////////////////
  // getters
  ///////////////////

  std::string getName()      const { return m_name; }

  int getNDimensions()       const { return m_ndim; }
  int getNComp()             const { return m_ncomp; }
  int getRandomSeed()        const { return m_ranseed; }
  int getMinEval()           const { return m_mineval; }
  int getMaxEval()           const { return m_maxeval; }
  int getFlags()             const { return m_flags; }

  double getEpsilonRel()     const { return m_epsrel; }
  double getEpsilonAbs()     const { return m_epsabs; }

  unsigned getTimeout( )     const { return m_timeout; }

  integrand_t getIntegrand() const { return m_integrand; }

protected:
  integrand_t m_integrand; //< pointer to function that calculates integrand at a given point in phase space
  void *m_userdata;
  std::string m_name; //< the name of the integrator
  int m_ndim; //< number of integration dimensions 
  int m_ncomp; //< number of components of the integrand (scalar = 1)
  double m_epsrel; //< relative precision desired
  double m_epsabs; //< absolute precision
  int m_flags; //< see note below @FLAGS
  int m_mineval; //< the minimum number of integrand evaluations required
  int m_maxeval; //< ... and the maximum
  int m_ranseed; //< random seed
  unsigned m_timeout; //< timeout in seconds
};

/*
 * FLAGS governing the integration:
 *
 *
 * - Bits 0 and 1 encode the verbosity level, i.e. 0 to 3. Level 0 does not print any 
 *   output, level 1 prints ‘reasonable’ information on the progress of the integration, 
 *   level 2 also echoes the input parameters, and level 3 further prints the subregion 
 *   results (if applicable).
 *
 * - Bit 2 = 0, all sets of samples collected on a subregion during the various itera-
 *   tions or phases contribute to the final result.
 * - Bit 2 = 1, only the last (largest) set of samples is used in the final result.
 *
 * - Bit 3 = 0, Sobol quasi-random numbers are used for sampling,
 * - Bit 3 = 1, Mersenne Twister pseudo-random numbers are used for sampling.
 *
 */

void Integrator::setVerbose( int input ) 
{ // 0x0011 mask
  m_flags &= 8;
  m_flags |= input & 3;
}

void Integrator::setSampleSet( bool input )
{ // 0x0100 mask
  m_flags &= 11;
  m_flags |= ( input << 2 ) & 4;
}

void Integrator::setPseudoRandom( bool input )  
{ // 0x1000 mask
  m_flags &= 7;
  m_flags |= ( input << 3 ) & 8;
}

void Integrator::setTimeout( unsigned timeout )
{
  m_timeout = timeout;
}

//////////////////////////////////
//
// Dummy integrator class
//
//////////////////////////////////

class NullIntegrator : public Integrator
{ 
public:
  NullIntegrator();

  virtual ~NullIntegrator() {  }

  virtual void doIntegral( double returnVal[], double error[], int* fail,
			  int* neval, double prob[] ) const;

  void setParam( unsigned param, double value );

private:
  std::vector<double> m_params;
};

#endif
