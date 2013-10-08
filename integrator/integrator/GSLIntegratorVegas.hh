//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#ifndef GSLINTEGRATORVEGAS_HH
#define GSLINTEGRATORVEGAS_HH

#include "integrator/Integrator.hh"

class GSLIntegratorVegas : public Integrator
{
public:
  GSLIntegratorVegas();
  virtual ~GSLIntegratorVegas() {}
  
  virtual void doIntegral(double returnVal[], double error[], int* fail,
			  int* neval, double prob[]) const;

  void setNSteps( unsigned nsteps ) { m_nsteps = nsteps; }
  void setMaxChiSquare( float chisqr ) { m_max_chisqr = chisqr; }
  void setAllowNull( bool flag ) { m_allow_null = flag; }
private:
  unsigned m_nsteps;
  float m_max_chisqr;
  bool m_allow_null;
};

#endif
