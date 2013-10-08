//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#ifndef GSLINTEGRATORMISER_HH
#define GSLINTEGRATORMISER_HH

#include "integrator/Integrator.hh"

class GSLIntegratorMiser : public Integrator
{
public:
  GSLIntegratorMiser();
  virtual ~GSLIntegratorMiser() {}
  
  virtual void doIntegral(double returnVal[], double error[], int* fail,
			  int* neval, double prob[]) const;
};

#endif
