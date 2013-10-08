//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#ifndef ROOTINTEGRATOR_HH
#define ROOTINTEGRATOR_HH

#include "integrator/Integrator.hh"

class RootIntegrator : public Integrator
{
public:
  RootIntegrator();
  virtual ~RootIntegrator() {}
  
  virtual void doIntegral(double returnVal[], double error[], int* fail,
			  int* neval, double prob[]) const;
  
private:
  
};

#endif
