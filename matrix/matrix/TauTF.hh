#ifndef TAU_TF
#define TAU_TF

//Standard includes
#include <string>

//Matrix includes
#include "matrix/TransferFunction.hh"

class TauTF : public TransferFunction
{
public:
  
  /* Construction/Destruction */
  
  TauTF();
  virtual ~TauTF();
  
  
  /* Public Limit/Evaluation Methods */
  
  virtual bool limits(const char * name, const double & val, const double par[], double vlim[]);

  virtual double operator()( const char * name, const double & eratio, /* hypothetical ratio of tauon to lepton energy */  
			     const double parameters[] /* measured lepton energy, charge */ );
  
private:
  
  /* Private (Disabled) Constructors */
  
  TauTF(const std::string & defn);
  
  
  /* Actual TF implementation */
  
  double _fragmentation_probability(double hypothesized_tau_energy, double measured_lepton_energy, bool antipaticle);
  
  
  /* Private Variables */
  double _tau_mass; //GeV
  double _michel_param_rho;
  double _michel_param_eta;
  double _michel_param_xi;
  double _michel_param_delta;
};

#endif
