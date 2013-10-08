//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef JETENERGYRESOLUTION_TF
#define JETENERGYRESOLUTION_TF

#include "matrix/TransferFunction.hh"

//////////////////////////////////////////////////////////
//
// Transfer function for jets 
//
// - provides probability distribution for parton energy
//   given observed energy
//
// Version 1.0
// - energy resolution provided for single pseudorapidity bin 
//
// Version 2.0
// - using TF from Bernd, double Gaussian fit
// 
//////////////////////////////////////////////////////////

class JetEnergyResolutionTF : public TransferFunction
{
public:
  JetEnergyResolutionTF( const std::string& defn );
  virtual ~JetEnergyResolutionTF() { }

  // return limits [a,b] such that 1 - int[a,b] f(E) dE <= 0.01 
  // parameters
  //    char : 'E' or 'p' -> E smearing
  //    double : E -> E of the jet
  //    double [] : parameters
  //    double [] : limits [a,b]
  virtual bool limits( const char*, const double&, const double[], double[] );
  virtual double operator()( const char*, const double&, const double[] );

  void setParametersOverride( bool, double, double, double );

private:
  JetEnergyResolutionTF( ) { }
  
  double m_sigma_e; //<< for version 1.0 TF specifications, a simple energy dependent Gaussian smearing
  
  std::vector<double> m_par; //<< for version 2.0 TF, an energy dependent double-Gaussian 
  double m_lower_bound;
  double m_upper_bound;
  double m_pt_cutoff;
  double m_ps_frac; //<< fraction of TF phase space to integrate over
  
  bool m_parameters_override; //<< set to override the parameters used to define the TF
  double m_par_emeas;
  double m_par_eta;

  unsigned m_version; //<< version of the TF defined

  bool m_norm_flag; //<< flag to ignore normalization

  double integral( double emeas, double eta );
};

#endif
