//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/TauTF.hh"

//Standard includes
#include <iostream>

//Standard namespaces
using namespace std;


/* Construction/Destruction */

TauTF::TauTF() : TransferFunction()
{
  _tau_mass = hepstd::tauMass;
  _michel_param_rho = 0.75;
  _michel_param_eta = 0.0;
  _michel_param_xi = 1.0;
  _michel_param_delta = 0.75;
}

TauTF::TauTF(const std::string & defn) : TransferFunction(defn)
{
  //Disabled private constructor
}

TauTF::~TauTF()
{
  
}


/* Public Limit/Evaluation Methods */

bool TauTF::limits(const char * name, const double & val, const double par[], double vlim[])
{
  vlim[0] = 0.5;
  vlim[1] = 1.0;
  return true;
}

double TauTF::operator()( const char * name, const double & eratio, /* hypothetical ratio of tauon to lepton energy */  
			  const double parameters[] /* measured lepton energy, charge */ )
{
  double emeas  = parameters[0];
  double charge = parameters[1];
  if( eratio < TINY )
  {
    return 0;
  }
  return _fragmentation_probability( emeas/eratio, emeas, (charge > 0) );
}


double TauTF::_fragmentation_probability(double hypothesized_tau_energy, double measured_lepton_energy, bool antiparticle)
{
  if(hypothesized_tau_energy < measured_lepton_energy
     || hypothesized_tau_energy > (2 * measured_lepton_energy))
  {
    //Not physically possible
    return 0.0;
  }
  
  //Calculate the energy fraction (x in the PDG)
  double x = measured_lepton_energy/hypothesized_tau_energy;
  
  //Calculate powers for efficiency
  double x2 = x*x;
  double x3 = x2*x;
  
  //Calculate functions
  double f0 = 2.0 - 6.0*x2 + 4.0*x3;
  double f1 = -4.0/9.0 + 4.0*x2 - 32.0/9.0*x3;
  double f2 = 12.0 * (1.0-x)*(1.0-x);
  
  double g1 = -2.0/3.0 +4.0*x -6.0*x2 +8.0*x3/3.0;
  double g2 =  4.0/9.0 -16.0*x/3.0 + 12.0*x2 - 64.0*x3/9.0;
  
  //Calculate the polarization term (normal polarization + small weak mixing)
  double anti_factor = antiparticle ? -1 : 1;
  double pol_tau = (-1 * anti_factor) + anti_factor * (_tau_mass / hypothesized_tau_energy);
  double polarization_term = g1 + _michel_param_delta * g2;
  polarization_term *= pol_tau * _michel_param_xi;
  
  //Calculate the final result
  return f0 + _michel_param_rho * f1 + _michel_param_eta * 0.055 * f2 - polarization_term;
}
