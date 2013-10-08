#ifndef CONSTANTS_HH
#define CONSTANTS_HH

#include <cmath>
#include <limits>

#include "matrix/GlobalFlags.hh"

#define DEBUG // allows to print debug statements

// #define BENCHMARK // benchmark components of the integrand

#define MAX_DIMENSIONS 10 // maximum number of integration dimensions

#define USECTEQ // use CTEQ pdf routines directly (instead of LHAPDF) << faster

#define NCOMBSTORED 512 // maximum # of helicity combinations to keep track of

#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif 

#ifndef M_PI_2 
#define M_PI_2 9.869604401089358
#endif 

#ifndef M_2PI 
#define M_2PI 6.283185307179586
#endif 

#ifndef M_SQRT_2
#define M_SQRT_2 1.4142135623730951
#endif

#ifndef M_SQRT_2PI
#define M_SQRT_2PI 2.5066282746310002
#endif

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

#ifndef __CINT__ // ugh, CINT ... WTF?
#define EPSILON std::numeric_limits<double>::epsilon()
#else
#define EPSILON 1.0e-16
#endif

#define TINY 1.0e-6

#include "matrix/GlobalFlags.hh"

namespace hepstd
{
  enum PartonType
    {
      kunknown = 100,
      kgluon = 0,
      kd = 1,
      ku = 2,
      ks = 3,
      kc = 4,
      kb = 5,
      kdbar = -1,
      kubar = -2,
      ksbar = -3,
      kcbar = -4,
      kbbar = -5
    };

  const int ktau      = 15;
  const int kmuon     = 13;
  const int kelectron = 11;

  // all units are in GeV ...
  
  const double beamEnergy = 4000.0;

  //
  //  << MG 5 values 
  //

  const double tauMass = 1.77682;
  const double bMass   = 4.7000; 
  const double cMass   = 1.4000; 
  const double tMass   = 173.00;
  const double wMass   = 80.4190024;
  const double wWidth  = 2.04760000;
  const double zMass   = 91.1876000;
  const double zWidth  = 2.44140000;
}

#endif
