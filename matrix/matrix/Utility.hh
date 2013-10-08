//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef UTILITY_HH
#define UTILITY_HH

#include <iostream>
#include <iomanip>
#include <exception>
#include <stdexcept>

#include <TLorentzVector.h>

#include "matrix/Constants.hh"

namespace hepstd {
  
  class negativeroot_exception : public std::exception
  {
  public:
    negativeroot_exception( ) : std::exception() { }

  };

  double topDecayWidth( double );

  double higgsDecayWidth( double );

  double fexp( double );

  int sign( double );

  double erf( double );

  double dP( const TLorentzVector& );

  void tlvCopy( TLorentzVector&, const TLorentzVector& );

  void tlvToArray( const TLorentzVector&, double[] );

  double norm( double[], unsigned );

  double norm( double, double );

  double norm( double, double, double );

  double norm( double, double, double, double );

  void swap( TLorentzVector&, TLorentzVector& );
  
  double phiMPiPi( double );

  double safesqrt( double );

  double safesqrt( long double );
}

#ifdef BENCHMARK
#define STRT_TIMER(t) t.Start( false );
#define STOP_TIMER(t) t.Stop( );
#else
#define STRT_TIMER(t) 
#define STOP_TIMER(t) 
#endif

#ifdef DEBUG
#define SHOW_DEBUG(expression) if( GlobalFlags::debug && ncalls() <= GlobalFlags::max_debug ) { expression; }
#else
#define SHOW_DEBUG(expression) 
#endif

namespace hepstd 
{
  template<int NDIM>
  void swaparry( double array[][NDIM], unsigned i, unsigned j ) 
  {
    double buff[NDIM];
    for( unsigned k = 0; k < NDIM; ++k )
    {
      buff[k] = array[i][k];
      array[i][k] = array[j][k];
      array[j][k] = buff[k];
    }
  }
}

#endif

