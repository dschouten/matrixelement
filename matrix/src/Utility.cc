//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#define EXP_A 1512775
#define EXP_C 60801
#define EXP_D 1072693248

#include "matrix/Utility.hh"
#include "matrix/Fortran.hh"

#include <stdexcept>

namespace hepstd {
  const unsigned int iE  = 0;
  const unsigned int iPX = 1;
  const unsigned int iPY = 2;
  const unsigned int iPZ = 3;

  // ------------------------- ======= ------------------------- ======= -------------------------
  double phiMPiPi( double x )
  {
    if( !std::isnan( x ) )
    {
      while(x >  M_PI) x -= M_2PI;
      while(x < -M_PI) x += M_2PI;
    }
    return x;
  }
  
  // ------------------------- ======= ------------------------- ======= -------------------------
  double topDecayWidth( double mass )
  {
    const double factor = 1.7e-05;
    const double wMassWidth = hepstd::wMass * hepstd::wWidth;
    
    double ffact = factor * mass * mass * mass;
    double alpha = (hepstd::wMass / mass) * (hepstd::wMass / mass);
    double massFactor = (1 - 3 * alpha * alpha + 2 * alpha * alpha * alpha) / wMassWidth;
    double piFactor = atan( ( (mass)*(mass) - (hepstd::wMass*hepstd::wMass) ) / wMassWidth )
      - atan( -hepstd::wMass / hepstd::wWidth );
    
    return ( ffact * piFactor * massFactor );
  }
  
  // ------------------------- ======= ------------------------- ======= -------------------------
  double higgsDecayWidth( double mass )
  {
    const double factor = 4.8e-2;
    return ( factor * mass / (4.0 * M_SQRT_2 * M_PI) );
  }
  
  // ------------------------- ======= ------------------------- ======= -------------------------
  double dP( const TLorentzVector& vec )
  {
    return vec.Px() * vec.Py() * vec.Pz();
  }
  
  // ------------------------- ======= ------------------------- ======= -------------------------
  void tlvCopy( TLorentzVector& l, const TLorentzVector& r )
  {
    l.SetPx( r.Px() );
    l.SetPy( r.Py() );
    l.SetPz( r.Pz() );
    l.SetE( r.E() );
  }
  
  // ------------------------- ======= ------------------------- ======= -------------------------
  void tlvToArray( const TLorentzVector& vec, double array[] )
  {
    array[iE]  = vec.E();
    array[iPX] = vec.Px();
    array[iPY] = vec.Py();
    array[iPZ] = vec.Pz();
  }
  
  // ------------------------- ======= ------------------------- ======= -------------------------
  void swap( TLorentzVector& l, TLorentzVector& r )
  {
    static TLorentzVector b;
    b.SetPxPyPzE( r.Px(), r.Py(), r.Pz(), r.E() );
    r.SetPxPyPzE( l.Px(), l.Py(), l.Pz(), l.E() );
    l.SetPxPyPzE( b.Px(), b.Py(), b.Pz(), b.E() );
  } 

  // ------------------------- ======= ------------------------- ======= -------------------------
  double fexp( double y ) 
  {  

    return std::exp( y );

    /*
     * Not sure how this approximation goes with -O3 CXX flag
     *

    // A Fast, Compact Approximation of the Exponential Function
    //   by  Nicol N. Schraudolph
    union
    {
      double d;
#ifdef LITTLE_ENDIAN
      struct { int j, i; } n;
#else
      struct { int i, j; } n;
#endif
    } _eco;
    
    _eco.n.i = (int)(EXP_A*(y)) + (EXP_D - EXP_C);
    _eco.n.j = 0;
    return _eco.d;

    */
  }

  int sign( double x )
  {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
  }
  
  // ------------------------- ======= ------------------------- ======= -------------------------
  double erf( double y )
  {
    static double A = 8. * (M_PI - 3.) / (3. * M_PI * (4. - M_PI));
    static double B = 4. / M_PI;
    
    int S = y >= 0 ? 1 : -1;
    
    // Winitzki, Sergei (6 February 2008).
    // "A handy approximation for the error function and its inverse" 
    // Retrieved 2011-10-03.
    
    // max |error| < 3.5e-4 for all y
    
    double ysqr = y*y;
    double arg = -ysqr * ((B + A*ysqr) / (1 + A*ysqr));
    return S * sqrt( 1.0 - exp( arg ) );
  }
  
  // ------------------------- ======= ------------------------- ======= -------------------------
  double norm( double array[], unsigned size )
  {
    double sum = 0;
    for( unsigned i = 0; i < size; ++i )
    {
      sum += array[i] * array[i];
    }
    if( sum >= 0 and sum == sum )
    {
      return sqrt( sum );
    }
    else
    {
      std::cout << "argument = " << sum << " terms = ";
      for( unsigned i = 0; i < size; ++i )
      {
	std::cout << array[i] << " ";
      }
      std::cout << std::endl;
      throw std::runtime_error( "sqrt() of negative number is not real" );
    }
  }
  
  // ------------------------- ======= ------------------------- ======= -------------------------
  double norm( double x, double y )
  {
    double array[2] = { x, y };
    return norm( array, 2 );
  }
  
  // ------------------------- ======= ------------------------- ======= -------------------------
  double norm( double x, double y, double z )
  {
    double array[3] = { x, y, z };
    return norm( array, 3 );
  }
  
  // ------------------------- ======= ------------------------- ======= -------------------------
  double norm( double x, double y, double z, double t )
  {
    double array[4] = { x, y, z, t };
    return norm( array, 4 );
  }

  // ------------------------- ======= ------------------------- ======= -------------------------
  double safesqrt( double arg )
  {
    if( arg < -TINY )
    {
      throw negativeroot_exception( );
    }
    return sqrt( fabs(arg) );
  }

  // ------------------------- ======= ------------------------- ======= -------------------------
  double safesqrt( long double arg )
  {
    if( arg < -TINY)
    {
      throw negativeroot_exception( );
    }
    return sqrt( fabs(arg) );
  }
    
}
