
// author: Doug Schouten <doug dot schouten at triumf dot ca>

#ifndef COMPLEX_HH
#define COMPLEX_HH

#include <iostream>

typedef struct{ double real, imag; } dbl_cpx;

class Complex
{
public:
  Complex( double real, double imag )
  {
    this->real = real;
    this->imag = imag;
  }

  Complex( double scalar )
  {
    this->real = scalar;
    this->imag = 0.;
  }

  Complex( )
  {
    this->real = 0;
    this->imag = 0;
  }

  Complex( const Complex& cl )
  {
    this->real = cl.real;
    this->imag = cl.imag;
  }

  ~Complex() { }
  
  Complex& operator=( const Complex& cl )
  {
    this->real = cl.real;
    this->imag = cl.imag;
    return (*this);
  }

  // <<<<<<<<<<< accessor as anonymous struct >>>>>>>>>>>>>>
  
  dbl_cpx operator()( ) 
  {
    dbl_cpx buffer = { real, imag };
    return buffer;
  }

  // <<<<<<<<<<< arithmetic operations with complex #'s >>>>>>>>>>>>>>

  Complex operator+( const Complex& cl ) const
  {
    return Complex( real + cl.real, imag + cl.imag );
  }

  Complex operator-( const Complex& cl ) const
  {
    return Complex( real - cl.real, imag - cl.imag );
  }

  Complex operator-() const
  {
    return Complex( -real, -imag );
  }

  Complex operator*( const Complex& cl ) const
  {
    return Complex( real*cl.real - imag*cl.imag, 
		    real*cl.imag + imag*cl.real );
  }
  
  Complex operator/( const Complex& cl ) const
  {
    double norm = cl.real*cl.real + cl.imag*cl.imag;
    return Complex( (real*cl.real + imag*cl.imag)/norm, 
		    (imag*cl.real - real*cl.imag)/norm );
  }

  Complex& operator+=( const Complex& cl )
  {
    real += cl.real;
    imag += cl.imag;
    return (*this);
  }

  Complex& operator-=( const Complex& cl )
  {
    real -= cl.real;
    imag -= cl.imag;
    return (*this);
  }

  Complex& operator*=( const Complex& cl )
  {
    Complex buffer = (*this) * cl;
    (*this) = buffer;
    return (*this);
  }

  Complex& operator/=( const Complex& cl )
  {
    Complex buffer = (*this) / cl;
    (*this) = buffer;
    return (*this);
  }

  // <<<<<<<<<<< operations with real #'s >>>>>>>>>>>>>>

  Complex operator+( double cl ) const
  {
    return Complex( real + cl, imag );
  }

  Complex operator-( double cl ) const
  {
    return Complex( real - cl, imag );
  }

  Complex operator*( double cl ) const
  {
    return Complex( real*cl, imag*cl );
  }
  
  Complex operator/( double cl ) const
  {
    return Complex( real/cl, imag/cl );
  }

  Complex& operator+=( double cl )
  {
    real += cl;
    return (*this);
  }

  Complex& operator-=( double cl )
  {
    real -= cl;
    return (*this);
  }

  Complex& operator*=( double cl )
  {
    real *= cl;
    imag *= cl;
    return (*this);
  }

  Complex& operator/=( double cl )
  {
    real /= cl;
    imag /= cl;
    return (*this);
  }

  //  <<<<<<<<<<< member data >>>>>>>>>>>>>>

  double real;
  double imag;
};

Complex operator*( double sc, const Complex& cl );
Complex operator/( double sc, const Complex& cl );
Complex operator+( double sc, const Complex& cl );
Complex operator-( double sc, const Complex& cl );

std::ostream& operator<<( std::ostream& out, const Complex& c );

#endif
