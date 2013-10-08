//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/Complex.hh"

Complex operator*( double sc, const Complex& cl )
{
  return Complex( sc * cl.real, sc * cl.imag );
}

Complex operator/( double sc, const Complex& cl )
{
  return Complex( cl.real / sc, cl.imag / sc );
}

Complex operator+( double sc, const Complex& cl )
{
  return Complex( cl.real + sc, cl.imag );
}

Complex operator-( double sc, const Complex& cl )
{
  return Complex( cl.real - sc, cl.imag );
}

std::ostream& operator<<( std::ostream& out, const Complex& c )
{
  out << "(" << c.real << "," << c.imag << ")";
  return out;
}
