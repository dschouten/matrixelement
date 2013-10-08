//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef TRANSFERFUNCTION_HH
#define TRANSFERFUNCTION_HH

#include <exception>
#include <stdexcept>
#include <map>
#include <cstring>
#include <string>

#include "matrix/Constants.hh"

// ------------------------- ======= ------------------------- ======= -------------------------
class TransferFunctionReadException : public std::exception
{
public:
  TransferFunctionReadException() : std::exception() { }

  const char* what() {
    return "FAIL to initialize transfer function"; 
  }
};

// ------------------------- ======= ------------------------- ======= -------------------------
class TransferFunction {
public:
  TransferFunction( );
  TransferFunction( const std::string& defn );
  virtual ~TransferFunction( ) { }

  // return a named parameter
  double parameter( const char* );
  double parameter( const char*, bool& );
  double parameter( const char*, bool&, double );

  std::string str_parameter( const char* );
  std::string str_parameter( const char*, bool& );

  // return limits [a,b] over which the transfer function is non-zero
  virtual bool limits( const char*, const double&, const double[], double[] ) = 0;

  // return probability for observed value 
  // @param par[] - array of parameters for transfer function
  // @param x - the value to be transformed
  virtual double operator()( const char*, const double&, const double[] ) = 0;

  long ncalls( ) const;
  
protected:
  std::map<std::string, double> m_parameters;
  std::map<std::string, std::string> m_parameters_raw;

private:
  mutable long m_ndebug;
};

inline long TransferFunction::ncalls( ) const 
{
  m_ndebug = m_ndebug <= GlobalFlags::max_debug ? m_ndebug+1 : m_ndebug; 
  return m_ndebug; 
}

#endif
