//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/TransferFunction.hh"

#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
#include <sstream>
#include <string>

// using namespace std;

namespace {
  // ------------------------- ======= ------------------------- ======= -------------------------
  bool fexist( const char *filename ) {
    struct stat buffer ;
    if ( stat( filename, &buffer ) != 0 ) return false ;
    return true ;
  }

  // ------------------------- ======= ------------------------- ======= -------------------------
  bool isnumeric( const char* token, int base )
  {
    for( unsigned ichar = 0; ichar < std::string( token ).size(); ++ichar )
    {
      if( std::string( token )[ichar] != '-' &&
	  std::string( token )[ichar] != '.' &&
	  std::string( token )[ichar] != '0' &&
	  std::string( token )[ichar] != '1' &&
	  std::string( token )[ichar] != '2' &&
	  std::string( token )[ichar] != '3' &&
	  std::string( token )[ichar] != '4' &&
	  std::string( token )[ichar] != '5' &&
	  std::string( token )[ichar] != '6' &&
	  std::string( token )[ichar] != '7' &&
	  std::string( token )[ichar] != '8' &&
	  std::string( token )[ichar] != '9' )
      {
	if( std::string( token )[ichar] == 'e' && ichar < std::string( token ).size() - 1 )
	{
	  for( unsigned jchar = ichar+1; jchar < std::string( token ).size(); ++jchar )
	  {
	    if( std::string( token )[jchar] != '-' &&
		std::string( token )[jchar] != '+' &&
		std::string( token )[jchar] != '0' &&
		std::string( token )[jchar] != '1' &&
		std::string( token )[jchar] != '2' &&
		std::string( token )[jchar] != '3' &&
		std::string( token )[jchar] != '4' &&
		std::string( token )[jchar] != '5' &&
		std::string( token )[jchar] != '6' &&
		std::string( token )[jchar] != '7' &&
		std::string( token )[jchar] != '8' &&
		std::string( token )[jchar] != '9' )
	      return false;
	  }
	  return true;
	}
	else
	{
	  return false;
	}
      }
    }
    return true;
  }
  
  // ------------------------- ======= ------------------------- ======= -------------------------
  template <class T>
  bool convert(T& t, const std::string& s, 
	       std::ios_base& (*f)(std::ios_base&))
  {
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
  }

}

// ------------------------- ======= ------------------------- ======= -------------------------
TransferFunction::TransferFunction( ) : m_ndebug( 0 )
{
  
}

// ------------------------- ======= ------------------------- ======= -------------------------
TransferFunction::TransferFunction( const std::string& defn ) : m_ndebug( 0 )
{
  if( fexist( defn.c_str() ) )
  {
    std::ifstream infile( defn.c_str() );
    if( infile.fail() ) {
      std::cout << "FATAL Unable to open [" << defn << "] for reading." << std::endl;
      throw TransferFunctionReadException();
    }
    bool flag = false;
    std::pair<std::string,std::string> param_def;
    std::string token;
    unsigned itoken = 0;
    while( infile >> token )
    {
      if( itoken % 2 == 0 )
      {
	param_def.first = token;
	if( GlobalFlags::debug ) std::cout << "DEBUG " << token;
	flag = false;
      }
      else
      {
	param_def.second = token;
	if( GlobalFlags::debug ) std::cout << " = " << token;
	flag = true;
      }
      itoken += 1;
      if( flag )
      {
	double val( 0. );
	if( isnumeric(param_def.second.c_str(), 10) )
	{
	  convert<double>(val, param_def.second, std::dec);
	  if( GlobalFlags::debug ) std::cout << " (number)" << std::endl;
	  m_parameters.insert( std::pair<std::string, double>( param_def.first, val ) );
	}
	else
	{
	  if( GlobalFlags::debug ) std::cout << " (string)" << std::endl;
	  m_parameters_raw.insert( param_def );
	}
      }
      else
      {
	continue;
      }
    }
  }
  else
  {
    std::cout << "FATAL Unable to open [" << defn << "] for reading. File does not exist." << std::endl;
    throw TransferFunctionReadException();
  }
}

// ------------------------- ======= ------------------------- ======= -------------------------
double TransferFunction::parameter( const char* par, bool& flag, double def )
{
  if( m_parameters.find( std::string(par) ) != m_parameters.end() )
  {
    flag = true;
    return m_parameters[std::string(par)];
  }
  flag = false;
  return def;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double TransferFunction::parameter( const char* par, bool& flag )
{
  if( m_parameters.find( std::string(par) ) != m_parameters.end() )
  {
    flag = true;
    return m_parameters[std::string(par)];
  }
  flag = false;
  return 0.0;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double TransferFunction::parameter( const char* par )
{
  if( m_parameters.find( std::string(par) ) != m_parameters.end() )
  {
    return m_parameters[std::string(par)];
  }
  throw std::runtime_error( (std::string("parameter [") + par + std::string("] not found")).c_str() );
}

// ------------------------- ======= ------------------------- ======= -------------------------
std::string TransferFunction::str_parameter( const char* par, bool& flag )
{
  if( m_parameters_raw.find( std::string(par) ) != m_parameters_raw.end() )
  {
    flag = true;
    return m_parameters_raw[std::string(par)];
  }
  flag = false;
  return std::string("");
}

// ------------------------- ======= ------------------------- ======= -------------------------
std::string TransferFunction::str_parameter( const char* par )
{
  if( m_parameters_raw.find( std::string(par) ) != m_parameters_raw.end() )
  {
    return m_parameters_raw[std::string(par)];
  }
  throw std::runtime_error( (std::string("parameter [") + par + std::string("] not found")).c_str() );
}
