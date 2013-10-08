//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "integrator/Integrator.hh"

#include <iostream>
#include <stdexcept>
#include <string>

// ------------------------- ======= ------------------------- ======= -------------------------
Integrator::Integrator( std::string name ):
  m_integrand( 0 ),
  m_userdata( NULL ),
  m_name( name ),
  m_ndim( 0 ),
  m_ncomp( 1 ),
  m_epsrel( 0 ),
  m_epsabs( 0 ),
  m_flags( 13 ),
  m_mineval( 0 ),
  m_maxeval( 100000 ),
  m_timeout( 600 )
{ }

// ------------------------- ======= ------------------------- ======= -------------------------
NullIntegrator::NullIntegrator() :
  Integrator( "Null integrator" )
{ }

// ------------------------- ======= ------------------------- ======= -------------------------
void NullIntegrator::doIntegral( double retval[], double[], int*, int*,
				 double[] ) const
{
  int dim(  getNDimensions()  );
  int comp(  getNComp()  );

  // int* dimPtr = &dim;
  // int* compPtr = &comp;

  double* params = new double[dim];

  if ( static_cast<int>( m_params.size() ) != getNDimensions() )
    throw std::runtime_error( "Parameters not set in NullIntegrator!" );
  for ( int i = 0; i < dim; ++i )
  {
    params[i] = m_params[i];
  }

  std::cout << "Params[0] : " << params[0] << std::endl;
  getIntegrand()( &dim, params, &comp, retval, m_userdata );

  delete[] params;
}

// ------------------------- ======= ------------------------- ======= -------------------------
void NullIntegrator::setParam( unsigned param, double value )
{
  if ( static_cast<int>( m_params.size() ) != getNDimensions() )
    m_params.resize( getNDimensions() );

  if ( static_cast<int>( param ) >= getNDimensions() )
    throw std::runtime_error( "Invalid setting in NullIntegrator::setParam!" );

  m_params[param] = value;
}
