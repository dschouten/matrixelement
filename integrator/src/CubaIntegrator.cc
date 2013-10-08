//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "integrator/CubaIntegrator.hh"

#include <iostream>
#include <stdexcept>

#include <signal.h>

static volatile sig_atomic_t alarmed = 0;

void alarm_handler( int sig ) 
{ 
  alarmed = 1; 
  
  sigset_t x;
  sigemptyset( &x );
  sigaddset( &x, SIGALRM );
  sigprocmask( SIG_UNBLOCK, &x, NULL );
  
  throw std::runtime_error( "timed out" ); 
}

using std::vector;

namespace CubaIntegrators
{
  Vegas::Vegas( ) :
    Integrator( "VEGAS" ),
    m_start( 10 ),
    m_increase( 10 )
  {
    // cuba::vegasgridno_ = 1;
  }

  void Vegas::doIntegral( double returnVal[], double error[], int* fail,
			  int* neval, double prob[] ) const
  {
    cuba::Vegas( getNDimensions(), getNComp(), getIntegrand(), m_userdata, getEpsilonRel(),
		 getEpsilonAbs(), getFlags(), getRandomSeed(), getMinEval(), getMaxEval(),
		 m_start, m_increase, 0, -1, "int.data", neval, fail, returnVal, error, prob );
  }

  // ------------------------- ======= ------------------------- ======= -------------------------
  Suave::Suave() :
    Integrator( "SUAVE" ),
    m_new( 10 ),
    m_flatness( 1 )
  {}

  void Suave::doIntegral( double returnVal[], double error[], int* fail,
			  int* neval, double prob[] ) const
  {
    int dummy = 0;
    int* dummy_ptr = &dummy;
    cuba::Suave( getNDimensions(), getNComp(), getIntegrand(), m_userdata, getEpsilonRel(),
		 getEpsilonAbs(), getFlags(), getRandomSeed(), getMinEval(), getMaxEval(), m_new,
		 m_flatness, dummy_ptr, neval, fail, returnVal, error, prob );
  }
  
  // ------------------------- ======= ------------------------- ======= -------------------------
  Divonne::Divonne() :
    Integrator( "DIVONNE" ),
    m_key_par( 10 ),
    m_key_int( 10 ),
    m_key_ref( 1 ),
    m_maxpass( 100 ),
    m_border( 0 ),
    m_chisqr( 10 ),
    m_minDev( .1 ),
    m_peakFinder( 0 )
  {
    m_points = new std::vector<VecDouble>();   
  }

  void Divonne::doIntegral( double returnVal[], double error[], int* fail,
			    int* neval, double prob[] ) const
  {
    int tmp = 0;

    unsigned size = m_points->size() * getNDimensions();
    double* pointList = size != 0 ? new double[size] : NULL;
    unsigned counter = 0;

    for ( vector<vector<double> >::const_iterator it = m_points->begin();
	  it != m_points->end(); ++it )
    {
      for ( vector<double>::const_iterator jit = it->begin();
	    jit != it->end(); ++jit )
      {
	pointList[counter++] = *jit;
      }
    }

    int maxPeaks = m_peakFinder ? m_maxPeaks : 0;

    try {  
      
      std::cout << "\tIntegration timeout will be set @ " << m_timeout << " seconds ..." << std::endl;
      std::cout << "\tEvaluating..." << std::endl;  
      
      // alarmed = 0;
      // alarm( m_timeout );
      
      // signal( SIGALRM, alarm_handler );   

      // double params  = 0;
      // double results = 0;
      // std::cout << getIntegrand()( NULL, &params, NULL, &results, NULL) << std::endl;

      cuba::Divonne( getNDimensions(), getNComp(), getIntegrand(), m_userdata, getEpsilonRel(),
		     getEpsilonAbs(), getFlags(), getRandomSeed(), getMinEval(), getMaxEval(),
		     m_key_par, m_key_int, m_key_ref, m_maxpass, m_border, m_chisqr,
		     m_minDev, m_points->size(), getNDimensions(), pointList,
		     maxPeaks, m_peakFinder,
		     &tmp, neval, fail, returnVal, error, prob );
    }
    catch( std::runtime_error e )
    {
      std::cout << "\tWARNING timeout reached!" << std::endl;
    }

    // alarm( 0 );

    // if( alarmed ) 
    // {
    //   (*fail) = 1;
    // }

    if( pointList != 0x0 )
      delete[] pointList;
  }

  void Divonne::addPoint( vector<double> point )
  {
    if ( point.size() != static_cast<unsigned>( getNDimensions() ) )
      throw std::invalid_argument( "CubaInt::Divonne::addPoint: Size of vector does not match dimension!" );

    m_points->push_back( point );
  }

  // ------------------------- ======= ------------------------- ======= -------------------------
  Cuhre::Cuhre() :
    Integrator( "CUHRE" ),
    m_rule( 9 )
  {}

  void Cuhre::doIntegral( double returnVal[], double error[], int* fail,
			  int* neval, double prob[] ) const
  {
    int dummy = 0;
               
    cuba::Cuhre( getNDimensions(), getNComp(), getIntegrand(), m_userdata, getEpsilonRel(),
		 getEpsilonAbs(), getFlags(), getMinEval(), getMaxEval(),
		 m_rule, &dummy, neval, fail, returnVal, error, prob );
  }
}
