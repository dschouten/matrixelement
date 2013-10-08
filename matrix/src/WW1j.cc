//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/WW1j.hh"
#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/Fortran.hh"

#include <iostream>
#include <cmath>
#include <vector>

extern "C"
{
  void* sdbardwwj_(double[][4], double*, int*);
  void* sddbarwwj_(double[][4], double*, int*);
  void* sdbargwwj_(double[][4], double*, int*);
  void* sdgwwj_(double[][4], double*, int*);
  void* sgdbarwwj_(double[][4], double*, int*);
  void* sgdwwj_(double[][4], double*, int*);

  void* subaruwwj_(double[][4], double*, int*);
  void* suubarwwj_(double[][4], double*, int*);
  void* subargwwj_(double[][4], double*, int*);
  void* sugwwj_(double[][4], double*, int*);
  void* sgubarwwj_(double[][4], double*, int*);
  void* sguwwj_(double[][4], double*, int*);
}

// ------------------------- ======= ------------------------- ======= -------------------------
WW1j::WW1j( Integrator* integrator,
	    int strategy,
	    TransferFunction* jfun  ) :
  WW1jLeptonsBaseIntegrand( "WW1j", integrator, jfun, (IntegrationStrategy) strategy ),
  m_doStrange( false )
{ 
  SHOW_DEBUG( std::cout << "\t integration strategy:" << strategy << std::endl );
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool WW1j::initialize()
{
  if( !WW1jLeptonsBaseIntegrand::initialize( ) )
    return false;
  
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WW1j::getME( double array[][4] ) const
{
  double me = 0.0;
  double me_buff = 0.0;
  unsigned idiagram = 0;

  using namespace hepstd;

  // checked order ... all OK
  
  // calculateDiagram( &sugwwj_, array, me_buff, idiagram );
  // me += me_buff * PDF( ku, kgluon );

  // calculateDiagram( &subargwwj_, array, me_buff, idiagram );
  // me += me_buff * PDF( kubar, kgluon );

  calculateDiagram( &suubarwwj_, array, me_buff, idiagram );
  me += me_buff * PDF( ku, kubar ); 

  // calculateDiagram( &sdgwwj_, array, me_buff, idiagram );
  // me += me_buff * PDF( kd, kgluon );
  // if( m_doStrange ) me += me_buff * PDF( ks, kgluon );

  // calculateDiagram( &sdbargwwj_, array, me_buff, idiagram );
  // me += me_buff * PDF( kdbar, kgluon );
  // if( m_doStrange) me += me_buff * PDF( ksbar, kgluon );

  calculateDiagram( &sddbarwwj_, array, me_buff, idiagram );
  me += me_buff * PDF( kd, kdbar );
  if( m_doStrange ) me += me_buff * PDF( ks, ksbar );
  
  // calculateDiagram( &sguwwj_, array, me_buff, idiagram );
  // me += me_buff * PDF( kgluon, ku );

  // calculateDiagram( &sgubarwwj_, array, me_buff, idiagram );
  // me += me_buff * PDF( kgluon, kubar );

  calculateDiagram( &subaruwwj_, array, me_buff, idiagram );
  me += me_buff * PDF( kubar, ku ); 

  // calculateDiagram( &sgdwwj_, array, me_buff, idiagram );
  // me += me_buff * PDF( kgluon, kd );
  // if( m_doStrange ) me += me_buff * PDF( kgluon, ks );

  // calculateDiagram( &sgdbarwwj_, array, me_buff, idiagram );
  // me += me_buff * PDF( kgluon, kdbar );
  // if( m_doStrange ) me += me_buff * PDF( kgluon, ksbar );

  calculateDiagram( &sdbardwwj_, array, me_buff, idiagram );
  me += me_buff * PDF( kdbar, kd ); 
  if( m_doStrange ) me += me_buff * PDF( ksbar, ks );

  SHOW_DEBUG( std::cout << "matrix element: " << me << std::endl );

  setDoClearHelicityCombinations( false );
  return me;
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WW1j::eventScale( double& sa, double& sb ) const
{
  sa = sb = std::sqrt( 2.0 * hepstd::wMass );
}
