//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/WW.hh"
#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/Fortran.hh"

#include <iostream>
#include <cmath>
#include <vector>

extern "C"
{
  void* sdbardww_(double[][4], double*, int*);
  void* sddbarww_(double[][4], double*, int*);
  void* subaruww_(double[][4], double*, int*);
  void* suubarww_(double[][4], double*, int*);
}

// ------------------------- ======= ------------------------- ======= -------------------------
WW::WW( Integrator* integrator, 
	int strategy,
	TransferFunction* tfun ) : 
  WWLeptonsBaseIntegrand( "WW", integrator, tfun, (IntegrationStrategy)strategy ), m_doStrange( false )
{ 
  SHOW_DEBUG( std::cout << "\t integration strategy:" << strategy << std::endl );
}

// ------------------------- ======= ------------------------- ======= -------------------------
WW::WW( Integrator* integrator ) : 
  WWLeptonsBaseIntegrand( "WW", integrator, 0x0, kNEUTRINO_4D ), m_doStrange( false )
{ 
  SHOW_DEBUG( std::cout << "\t integration strategy:" << kNEUTRINO_4D << std::endl );
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool WW::initialize()
{
  if( !WWLeptonsBaseIntegrand::initialize( ) )
    return false;
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WW::getME( double array[][4] ) const
{
  using namespace hepstd;

  // if( (partons.wplus() + partons.wminus()).M() < 6.0/5.0 * wMass )
  //   return 0;

  unsigned idiagram = 0;
  double me = 0.0;
  double me_buff = 0.0;
  
  calculateDiagram( &suubarww_, array, me_buff, idiagram );
  me += me_buff * PDF( ku, kubar );
  if( m_doStrange) me += me_buff * PDF( kc, kcbar );
    
  calculateDiagram( &sddbarww_, array, me_buff, idiagram );
  me += me_buff * PDF( kd, kdbar );
  if( m_doStrange) me += me_buff * PDF( ks, ksbar );
    
  calculateDiagram( &subaruww_, array, me_buff, idiagram );
  me += me_buff * PDF( kubar, ku );
  if( m_doStrange) me += me_buff * PDF( kcbar, kc );
    
  calculateDiagram( &sdbardww_, array, me_buff, idiagram );
  me += me_buff * PDF( kdbar, kd );
  if( m_doStrange) me += me_buff * PDF( ksbar, ks );
      
  SHOW_DEBUG( std::cout << "matrix element: " << me << std::endl );

  if( std::isnan( me ) )
    return 0;

  setDoClearHelicityCombinations( false );
  return me;
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WW::eventScale( double& sa, double& sb ) const
{
  sa = sb = std::sqrt( partons.total().M() ); 
}

// ------------------------- ======= ------------------------- ======= -------------------------

void WW::getArray( const WW0jFeynDiagram& t, double array[][4] ) const
{  
  using namespace hepstd;

  tlvToArray( t.pa,  array[0] );
  tlvToArray( t.pb,  array[1] );
  tlvToArray( t.lp,  array[2] );
  tlvToArray( t.lm,  array[3] );
  tlvToArray( t.nul, array[4] );
  tlvToArray( t.nur, array[5] );
}
