//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/DY.hh"
#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/Fortran.hh"

#include <iostream>
#include <cmath>
#include <vector>

extern "C"
{
  void* suubarzaj_(double[][4], double*, int*);
  void* sddbarzaj_(double[][4], double*, int*);
}

#define IR_CUTOFF 2.5

// ------------------------- ======= ------------------------- ======= -------------------------
DY::DY( Integrator* integrator, int strategy )
  : FeynIntegrand<DYFeynDiagram>( "DY", integrator, abs(strategy), 1 ), 
    m_strategy( static_cast<IntegrationStrategy>( strategy ) )
{ 
  SHOW_DEBUG( std::cout << "\t integration strategy:" << strategy << std::endl );
  configure();
}

// ------------------------- ======= ------------------------- ======= -------------------------
void DY::configure( )
{
  switch( m_strategy ) {
    case kNONE:
      IETA = -1;
      break;
    case kJET_1D:
      IETA = 0;
      break;
    case kJET_TT_3D:
      IETA = 0;
      IXP  = 1;
      IXM  = 2;
      break;
  }
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool DY::initialize()
{
  if( !FeynIntegrand<DYFeynDiagram>::initialize( ) )
    return false;
  
  measured.ida = hepstd::kunknown;
  measured.idb = hepstd::kunknown;
  
  hepstd::prepareCommonBlocks( );  

  double xlim[2];
  m_tauTF.limits( "", 0, 0, xlim );

  setIntegrationLimits( IXP, xlim[0], xlim[1] );
  setIntegrationLimits( IXM, xlim[0], xlim[1] );

  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool DY::finalize()
{
  return FeynIntegrand<DYFeynDiagram>::finalize( );
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool DY::setKinematics( double parameters[] )
{  
  double p, pt;

  if( m_strategy == kJET_TT_3D )
  {

    //
    // use collinear approximation implicitly here
    // where energy fraction X is in [0.5,1], i.e. the lepton
    // carries between 1/2 and all the tauon energy
    //

    p = std::sqrt( pow(measured.lp.E(),2) - pow(hepstd::tauMass,2) );
    pt = p / cosh( measured.lp.Eta() ) / parameters[IXP];
    partons.lp.SetPtEtaPhiM( pt, measured.lp.Eta(), measured.lp.Phi(), hepstd::tauMass );

    p = std::sqrt( pow(measured.lm.E(),2) - pow(hepstd::tauMass,2) );
    pt = p / cosh( measured.lm.Eta() ) / parameters[IXM];
    partons.lm.SetPtEtaPhiM( pt, measured.lm.Eta(), measured.lm.Phi(), hepstd::tauMass );

    m_xp = parameters[IXP];
    m_xm = parameters[IXM];
  }

  if( m_strategy == kJET_1D || m_strategy == kJET_TT_3D )
  {
    partons.jet.SetPx( -(partons.lp + partons.lm).Px() );
    partons.jet.SetPy( -(partons.lp + partons.lm).Py() );
    
    partons.jet.SetPtEtaPhiM( std::max( partons.jet.Pt(), IR_CUTOFF ), parameters[IETA], partons.jet.Phi(), 0.0 );
  }
  
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double DY::matrixElement() const
{
  double array[5][4];
  
  double me = 0.0;
  double me_buff = 0.0;
  
  unsigned idiagram = 0;
  
  using namespace hepstd;
  
  tlvToArray( partons.pa , array[0] );
  tlvToArray( partons.pb , array[1] );
  tlvToArray( partons.lp , array[2] );
  tlvToArray( partons.lm , array[3] );
  tlvToArray( partons.jet, array[4] );
  
  calculateDiagram( &suubarzaj_, array, me_buff, idiagram );
  me += me_buff * PDF( ku, kubar );  
  // calculateDiagram( &suubarzaj_, array, me_buff, idiagram );
  // me += me_buff * PDF( kc, kcbar ); 

  calculateDiagram( &sddbarzaj_, array, me_buff, idiagram );
  me += me_buff * PDF( kd, kdbar );
  // calculateDiagram( &sddbarzaj_, array, me_buff, idiagram );
  // me += me_buff * PDF( ks, ksbar );
  
  SHOW_DEBUG( partons.print() );
  SHOW_DEBUG( std::cout << "matrix element: " << me << std::endl );
  
  setDoClearHelicityCombinations( false );
  return me;
}

// ------------------------- ======= ------------------------- ======= -------------------------
void DY::eventScale( double& sa, double& sb ) const
{
  sa = sb = std::sqrt( hepstd::zMass );
}

// ------------------------- ======= ------------------------- ======= -------------------------
double DY::totalTF( ) const 
{
  if( m_strategy == kJET_TT_3D )
  {
    double tf = 0.;
    double parameters[2] = { 0 };

    parameters[0] = measured.lp.E();
    parameters[1] = +1;
    tf += m_tauTF( "", m_xp, parameters );

    parameters[0] = measured.lm.E();
    parameters[1] = -1;
    tf += m_tauTF( "", m_xm, parameters );
      
    return tf / 2.0;
  }

  return 1.0;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double DY::phaseSpace() const 
{
  double ps = partons.phaseSpace();  
  if( m_strategy == kJET_1D )
  {
    return ps * partons.jet.Rho() * partons.jet.Pt();
  } 
  else
  { 
    return ps;
  }
}
