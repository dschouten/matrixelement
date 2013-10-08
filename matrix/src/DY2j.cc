//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/DY2j.hh"
#include "matrix/Constants.hh"
#include "matrix/TauTF.hh"
#include "matrix/JetEnergyResolutionTF.hh"
#include "matrix/Utility.hh"
#include "matrix/Fortran.hh"

#include <iostream>
#include <cmath>
#include <vector>

extern "C"
{
  void* sdddbardyjj_(double[][4], double*, int*);
  void* sdgdbardyjj_(double[][4], double*, int*);
  void* sdgddyjj_(double[][4], double*, int*);
  void* sdggdyjj_(double[][4], double*, int*);
  void* suggdyjj_(double[][4], double*, int*);
  void* sugubardyjj_(double[][4], double*, int*);
  void* sugudyjj_(double[][4], double*, int*);
  void* suuubardyjj_(double[][4], double*, int*);
}

// ------------------------- ======= ------------------------- ======= -------------------------
DY2j::DY2j( Integrator* integrator, int strategy )
  : FeynIntegrand<DY2jFeynDiagram>( "DY2j", integrator, abs(strategy), 1 ), 
    m_strategy( static_cast<IntegrationStrategy>( strategy ) ),
    m_tau_tf( 0x0 ),
    m_jet_tf( 0x0 )
{ 
  SHOW_DEBUG( std::cout << "\t integration strategy:" << strategy << std::endl );
  configure();
  
  m_tau_tf = new TauTF();
}

// ------------------------- ======= ------------------------- ======= -------------------------
void DY2j::configure( )
{
  switch( m_strategy ) {
    case kNONE:
      break;
    case k2D:
      IXP  = 0;
      IXM  = 1;
      IJETA = -1;
      IJETB = -1;
      break;
    case k4D:
      IXP  = 0;
      IXM  = 1;
      IJETA = 2;
      IJETB = 3;
      break;
  }
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool DY2j::initialize()
{
  if( !FeynIntegrand<DY2jFeynDiagram>::initialize( ) )
    return false;
  
  measured.ida = hepstd::kunknown;
  measured.idb = hepstd::kunknown;
  
  hepstd::prepareCommonBlocks( );  

  if( m_strategy != kNONE )
  {  
    double xlim[2];

    if( m_strategy == k4D )
    {
      double parameter;

      parameter = measured.jeta.Eta();    
      if( m_jet_tf->limits( "e", measured.jeta.E(), &parameter, xlim ) )
      {
	setIntegrationLimits( IJETA, xlim[0], xlim[1] ); 
	std::cout << "\tjet E limits: " << xlim[0] << " - " << xlim[1] << std::endl;
      }
      parameter = measured.jetb.Eta();    
      if( m_jet_tf->limits( "e", measured.jetb.E(), &parameter, xlim ) )
      {
	setIntegrationLimits( IJETB, xlim[0], xlim[1] ); 
	std::cout << "\tjet E limits: " << xlim[0] << " - " << xlim[1] << std::endl;
      }
    }
    
    m_tau_tf->limits( "", 0, 0, xlim );
    
    setIntegrationLimits( IXP, xlim[0], xlim[1] );
    setIntegrationLimits( IXM, xlim[0], xlim[1] );
  }
  
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool DY2j::finalize()
{
  return FeynIntegrand<DY2jFeynDiagram>::finalize( );
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool DY2j::setKinematics( double parameters[] )
{  
  if( m_strategy != kNONE )
  {
    partons.lp.SetPtEtaPhiM( std::sqrt( pow(measured.lp.E(),2) - pow(hepstd::tauMass,2) ) / cosh( measured.lp.Eta() ) / parameters[IXP],
			     measured.lp.Eta(), measured.lp.Phi(), hepstd::tauMass );
    
    partons.lm.SetPtEtaPhiM( std::sqrt( pow(measured.lm.E(),2) - pow(hepstd::tauMass,2) ) / cosh( measured.lm.Eta() ) / parameters[IXM],
			     measured.lm.Eta(), measured.lm.Phi(), hepstd::tauMass );
    
    m_xp = parameters[IXP];
    m_xm = parameters[IXM];
    
    if( m_strategy == k4D )
    {
      partons.jeta.SetPtEtaPhiE( parameters[IJETA] / cosh( measured.jeta.Eta() ),
				 measured.jeta.Eta(),
				 measured.jeta.Phi(),
				 parameters[IJETA] );

      partons.jetb.SetPtEtaPhiE( parameters[IJETB] / cosh( measured.jetb.Eta() ),
				 measured.jetb.Eta(),
				 measured.jetb.Phi(),
				 parameters[IJETB] );
    }
    
    TLorentzVector system = partons.total();
    
    double boostarry[3] = { -system.BoostVector().X(),
			    -system.BoostVector().Y(), 0 };
    
    partons.lp.Boost( boostarry );
    partons.lm.Boost( boostarry );
    partons.jeta.Boost( boostarry );
    partons.jetb.Boost( boostarry );
    
    if( partons.total().Pt() > TINY )
    {
      throw std::runtime_error( "transverse momentum imbalance" );
    }

  }
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double DY2j::matrixElement() const
{
  double array[6][4];
  
  double me = 0.0;
  double me_buff = 0.0;
  
  unsigned idiagram = 0;
  
  using namespace hepstd;
  
  tlvToArray( partons.pa , array[0] );
  tlvToArray( partons.pb , array[1] );
  tlvToArray( partons.lp , array[2] );
  tlvToArray( partons.lm , array[3] );

  for( unsigned icomb=0; icomb < 2; icomb++ )
  {
    if( icomb == 0 )
      hepstd::swap( partons.jeta, partons.jetb );

    tlvToArray( partons.jeta, array[4] );
    tlvToArray( partons.jetb, array[5] );
    
    calculateDiagram( &suggdyjj_, array, me_buff, idiagram );
    me += me_buff * PDF( kgluon, kgluon );  
    calculateDiagram( &sugubardyjj_, array, me_buff, idiagram );
    me += me_buff * PDF( kgluon, kubar );  
    calculateDiagram( &sugudyjj_, array, me_buff, idiagram );
    me += me_buff * PDF( kgluon, ku );  
    calculateDiagram( &suuubardyjj_, array, me_buff, idiagram );
    me += me_buff * PDF( ku, kubar );  
    
    calculateDiagram( &sdggdyjj_, array, me_buff, idiagram );
    me += me_buff * PDF( kgluon, kgluon );  
    calculateDiagram( &sdgdbardyjj_, array, me_buff, idiagram );
    me += me_buff * PDF( kgluon, kdbar );  
    calculateDiagram( &sdgddyjj_, array, me_buff, idiagram );
    me += me_buff * PDF( kgluon, kd );  
    calculateDiagram( &sdddbardyjj_, array, me_buff, idiagram );
    me += me_buff * PDF( kd, kdbar );  

    if( icomb == 0 )
      hepstd::swap( partons.jeta, partons.jetb );
  }

  SHOW_DEBUG( partons.print() );
  SHOW_DEBUG( std::cout << "matrix element: " << me << std::endl );

  setDoClearHelicityCombinations( false );
  return me;
}

// ------------------------- ======= ------------------------- ======= -------------------------
void DY2j::eventScale( double& sa, double& sb ) const
{
  sa = sb = std::sqrt( hepstd::zMass );
}

// ------------------------- ======= ------------------------- ======= -------------------------
double DY2j::totalTF( ) const 
{
  if( m_strategy != kNONE )
  {
    double tf = 0.;
    double parameters[2] = { 0., 0. };
    
    parameters[0] = measured.lp.E();
    parameters[1] = +1;
    tf += (*m_tau_tf)( "", m_xp, parameters );
    
    parameters[0] = measured.lm.E();
    parameters[1] = -1;
    tf += (*m_tau_tf)( "", m_xm, parameters );
    
    tf /= 2;
    
    if( m_strategy == k4D )
    {
      parameters[0] = measured.jeta.E();
      parameters[1] = measured.jeta.Eta();
      tf = tf * (*m_jet_tf)( "e", partons.jeta.E(), parameters );
      
      parameters[0] = measured.jetb.E();
      parameters[1] = measured.jetb.Eta();
      tf = tf * (*m_jet_tf)( "e", partons.jetb.E(), parameters );
    }
  }
  
  return 1.0;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double DY2j::phaseSpace() const 
{
  double ps = partons.phaseSpace();  
  if( m_strategy == k4D )
  {
    return ( ps * 
	     pow( partons.jeta.E(), 2 ) * fabs(sin(partons.jeta.Theta())) * 
	     pow( partons.jetb.E(), 2 ) * fabs(sin(partons.jetb.Theta())) );
  } 
  else
  { 
    return ps;
  }
}
