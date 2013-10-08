//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/HWW1j.hh"
#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/Fortran.hh"

#include <iostream>
#include <cmath>
#include <vector>

extern "C"
{
  void* sgghwwj_(double[][4], double*, int*); ///< HEFT gg -> Hg
  void* sgqhwwj_(double[][4], double*, int*); ///< HEFT gq -> Hq
  void* sgqbarhwwj_(double[][4], double*, int*); ///< HEFT gqbar -> Hqbar
}

// ------------------------- ======= ------------------------- ======= -------------------------
HWW1j::HWW1j( Integrator* integrator, 
	      double mass, 
	      int strategy,
	      TransferFunction* jfun, bool useNarrowWidth = false ) :
  WW1jLeptonsBaseIntegrand( "HWW1j", integrator, jfun, (IntegrationStrategy)strategy ),
  m_mass( mass ),
  m_doStrange( false ),
  m_useNW( useNarrowWidth )
{ 
  m_width = hepstd::higgsDecayWidth( m_mass );
  
  SHOW_DEBUG( std::cout << "\tmH = " << m_mass << std::endl );
  SHOW_DEBUG( std::cout << "\tintegration strategy:" << strategy << std::endl );
  reconfigure();
}

// ------------------------- ======= ------------------------- ======= -------------------------
void HWW1j::reconfigure( )
{
  if( m_strategy == kWMASS_4D )
  {
    IWP	= 0;
    IWM	= 1;
    IH	= 2;
    IY	= 3;
  }
  if( m_strategy == kWMASS_5D )
  {
    IWP	= 0;
    IWM	= 1;
    IH	= 2;
    IY	= 3;
    IE  = 4;
  }
  
  if( m_useNW && ( m_strategy == kWMASS_4D || m_strategy == kWMASS_5D ) ) 
  {
    SHOW_DEBUG( std::cout
		<< "\t using narrow width approximation for mH integration, phase space dimension = "
		<< abs( static_cast<int>( m_strategy ) ) - 1 << std::endl );
    
    // shift the indices of integration if collapsing the H mass integral
    IH     = -1;
    IY     = IY - 1;
    IE     = IE > 0 ? IE - 1 : IE;
    
    // change the number of integration dimensions
    setNumIntegrationDimensions( abs( static_cast<int>( m_strategy ) ) - 1 );
  }
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool HWW1j::initialize()
{
  if( !FeynIntegrand<WW1jFeynDiagram>::initialize( ) )
    return false;
  
  hepstd::prepareCommonBlocks( mass(), hepstd::tMass, width() );
  
  setDynamicLimits( );
  
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool HWW1j::setKinematics( double parameters[] )
{
  using namespace hepstd;
  
  if( static_cast<int>( m_strategy ) > 0 )
    return WW1jLeptonsBaseIntegrand::setKinematics( parameters );
  
  double qwp_sqr = 0;
  double qwm_sqr = 0;
  double qh_sqr = pow(mass(),2);
  
  if( static_cast<int>( m_strategy ) < 0 )
  {
    qwp_sqr = pow(wMass,2) + wMass * wWidth * tan(parameters[IWP]);
    qwm_sqr = pow(wMass,2) + wMass * wWidth * tan(parameters[IWM]);
    
    if( !m_useNW )
    {
      qh_sqr = pow(mass(),2) + mass() * width() * tan(parameters[IH]);
    }
    
    if( qwp_sqr <= 0 || qwm_sqr <= 0 || qh_sqr <= 0 )
    {
      std::cout << "ERROR q < 0, need to adjust integration volume" << std::endl;
      throw std::runtime_error( "q < 0" );
    }
  }
  
  TLorentzVector recoil;
  
  double rho = 0;
  if( m_strategy == kNEUTRINO_5D || m_strategy == kWMASS_5D )
  {
    rho = parameters[IE]; 
  }
  
  m_kine_solutions.clear();
  
  if( m_strategy == kWMASS_4D )
  {
    recoil = partons.jet;
    partons.lp = measured.lp;
    partons.lm = measured.lm;
    qwqwqhgy<WW1jFeynDiagram>( partons, 0, 0, 
			       qh_sqr, qwp_sqr, qwm_sqr, parameters[IY],
			       m_kine_solutions, recoil ); 
  }
  if( m_strategy == kWMASS_5D )
  {
    partons.jet.SetPtEtaPhiM( rho / cosh( measured.jet.Eta() ), 
			      measured.jet.Eta(),
			      measured.jet.Phi(),
			      measured.jet.M() );
    partons.lp = measured.lp;
    partons.lm = measured.lm;
    recoil = partons.jet;
    qwqwqhgy<WW1jFeynDiagram>( partons, 0, 0,
			       qh_sqr, qwp_sqr, qwm_sqr, parameters[IY],
			       m_kine_solutions, recoil );
  }
  
  for( unsigned int isol = 0; isol < m_kine_solutions.size(); ++isol )
  {
    if( m_kine_solutions[isol].weight <= 0 )      
      continue;
    
    m_kine_solutions[isol].jet = partons.jet;
    if( m_kine_solutions[isol].total().Pt() > TINY )
    {
      m_kine_solutions[isol].weight = -1;
      continue;
    }
    
    m_kine_solutions[isol].weight *= ( pow( (qwp_sqr - pow(wMass,2)), 2 ) + pow( wMass * wWidth, 2 ) ) / ( pow( wMass * wWidth, 2 ) );
    m_kine_solutions[isol].weight *= ( pow( (qwm_sqr - pow(wMass,2)), 2 ) + pow( wMass * wWidth, 2 ) ) / ( pow( wMass * wWidth, 2 ) );
    
    if( !m_useNW )
    {
      m_kine_solutions[isol].weight *= ( pow( (qh_sqr - pow(mass(),2)), 2 ) + pow( mass() * width(), 2 ) ) / ( pow( mass() * width(), 2 ) );
    }
  }
  
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double HWW1j::getME( double array[][4] ) const
{
  double me = 0.0;
  double me_buff = 0.0;
  unsigned idiagram = 0;
  
  using namespace hepstd;
  
  calculateDiagram( &sgghwwj_, array, me_buff, idiagram );
  me += me_buff * PDF( kgluon, kgluon );
  
  // calculateDiagram( &sgqhwwj_, array, me_buff, idiagram );
  // me += me_buff * PDF( kgluon, ku );
  // me += me_buff * PDF( kgluon, kd );
  // if( m_doStrange ) me += me_buff * PDF( kgluon, ks );
  
  // calculateDiagram( &sgqbarhwwj_, array, me_buff, idiagram );
  // me += me_buff * PDF( kgluon, kubar );
  // me += me_buff * PDF( kgluon, kdbar );
  // if( m_doStrange ) me += me_buff * PDF( kgluon, ksbar );
  
  SHOW_DEBUG( std::cout << "matrix element: " << me << std::endl );
  
  if( std::isnan( me ) )
    return 0;
  
  setDoClearHelicityCombinations( false );
  return me;
}

// ------------------------- ======= ------------------------- ======= -------------------------
void HWW1j::setMass( double mass )
{
  m_mass  = mass;
  m_width = hepstd::higgsDecayWidth( mass );
}

// ------------------------- ======= ------------------------- ======= -------------------------
void HWW1j::setWidth( double width )
{
  m_width = width;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool HWW1j::isPossible( const TLorentzVector& t ) const 
{
  return ( WW1jLeptonsBaseIntegrand::isPossible( t ) );
}

// ------------------------- ======= ------------------------- ======= -------------------------

void HWW1j::getArray( const WW1jFeynDiagram& t, double array[][4] ) const
{  
  using namespace hepstd;
  
  tlvToArray( t.pa,  array[0] );
  tlvToArray( t.pb,  array[1] );
  tlvToArray( t.lp,  array[2] );
  tlvToArray( t.nul, array[3] );
  tlvToArray( t.lm,  array[4] );
  tlvToArray( t.nur, array[5] );
  tlvToArray( t.jet, array[6] );
}

