//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/HWW.hh"
#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/Fortran.hh"

#include "integrator/CubaIntegrator.hh" 

#include <iostream>
#include <cmath>
#include <vector>

extern "C"
{
  void* sgghww_(double[][4], double*, int*); // << HEFT model
  void* sgghwwsm_(double[][4], double*, int*); // << SM bb->H kludge
  void* sggyww_(double[][4], double*, int*); // << RS model spin-2
}

#include "matrix/jhudefs.icc" // << JHU ME calculations

typedef void* (*fptr_t)(double[][4], double*, int*) ;

// ------------------------- ======= ------------------------- ======= -------------------------
HWW::HWW( Integrator* integrator, 
	  double mass, 
	  int strategy,
	  TransferFunction* tfun,  bool useNarrowWidth ) :
  WWLeptonsBaseIntegrand( "HWW", integrator, tfun, (IntegrationStrategy)(strategy) ),
  m_jhuspin( 0 ), 
  m_jhucp( +1 ),
  m_jhufracqq( -1 ),
  m_useSM( false ),
  m_useRS( false ),
  m_useHEFT( true ),
  m_useJHU( false ),
  m_useNW( useNarrowWidth ), 
  m_useBoxPS( true )
{ 
  setMass( mass );
  
  SHOW_DEBUG( std::cout << "\t mH = " << mass << std::endl );
  SHOW_DEBUG( std::cout << "\t integration strategy:" << strategy << std::endl );
  
  reconfigure();
}

// ------------------------- ======= ------------------------- ======= -------------------------
void HWW::reconfigure( )
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
    ISPT	= 4;
  }
  if( m_strategy == kWMASS_6D ) 
  {
    IWP	= 0;
    IWM	= 1;
    IH	= 2;
    IY	= 3;
    ISPT	= 4;
    ISDPHI	= 5;
  }
  
  if( m_useNW && ( m_strategy == kWMASS_4D || m_strategy == kWMASS_5D || m_strategy == kWMASS_6D ) ) 
  {
    SHOW_DEBUG( std::cout
		<< "\t using narrow width approximation for mH integration, phase space dimension = "
		<< abs( static_cast<int>( m_strategy ) ) - 1 << std::endl );
    
    // shift the indices of integration if collapsing the H mass integral
    IH     = -1;
    IY     = IY - 1;
    ISPT   = ISPT > 0 ? ISPT - 1 : -1;
    ISDPHI = ISDPHI > 0 ? ISDPHI - 1 : -1;
    
    // change the number of integration dimensions
    setNumIntegrationDimensions( abs( static_cast<int>( m_strategy ) ) - 1 );
  }
}

// ------------------------- ======= ------------------------- ======= -------------------------
void HWW::getPeaks( dvectorlist& answer, const double bounds[] )
{
  
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool HWW::setKinematics( double parameters[] )
{
  
  using namespace hepstd;
  
  // only override the base class setKinematics() for the transformed integration
  // strategy, using qw+, qw+, qh, y coordinates instead of qw+, qw-, xa, xb
  
  if( static_cast<int>( m_strategy ) > 0 )
    return WWLeptonsBaseIntegrand::setKinematics( parameters );
  
  double qwp_sqr = 0;
  double qwm_sqr = 0;
  double qh_sqr  = pow( mass(), 2 );
  
  if( m_useBoxPS )
  {
    qwp_sqr = pow(wMass,2) + wMass * wWidth * tan(parameters[IWP]);
    qwm_sqr = pow(wMass,2) + wMass * wWidth * tan(parameters[IWM]);
    if( !m_useNW )
    {
      qh_sqr = pow(mass(),2) + mass() * m_width * tan(parameters[IH]);
    }
  }
  else
  {
    qwp_sqr = pow( parameters[IWP], 2 );
    qwm_sqr = pow( parameters[IWM], 2 );
    if( !m_useNW )
    {
      qh_sqr = pow( parameters[IH], 2 );
    }
  }
  if( qwp_sqr <= 0 || qwm_sqr <= 0 || qh_sqr <= 0 )
  {
    std::cout << "ERROR q < 0, need to adjust integration volume" << std::endl;
    throw std::runtime_error( "q < 0" );
  }
  
  m_kine_solutions.clear();
  
  if( m_strategy == kWMASS_4D )
  {
    partons.lp = measured.lp;
    partons.lm = measured.lm;
    qwqwqhgy<WW0jFeynDiagram>( partons, 0, 0, 
			       qh_sqr, qwp_sqr, qwm_sqr, parameters[IY],
			       m_kine_solutions );
  }
  if( m_strategy == kWMASS_5D )
  {
    partons.lp = measured.lp;
    partons.lm = measured.lm;
    partons.recoil.SetPx( parameters[ISPT] * cos(measured.recoil.Phi()) );
    partons.recoil.SetPy( parameters[ISPT] * sin(measured.recoil.Phi()) ); 
    partons.recoil.SetPz( 0 );
    partons.recoil.SetE( partons.recoil.Pt() );
    SHOW_DEBUG( std::cout << "recoil: " << partons.recoil.Px() << ", " << partons.recoil.Py() << std::endl );
    qwqwqhgy<WW0jFeynDiagram>( partons, 
			       0, 0,
			       qh_sqr, qwp_sqr, qwm_sqr, parameters[IY],
			       m_kine_solutions, partons.recoil );
  }
  if( m_strategy == kWMASS_6D )
  {
    partons.lp = measured.lp;
    partons.lm = measured.lm;
    partons.recoil.SetPx( parameters[ISPT] * cos(measured.recoil.Phi() + parameters[ISDPHI]) );
    partons.recoil.SetPy( parameters[ISPT] * sin(measured.recoil.Phi() + parameters[ISDPHI]) ); 
    partons.recoil.SetPz( 0 );
    partons.recoil.SetE( partons.recoil.Pt() );
    SHOW_DEBUG( std::cout << "recoil: " << partons.recoil.Px() << ", " << partons.recoil.Py() << std::endl );
    qwqwqhgy<WW0jFeynDiagram>( partons, 
			       0, 0,
			       qh_sqr, qwp_sqr, qwm_sqr, parameters[IY],
			       m_kine_solutions, partons.recoil );
  }
  
  for( unsigned int isol = 0; isol < m_kine_solutions.size(); ++isol )
  {
    if( m_kine_solutions[isol].weight <= 0 )      
      continue;
    SHOW_DEBUG( std::cout << "weight: " << m_kine_solutions[isol].weight << std::endl );
    if( m_strategy == kWMASS_5D || m_strategy == kWMASS_6D ) // need to boost from LAB to CM frame after solving for neutrino momenta ... 
    {
      m_kine_solutions[isol].recoil = partons.recoil;
      applyRecoilBoost( m_kine_solutions[isol] );
      SHOW_DEBUG( std::cout << "ll pT (boosted): " << (m_kine_solutions[isol].lp + m_kine_solutions[isol].lm).Pt() << std::endl );
    }    
    if( m_useBoxPS )
    {
      m_kine_solutions[isol].weight *= fabs( pow( (qwp_sqr - pow(wMass,2)), 2 ) + 
					     pow( wMass * wWidth, 2 ) ) / ( ( wMass * wWidth ) );
      m_kine_solutions[isol].weight *= fabs( pow( (qwm_sqr - pow(wMass,2)), 2 ) + 
					     pow( wMass * wWidth, 2 ) ) / ( ( wMass * wWidth ) );
      if( !m_useNW )
      {
	m_kine_solutions[isol].weight *= fabs( pow( (qh_sqr - pow(mass(),2)), 2 ) + pow( mass() * m_width, 2 ) ) / ( mass() * m_width );
      }
      else
      {
	m_kine_solutions[isol].weight *= 1.0; // < placeholder for narrow width BW factor
      }
    }
  }
  
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double HWW::getME( double array[][4] ) const
{ 
  unsigned idiagram = 0;
  double me = 0.0;
  double me_buff = 0.0;
  
  using namespace hepstd;
  
  fptr_t mgfun = &sgghww_; // < overridden by flags in following lines 
  
  if( m_useSM )
  {
    mgfun = &sgghwwsm_;
  }
  if( m_useRS )
  {
    mgfun = &sggyww_;
  }
  if( m_useHEFT )
  {
    mgfun = &sgghww_;
  }
  
  if( m_useJHU )
  {
    initParametersJHU( this->mass(), this->width(), false, true, ihelicity() );
    if( m_jhuspin == 0 )
    {
      mgfun = &sgghww_jhu_ ;
      calculateDiagram( mgfun, array, me_buff, idiagram );
      me += PDF( kgluon, kgluon ) * me_buff;
    }
    else
    {
      if( (1.0-m_jhufracqq) > TINY )
      {
	mgfun = &sggyww_jhu_ ;
	calculateDiagram( mgfun, array, me_buff, idiagram );
	me += (1.0 - m_jhufracqq) * PDF( kgluon, kgluon ) * me_buff;
      }
      if( m_jhufracqq > TINY )
      {
	mgfun = &sqqyww_jhu_ ; 
	calculateDiagram( mgfun, array, me_buff, idiagram );
	me += m_jhufracqq * (PDF( kubar, ku ) + PDF( kdbar, kd )) * me_buff;
      }
    }
  }
  else
  {
    calculateDiagram( mgfun, array, me_buff, idiagram );
    me += PDF( kgluon, kgluon ) * me_buff;
  }
  
  SHOW_DEBUG( std::cout << "matrix element: " << me << std::endl );
  
  if( std::isnan( me ) )
    return 0;
  
  if( me > 0 )
    m_num_ps_nonnull += 1;
  
  setDoClearHelicityCombinations( false );
  return me;
}

// ------------------------- ======= ------------------------- ======= -------------------------
void HWW::setMass( double mass )
{
  m_mass  = mass;
  m_width = hepstd::higgsDecayWidth( mass );
  
  initParametersJHU( this->mass(), this->width(), false, true, ihelicity() );
}

// ------------------------- ======= ------------------------- ======= -------------------------
void HWW::setWidth( double width )
{
  m_width = width;
  
  initParametersJHU( this->mass(), this->width(), false, true, ihelicity() );
}

// ------------------------- ======= ------------------------- ======= -------------------------

void HWW::getArray( const WW0jFeynDiagram& t, double array[][4] ) const
{  
  using namespace hepstd;
  
  tlvToArray( t.pa,  array[0] );
  tlvToArray( t.pb,  array[1] );
  tlvToArray( t.lp,  array[2] );
  tlvToArray( t.nul, array[3] );
  tlvToArray( t.lm,  array[4] );
  tlvToArray( t.nur, array[5] );
}
