//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/WW1jLeptonsBaseIntegrand.hh"
#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/Fortran.hh"

#include <iostream>
#include <cmath>
#include <vector>

// ------------------------- ======= ------------------------- ======= -------------------------
WW1jLeptonsBaseIntegrand::WW1jLeptonsBaseIntegrand( const std::string& name,
						    Integrator* integrator, 
						    TransferFunction* jfun, 
						    IntegrationStrategy strategy ) :
  FeynIntegrand<WW1jFeynDiagram>( name, integrator, abs( static_cast<int>( strategy) ), 1 ), 
  m_strategy( strategy ), 
  m_tf_jet( jfun ),
  m_me_flag( false )
{ 
  configure( );
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool WW1jLeptonsBaseIntegrand::setDynamicLimits(  )
{
  double xlim[MAX_DIMENSIONS];
  double rapidity;
  
  bool flag = false;
  
  if( m_strategy == kNEUTRINO_5D || m_strategy == kWMASS_5D )
  {
    //
    // set E limits for jet based on jet resolution TF
    //
    
    rapidity = measured.jet.Eta();
    
    if( m_tf_jet->limits( "e", measured.jet.E(), &rapidity, xlim ) )
    {
      setIntegrationLimits( IE, xlim[0], xlim[1] );
      std::cout << "\tjet E limits: " << xlim[0] << " - " << xlim[1] << std::endl;
    }
    else
    {
      throw std::runtime_error( "ERROR cannot set integration limits" );
    }
    
    flag = true;
  }
  
  return flag;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool WW1jLeptonsBaseIntegrand::setKinematics( double parameters[] )
{
  using namespace hepstd;
  
  /////////////////////////////
  // \bar{\nu} = A
  // \nu       = B
  /////////////////////////////
  
  double E, rho;
  double qwp_sqr(-1), qwm_sqr(-1);
  double xa(0), xb(0);
  
  if( static_cast<int>( m_strategy ) < 0 )
  {
    qwp_sqr = pow(wMass,2) + wMass * wWidth * tan(parameters[IWP]);
    qwm_sqr = pow(wMass,2) + wMass * wWidth * tan(parameters[IWM]);
    if( qwp_sqr <= 0 || qwm_sqr <= 0 )
    {
      std::cout << "ERROR q < 0, need to adjust integration volume" << std::endl;
      throw std::runtime_error( "q < 0" );
    }  
    xa = std::sqrt( parameters[IT] * std::exp(-2*parameters[IY] ) );
    xb = parameters[IT] / xa;
  }
  
  m_kine_solutions.clear();
  
  if( m_strategy == kNEUTRINO_5D || m_strategy == kWMASS_5D )
  {
    rho = parameters[IE]; 
    partons.jet.SetPtEtaPhiM( rho / cosh( measured.jet.Eta() ), 
			      measured.jet.Eta(),
			      measured.jet.Phi(),
			      measured.jet.M() );
  }
  
  TLorentzVector recoil;
  
  switch( m_strategy ) 
  {
    case kNEUTRINO_2D : 
      partons.nur.SetPz( parameters[IPZA] );
      partons.nul.SetPz( parameters[IPZB] );
      partons.nur.SetE( partons.nur.Rho() );
      partons.nul.SetE( partons.nul.Rho() );
      break;
    case kNEUTRINO_4D : 
      E = norm( parameters[IPTA], parameters[IPZA] );
      partons.nur.SetPxPyPzE( parameters[IPTA] * cos( parameters[IPHIA] ),
			      parameters[IPTA] * sin( parameters[IPHIA] ),
			      parameters[IPZA], E );
      E = norm( -( partons.lp + partons.lm + partons.jet ).Px() - partons.nur.Px(),
		-( partons.lp + partons.lm + partons.jet ).Py() - partons.nur.Py(),
		parameters[IPZB] );
      partons.nul.SetPxPyPzE( -( partons.lp + partons.lm + partons.jet ).Px() - partons.nur.Px(),
			      -( partons.lp + partons.lm + partons.jet ).Py() - partons.nur.Py(),
			      parameters[IPZB], E );  
      break;
    case kNEUTRINO_5D : 
      E = norm( parameters[IPTA], parameters[IPZA] );
      partons.nur.SetPxPyPzE( parameters[IPTA] * cos( parameters[IPHIA] ),
			      parameters[IPTA] * sin( parameters[IPHIA] ),
			      parameters[IPZA], E );
      E = norm( -( partons.lp + partons.lm + partons.jet ).Px() - partons.nur.Px(),
		-( partons.lp + partons.lm + partons.jet ).Py() - partons.nur.Py(),
		parameters[IPZB] );
      partons.nul.SetPxPyPzE( -( partons.lp + partons.lm + partons.jet ).Px() - partons.nur.Px(),
			      -( partons.lp + partons.lm + partons.jet ).Py() - partons.nur.Py(),
			      parameters[IPZB], E );    
      break;
    case kWMASS_4D : 
      partons.lp = measured.lp;
      partons.lm = measured.lm;
      recoil = partons.jet;
      qwqwxaxb<WW1jFeynDiagram>( partons, 0, 0,
				 qwp_sqr, qwm_sqr,
				 xa, xb,
				 m_kine_solutions, recoil );
      break;
    case kWMASS_5D :
      partons.lp = measured.lp;
      partons.lm = measured.lm;
      recoil = partons.jet;
      qwqwxaxb<WW1jFeynDiagram>( partons, 0, 0,
				 qwp_sqr, qwm_sqr,
				 xa, xb,
				 m_kine_solutions, recoil );
      break;
  }
  
  for( unsigned int isol = 0; isol < m_kine_solutions.size(); ++isol )
  {
    m_kine_solutions[isol].jet = partons.jet;
    m_kine_solutions[isol].weight *= fabs( pow( (qwp_sqr - pow(wMass,2)), 2 ) + pow( wMass * wWidth, 2 ) ) / ( ( wMass * wWidth ) );
    m_kine_solutions[isol].weight *= fabs( pow( (qwm_sqr - pow(wMass,2)), 2 ) + pow( wMass * wWidth, 2 ) ) / ( ( wMass * wWidth ) );
  }
  
  return true;
}
// ------------------------- ======= ------------------------- ======= -------------------------
bool WW1jLeptonsBaseIntegrand::initialize()
{
  if( !FeynIntegrand<WW1jFeynDiagram>::initialize( ) )
    return false;
  
  hepstd::prepareCommonBlocks( );  
  
  this->setDynamicLimits( );
  
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool WW1jLeptonsBaseIntegrand::finalize()
{
  if( !FeynIntegrand<WW1jFeynDiagram>::finalize( ) )
    return false;
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WW1jLeptonsBaseIntegrand::totalTF( ) const 
{
  double parameters[2];
  
  double tf = 1.0;
  
  if( m_strategy == kNEUTRINO_5D || m_strategy == kWMASS_5D )
  {
    parameters[0] = measured.jet.E();
    parameters[1] = measured.jet.Eta();
    tf *= (*m_tf_jet)( "e", partons.jet.E(), parameters );
  }
  
  return tf;
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WW1jLeptonsBaseIntegrand::eventScale( double& sa, double& sb ) const
{
  double scale = 2.0 * hepstd::wMass; 
  if (scale < 0)
    sa = sb = 0;
  else
    sa = sb = std::sqrt( scale );
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WW1jLeptonsBaseIntegrand::phaseSpace() const
{
  if( (!m_me_flag && static_cast<int>( m_strategy ) < 0) )
  {
    return 1.0;
  }
  
  double p = partons.phaseSpace();
  
  if( static_cast<int>(m_strategy) > 0 ) // Cartesian coordinates system
  {
    if( m_strategy == kNEUTRINO_2D )
    {
      return p;
    }
    if( m_strategy == kNEUTRINO_4D )
    {
      p *= partons.nur.Pt();
      return p;
    }
    if( m_strategy == kNEUTRINO_5D )
    {
      p *= partons.nur.Pt();
      p *= pow(partons.jet.Rho(),2) * fabs(sin(partons.jet.Theta()));
      return p;
    }
  }
  
  if( static_cast<int>(m_strategy) < 0 ) // integration over W virtuality, parton x
  {     
    throw std::runtime_error( "phase space factor not implemented for kWMASS_4D" );
  }
  
  return p;
}

// ------------------------- ======= ------------------------- ======= -------------------------

bool WW1jLeptonsBaseIntegrand::isPossible( const TLorentzVector& total ) const
{
  return ( (!m_me_flag && static_cast<int>( m_strategy ) < 0) || ( FeynIntegrand<WW1jFeynDiagram>::isPossible( total ) ) );
}

// ------------------------- ======= ------------------------- ======= -------------------------

bool WW1jLeptonsBaseIntegrand::getPartonEnergies(double& ea, double& eb, const TLorentzVector& total ) const
{
  return ( (!m_me_flag && static_cast<int>( m_strategy ) < 0) || ( FeynIntegrand<WW1jFeynDiagram>::getPartonEnergies( ea, eb, total ) ) );
}

// ------------------------- ======= ------------------------- ======= -------------------------

void WW1jLeptonsBaseIntegrand::getArray( const WW1jFeynDiagram& t, double array[][4] ) const
{  
  using namespace hepstd;
  
  tlvToArray( t.pa,  array[0] );
  tlvToArray( t.pb,  array[1] );
  tlvToArray( t.lp,  array[2] );
  tlvToArray( t.lm,  array[3] );
  tlvToArray( t.nul, array[4] );
  tlvToArray( t.nur, array[5] );
  tlvToArray( t.jet, array[6] );
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WW1jLeptonsBaseIntegrand::matrixElement() const
{
  double array[7][4];
  getArray( partons, array );
  if( static_cast<int>( m_strategy ) > 0 )
  {
    SHOW_DEBUG( partons.print() );
    return getME( array );
  }
  else
  {
    TLorentzVector total;
    double me = 0.0;
    double ea, eb;
    m_me_flag = true;
    for( unsigned int isol = 0; isol < m_kine_solutions.size(); ++isol )
    {
      if( m_kine_solutions[isol].weight < 0 )
	continue;
      
      if( !isPossible( m_kine_solutions[isol].total() ) )
      {
	continue;
      }
      
      getPartonEnergies( ea, eb, m_kine_solutions[isol].total() );
      partons.pa.SetPxPyPzE( 0, 0,  ea, ea ); // this is needed for querying the PDF tables 
      partons.pb.SetPxPyPzE( 0, 0, -eb, eb );
      
      m_kine_solutions[isol].pa = partons.pa;
      m_kine_solutions[isol].pb = partons.pb;
    
      if( m_kine_solutions[isol].total().Pt() > TINY )
      {
	continue;
      }
      
      getArray( m_kine_solutions[isol], array );
      
      me += getME( array ) * m_kine_solutions[isol].weight;
    }
    m_me_flag = false;
    return me;
  }
}
