//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/TT1jSimple.hh"
#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/Fortran.hh"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>

extern "C"
{
  void* sggttbar_(double[][4], double*, int*);
  void* sqqbarttbar_(double[][4], double*, int*);
}

// ------------------------- ======= ------------------------- ======= -------------------------
TT1jSimple::TT1jSimple( Integrator* integrator, 
			int strategy ) :
  FeynIntegrand<WW2jFeynDiagram>( "TT1jSimple", integrator, abs( strategy ), 1 ),
  m_strategy( static_cast<IntegrationStrategy>( strategy ) ),
  m_mass( hepstd::tMass ),
  m_cutoffScale( 1.0 )
{
  m_inv_jet_tf = NULL;
  m_obs_jet_tf = NULL;
  
  m_gg_only = true;
  
  configure();
}

// ------------------------- ======= ------------------------- ======= -------------------------
TT1jSimple::TT1jSimple( Integrator* integrator, 
			int strategy,
			double mass ) :
  FeynIntegrand<WW2jFeynDiagram>( "TT1jSimple", integrator, abs( strategy ), 1 ),
  m_strategy( static_cast<IntegrationStrategy>( strategy ) ),
  m_mass( mass ),
  m_cutoffScale( 1.0 )
{
  m_inv_jet_tf = NULL;
  m_obs_jet_tf = NULL;
  
  m_gg_only = true;
  
  configure();  
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool TT1jSimple::initialize()
{
  if( !FeynIntegrand<WW2jFeynDiagram>::initialize( ) )
    return false;
  
  measured.ida = hepstd::kunknown;
  measured.idb = hepstd::kunknown;
  
  hepstd::prepareCommonBlocks( );  
  
  this->setDynamicLimits( );
  
  return true;
}


// ------------------------- ======= ------------------------- ======= -------------------------
bool TT1jSimple::finalize()
{
  return FeynIntegrand<WW2jFeynDiagram>::finalize( );
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool TT1jSimple::setDynamicLimits( )
{
  double xlim[MAX_DIMENSIONS];
  
  double parameter( 0 );
  
  bool flag = false;
  
  // 
  // note: it is entirely valid if the efficiency limits are not defined 
  //       because these limits only make sense when ignoring detector
  //       fiducial acceptance (i.e., setting them implies a simplification
  //       in the integration volume)
  //
  
  if( m_inv_jet_tf != NULL )
  {
    if( m_inv_jet_tf->limits( "eff", m_cutoffScale, &parameter, xlim ) )
    {
      // find limits CUTOFF < pT < X where X is the pT at which 1 - eff(X,eta) < 0.1e-1
      if( xlim[1] < m_cutoffScale )
	throw std::runtime_error( "ERROR cannot set integration limits XL > XH" );
      setIntegrationLimits( IPTJB, m_cutoffScale, xlim[1] ); 
      
      SHOW_DEBUG( std::cout << "jet pT limits:" << m_cutoffScale << ":" << xlim[1] << std::endl );
    } 
    else
    {
      throw std::runtime_error( "ERROR cannot set integration limits" );
    }
    flag = true;
  }
  
  if( m_strategy == kJET_8D )
  { 
    // find limits L < pT < R where L, R define 3 \sigma sidebands
    parameter = measured.jeta.Eta();
    
    if( m_obs_jet_tf->limits( "e", measured.jeta.E(), &parameter, xlim ) )
    {
      setIntegrationLimits( IEJA, xlim[0], xlim[1] ); 
      SHOW_DEBUG( std::cout << "jet E limits: " << xlim[0] << " - " << xlim[1] << std::endl );
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
bool TT1jSimple::setKinematics( double parameters[] ) 
{
  double rho = 0.;
  
  if( m_strategy == kJET_8D )
  {
    rho = sqrt( pow( parameters[IEJA], 2 ) - pow( measured.jeta.M(), 2 ) );
  }
  
  double dx, dy, E;
  
  switch( m_strategy ) 
  {
    case kJET_7D:
      partons.jetb.SetPtEtaPhiM( parameters[IPTJB],
				 parameters[IETAJB],
				 parameters[IPHIJB],
				 hepstd::bMass );
      E = hepstd::norm( parameters[IPTA], parameters[IPZA] );
      partons.nur.SetPxPyPzE( parameters[IPTA] * cos( parameters[IPHIA] ),
			      parameters[IPTA] * sin( parameters[IPHIA] ),
			      parameters[IPZA], E );
      E = hepstd::norm( -(partons.lp.Px() + partons.lm.Px() + partons.jeta.Px() + partons.jetb.Px()) - partons.nur.Px(),
			-(partons.lp.Py() + partons.lm.Py() + partons.jeta.Py() + partons.jetb.Py()) - partons.nur.Py(),
			parameters[IPZB] );
      partons.nul.SetPxPyPzE( -(partons.lp.Px() + partons.lm.Px() + partons.jeta.Px() + partons.jetb.Px()) - partons.nur.Px(),
			      -(partons.lp.Py() + partons.lm.Py() + partons.jeta.Py() + partons.jetb.Py()) - partons.nur.Py(),
			      parameters[IPZB], E );
      break;
    case kJET_8D:
      partons.jetb.SetPtEtaPhiM( parameters[IPTJB],
				 parameters[IETAJB],
				 parameters[IPHIJB],
				 hepstd::bMass );
      partons.jeta.SetPtEtaPhiM( rho / cosh( measured.jeta.Eta() ), 
      				 measured.jeta.Eta(), 
      				 measured.jeta.Phi(), 
      				 measured.jeta.M() );
      dx = measured.jeta.Px() - partons.jeta.Px();
      dy = measured.jeta.Py() - partons.jeta.Py();
      E = hepstd::norm( parameters[IPTA], parameters[IPZA] );
      partons.nur.SetPxPyPzE( parameters[IPTA] * cos( parameters[IPHIA] ),
			      parameters[IPTA] * sin( parameters[IPHIA] ),
			      parameters[IPZA], E );
      E = hepstd::norm( -(partons.lp.Px() + partons.lm.Px() + partons.jeta.Px() + partons.jetb.Px()) + dx - partons.nur.Px(),
			-(partons.lp.Py() + partons.lm.Py() + partons.jeta.Py() + partons.jetb.Py()) + dy - partons.nur.Py(),
			parameters[IPZB] );
      partons.nul.SetPxPyPzE( -(partons.lp.Px() + partons.lm.Px() + partons.jeta.Px() + partons.jetb.Px()) + dx - partons.nur.Px(),
			      -(partons.lp.Py() + partons.lm.Py() + partons.jeta.Py() + partons.jetb.Py()) + dy - partons.nur.Py(),
			      parameters[IPZB], E );
      break;
  }
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
void TT1jSimple::eventScale( double& sa, double& sb ) const
{
  double scale = partons.sHat();
  if (scale < 0)
    sa = sb = 0;
  else
    sa = sb = std::sqrt( scale );
}

// ------------------------- ======= ------------------------- ======= -------------------------
double TT1jSimple::totalTF() const
{
  double parameters[2];
  
  double tf = 1.0;
  
  if( m_inv_jet_tf != NULL )
  {
    parameters[0] = partons.jetb.Eta(); 
    tf *= (*m_inv_jet_tf)( "eff", partons.jetb.Pt(), parameters ); 
  }
  
  if( m_strategy == kJET_8D )
  { 
    parameters[0] = measured.jeta.E(); 
    parameters[1] = measured.jeta.Eta();
    tf *= (*m_obs_jet_tf)( "e", partons.jeta.E(), parameters );
  }
  
  return tf;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double TT1jSimple::phaseSpace() const
{
  double p = partons.phaseSpace();
  switch( m_strategy ) 
  {
    case kJET_7D:
      p *= partons.nur.Pt();
      p *= partons.jetb.E() * partons.jetb.Pt();
      return p;
      break;
    case kJET_8D:
      p *= partons.nur.Pt();
      p *= partons.jetb.E() * partons.jetb.Pt();
      p *= pow( partons.jeta.E(), 2 ) * fabs(sin(partons.jeta.Theta()));
      return p;
      break;
  }
  return p;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double TT1jSimple::matrixElement() const
{
  using namespace hepstd;
  
  unsigned idiagram = 0;
  
  double me = 0.0;
  double me_buff = 0.0;
  
  double array[8][4];
  
  tlvToArray( partons.pa,   array[0] );
  tlvToArray( partons.pb,   array[1] );
  tlvToArray( partons.lp,   array[2] );
  tlvToArray( partons.nul,  array[3] );
  tlvToArray( partons.lm,   array[4] );
  tlvToArray( partons.nur,  array[5] );
  tlvToArray( partons.jeta, array[6] );
  tlvToArray( partons.jetb, array[7] );
  
  calculateDiagram( &sggttbar_, array, me_buff, idiagram );
  me += me_buff * PDF( kgluon, kgluon );
  
  if( !m_gg_only )
  {
    calculateDiagram( &sqqbarttbar_, array, me_buff, idiagram );
    me += me_buff * PDF( ku, kubar );
    me += me_buff * PDF( kd, kdbar );
    
    //
    // note: d/dbar initial state truly can't be trivially added in 
    //       with different PDF weight b/c s-channel qq->Z->tt 
    //       diagram is different for up/down quarks
    //       however, expect this to be small since Z is far off-shell
    //       so qq->g->tt surely dominates ... (?)
    //
  }
  
  SHOW_DEBUG( partons.print() );
  SHOW_DEBUG( std::cout << "matrix element: " << me << std::endl );

  setDoClearHelicityCombinations( false );  
  return me;
}
