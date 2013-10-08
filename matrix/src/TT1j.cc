//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/TT1j.hh"
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
TT1j::TT1j( Integrator* integrator, 
	    int strategy ) :
  FeynIntegrand<WW2jFeynDiagram>( "TT1j", integrator, abs( strategy ), 1 ),
  m_strategy( static_cast<IntegrationStrategy>( strategy ) ),
  m_mass( hepstd::tMass ),
  m_cutoffScale( 10.0 ),
  m_me_flag( false )
{
  m_met_tf = NULL;
  m_inv_jet_tf = NULL;
  m_obs_jet_tf = NULL;
  
  configure();
  
  m_pnux.reserve(32);
  m_pnuy.reserve(32);
  m_pnuz.reserve(32);    
  
  m_pnubx.reserve(32);
  m_pnuby.reserve(32);
  m_pnubz.reserve(32);
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool TT1j::initialize()
{
  SHOW_DEBUG( std::cout << "\t initializing ... " << std::endl );
  
  if( !FeynIntegrand<WW2jFeynDiagram>::initialize( ) )
    return false;
  
  measured.ida = hepstd::kunknown;
  measured.idb = hepstd::kunknown;
  
  hepstd::prepareCommonBlocks( );  
  
  this->setDynamicLimits( );
  
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool TT1j::finalize()
{
  return FeynIntegrand<WW2jFeynDiagram>::finalize( );
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool TT1j::setDynamicLimits( )
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
      // find limits CUTOFF < pT < X where X is the pT at which 1 - eff(X,eta) < 1%
      if( xlim[1] < m_cutoffScale )
	throw std::runtime_error( "ERROR cannot set integration limits XL > XH" );
      setIntegrationLimits( IPTB, m_cutoffScale, xlim[1] ); 
      
      SHOW_DEBUG( std::cout << "\tjet pT limits: " << m_cutoffScale << " - " << xlim[1] << std::endl );
    } 
    else
    {
      throw std::runtime_error( "ERROR cannot set integration limits" );
    }
    flag = true;
  }
  
  if( m_strategy == kJET_4D || 
      m_strategy == kJET_WMASS_6D )
  { 
    // find limits L < pT < R where L, R define 3 \sigma sidebands
    parameter = measured.jeta.Eta();
    
    if( m_obs_jet_tf->limits( "e", measured.jeta.E(), &parameter, xlim ) )
    {
      setIntegrationLimits( IEA, xlim[0], xlim[1] ); 
      SHOW_DEBUG( std::cout << "\tjet E limits: " << xlim[0] << " - " << xlim[1] << std::endl );
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
bool TT1j::setKinematics( double parameters[] ) 
{
  m_wp_virt = hepstd::wMass;
  m_wm_virt = hepstd::wMass;
  
  double rho = 0.;
  
  if( m_strategy == kJET_4D ||
      m_strategy == kJET_WMASS_6D )
  {
    rho = sqrt( pow( parameters[IEA], 2 ) - pow( hepstd::bMass, 2 ) );
  }
  
  switch( m_strategy ) 
  {
    case kJET_3D:
      partons.jetb.SetPtEtaPhiM( parameters[IPTB],
				 parameters[IETAB],
				 parameters[IPHIB],
				 hepstd::bMass );
      partons.met.SetPx( -partons.lp.Px() - partons.lm.Px() - partons.jeta.Px() - partons.jetb.Px() );
      partons.met.SetPy( -partons.lp.Py() - partons.lm.Py() - partons.jeta.Py() - partons.jetb.Py() );
      break;
    case kJET_4D:
      partons.jetb.SetPtEtaPhiM( parameters[IPTB],
				 parameters[IETAB],
				 parameters[IPHIB],
				 hepstd::bMass );
      partons.jeta.SetPtEtaPhiM( rho / cosh( measured.jeta.Eta() ), 
      				 measured.jeta.Eta(), 
      				 measured.jeta.Phi(), 
      				 hepstd::bMass );
      partons.met.SetPx( -partons.lp.Px() - partons.lm.Px() - partons.jeta.Px() - partons.jetb.Px() );
      partons.met.SetPy( -partons.lp.Py() - partons.lm.Py() - partons.jeta.Py() - partons.jetb.Py() );
      break; 
    case kJET_WMASS_6D:
      partons.jetb.SetPtEtaPhiM( parameters[IPTB],
				 parameters[IETAB],
				 parameters[IPHIB],
				 hepstd::bMass );
      partons.jeta.SetPtEtaPhiM( rho / cosh( measured.jeta.Eta() ),  
      				 measured.jeta.Eta(), 
      				 measured.jeta.Phi(), 
      				 hepstd::bMass );
      partons.met.SetPx( -partons.lp.Px() - partons.lm.Px() - partons.jeta.Px() - partons.jetb.Px() );
      partons.met.SetPy( -partons.lp.Py() - partons.lm.Py() - partons.jeta.Py() - partons.jetb.Py() );
            
      m_wp_virt = sqrt(pow(hepstd::wMass,2) + hepstd::wMass * hepstd::wWidth * tan(parameters[IWP]));
      m_wm_virt = sqrt(pow(hepstd::wMass,2) + hepstd::wMass * hepstd::wWidth * tan(parameters[IWM]));
      
      if( m_wp_virt <= 0 || m_wm_virt <= 0 )
      {
        std::cerr << "ERROR q < 0, need to adjust integration volume" << std::endl;
        throw std::runtime_error( "q < 0" );
      }     
      
      break;      
  }
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double TT1j::matrixElement() const
{
  using namespace hepstd;
  
  double me = 0.0;
  unsigned idiagram = 0;
  
  m_me_flag = true;
  
  TLorentzVector total;
  double ea = 0.;
  double eb = 0.;

  for( m_icomb = 0; m_icomb < 2; ++m_icomb )
  { 
    //
    // loop over b -> t combinations
    //

    if( m_icomb == 0 )
      hepstd::swap( partons.jeta, partons.jetb );
    
    if( solveTTBar( ) )
    {        
      for( unsigned int isol = 0; isol < m_nsolutions; ++isol )
      { 
	
	//
	// loop over neutrino momenta solutions ...
	//
	
	if( m_pnux[isol]  != m_pnux[isol]  ||
	    m_pnuy[isol]  != m_pnuy[isol]  ||
	    m_pnuz[isol]  != m_pnuz[isol]  ||
	    m_pnubx[isol] != m_pnubx[isol] ||
	    m_pnuby[isol] != m_pnuby[isol] ||
	    m_pnubz[isol] != m_pnubz[isol] )
	{
	  SHOW_DEBUG( std::cerr << "WARNING caught NaN in TT1j::matrixElement() for neutrino solution" << std::endl );
	  SHOW_DEBUG( std::cerr << "     jet: " << partons.jeta.Pt()
		      << ", " << partons.jeta.Eta() 
		      << ", " << partons.jeta.Phi() << std::endl );
	  SHOW_DEBUG( std::cerr << "     jet: " << partons.jetb.Pt()
		      << ", " << partons.jetb.Eta() 
		      << ", " << partons.jetb.Phi() << std::endl );
	  continue;
	}
	
	partons.nul.SetPxPyPzE( m_pnux[isol], m_pnuy[isol], m_pnuz[isol],
				norm( m_pnux[isol], m_pnuy[isol], m_pnuz[isol] ) );
	partons.nur.SetPxPyPzE( m_pnubx[isol], m_pnuby[isol], m_pnubz[isol],
				norm( m_pnubx[isol], m_pnuby[isol], m_pnubz[isol] ) );	  
	
	STRT_TIMER(tQuarks);
	
	//
	// set the initial state (now that final state is completely solved ...)
	//
	
	partons.total( total );
	if( !isPossible( total ) )
	{
	  STOP_TIMER(tQuarks);
	  continue;
	}
	if( !getPartonEnergies( ea, eb, total ) )
	{
	  STOP_TIMER(tQuarks);
	  continue;
	}
	partons.pa.SetPxPyPzE(0, 0,  ea, ea);
	partons.pb.SetPxPyPzE(0, 0, -eb, eb);
	
	STOP_TIMER(tQuarks);
	
	if( m_strategy == kJET_WMASS_6D )
	{
	  double m_wp_virt_weight = ( pow( (m_wp_virt*m_wp_virt - pow(hepstd::wMass,2)), 2 ) + 
				      pow( hepstd::wMass * hepstd::wWidth, 2 ) ) / ( M_PI * hepstd::wMass * hepstd::wWidth );
	  double m_wm_virt_weight = ( pow( (m_wm_virt*m_wm_virt - pow(hepstd::wMass,2)), 2 ) + 
				      pow( hepstd::wMass * hepstd::wWidth, 2 ) ) / ( M_PI * hepstd::wMass * hepstd::wWidth );
	  me += getME( idiagram ) * m_wp_virt_weight * m_wm_virt_weight * phaseSpace();	    
	}
	else
	{
	  me += getME( idiagram ) * phaseSpace();
	} 
      }
    }
    
    if( m_icomb == 0 ) 
      hepstd::swap( partons.jeta, partons.jetb ); 
    
  }
  
  SHOW_DEBUG( partons.print() );
  SHOW_DEBUG( std::cout << "matrix element: " << me << std::endl );
  
  m_me_flag = false;
  
  setDoClearHelicityCombinations( false );
  return me;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool TT1j::solveTTBar( unsigned count ) const
{ 
  int cmplx(0);
  
  m_pnux.resize(0);
  m_pnuy.resize(0);
  m_pnuz.resize(0);
  
  m_pnubx.resize(0);
  m_pnuby.resize(0);
  m_pnubz.resize(0);
  
  m_cd_diff.clear();  
  
  double _lp[4]; hepstd::tlvToArray( partons.lp, _lp ); 
  double _lm[4]; hepstd::tlvToArray( partons.lm, _lm ); 
  double _b[4];  hepstd::tlvToArray( partons.jetb, _b );
  double _bb[4]; hepstd::tlvToArray( partons.jeta, _bb ); 
  double _met[2] = { partons.met.Px(), partons.met.Py() }; 
  
  m_solver.solve( _met, _b, _bb, _lp, _lm, 
		  m_wm_virt, m_wp_virt, 
		  m_mass, m_mass, 0., 0.,
		  &m_pnux, &m_pnuy, &m_pnuz, &m_pnubx, &m_pnuby, &m_pnubz,
		  &m_cd_diff, cmplx );
  
  if( m_pnux.size() != 0 )
  {
    m_nsolutions = m_pnux.size();
    return true;
  }
  
  return false;  
}

// ------------------------- ======= ------------------------- ======= -------------------------
void TT1j::eventScale( double& sa, 
		       double& sb ) const 
{
  double scale = 2 * hepstd::tMass; 
  if (scale < 0)
    sa = sb = 0;
  else
    sa = sb = std::sqrt( scale );
}

// ------------------------- ======= ------------------------- ======= -------------------------
double TT1j::totalTF() const
{
  double parameters[2];
  
  double tf = 1.0;

  if( m_inv_jet_tf != NULL )
  {
    parameters[0] = partons.jetb.Eta(); 
    tf *= (*m_inv_jet_tf)( "eff", partons.jetb.Pt(), parameters ); 
  }
  
  if( m_strategy == kJET_4D || m_strategy == kJET_WMASS_6D )
  { 
    parameters[0] = measured.jeta.E(); 
    parameters[1] = measured.jeta.Eta();
    tf *= (*m_obs_jet_tf)( "e", partons.jeta.E(), parameters );
  }
  
  return tf;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double TT1j::phaseSpace() const
{
  
  // ttbar ME has different phasespace factor, because of jet-parton permutation
  if (m_me_flag != true)
    return 1.0;
  
  double p = partons.phaseSpace();

  // evaluate determinant  of jacobian for variable transformation to top and W masses
  double j3x = 2*((partons.lp.E() + partons.jetb.E())*partons.nul.Px()/partons.nul.E() - (partons.lp.Px() + partons.jetb.Px()));
  double j3y = 2*((partons.lp.E() + partons.jetb.E())*partons.nul.Py()/partons.nul.E() - (partons.lp.Py() + partons.jetb.Py()));
  double j3z = 2*((partons.lp.E() + partons.jetb.E())*partons.nul.Pz()/partons.nul.E() - (partons.lp.Pz() + partons.jetb.Pz()));
  double j5x = 2*(partons.lp.E()*partons.nul.Px()/partons.nul.E() - partons.lp.Px());
  double j5y = 2*(partons.lp.E()*partons.nul.Py()/partons.nul.E() - partons.lp.Py());
  double j5z = 2*(partons.lp.E()*partons.nul.Pz()/partons.nul.E() - partons.lp.Pz());
      
  double j4x = 2*((partons.lm.E() + partons.jeta.E())*partons.nur.Px()/partons.nur.E() - (partons.lm.Px() + partons.jeta.Px()));
  double j4y = 2*((partons.lm.E() + partons.jeta.E())*partons.nur.Py()/partons.nur.E() - (partons.lm.Py() + partons.jeta.Py()));
  double j4z = 2*((partons.lm.E() + partons.jeta.E())*partons.nur.Pz()/partons.nur.E() - (partons.lm.Pz() + partons.jeta.Pz()));
  double j6x = 2*(partons.lm.E()*partons.nur.Px()/partons.nur.E() - partons.lm.Px());
  double j6y = 2*(partons.lm.E()*partons.nur.Py()/partons.nur.E() - partons.lm.Py());
  double j6z = 2*(partons.lm.E()*partons.nur.Pz()/partons.nur.E() - partons.lm.Pz());
      
  double DetJ = fabs(j3z*j4z*j5y*j6x - j3y*j4z*j5z*j6x - j3z*j4z*j5x*j6y + 
  		     j3x*j4z*j5z*j6y + j3z*j4y*j5x*j6z - j3z*j4x*j5y*j6z + j3y*j4x*j5z*j6z - 
  		     j3x*j4y*j5z*j6z);	
  
  // double DetJ = 1;
      
  // try to remove unphysical spikes in PS due to variable transformation
  if(DetJ < 1.0e-4)
  {
    p = 0.0;
  }
  else
  {
    p /= DetJ;
  }

  if( m_icomb == 0 ) 
    hepstd::swap( partons.jeta, partons.jetb ); 

  switch( m_strategy ) 
  {
    case kJET_3D:
      p *= partons.jetb.E() * partons.jetb.Pt();
      return p;
      break;
    case kJET_4D:
      p *= partons.jetb.E() * partons.jetb.Pt();
      p *= pow( partons.jeta.E(), 2 ) * fabs(sin(partons.jeta.Theta()));
      return p;
      break;
    case kJET_WMASS_6D:
      p *= partons.jetb.E() * partons.jetb.Pt();
      p *= pow( partons.jeta.E(), 2 ) * fabs(sin(partons.jeta.Theta()));      
      return p;
      break;      
  }

  if( m_icomb == 0 ) 
    hepstd::swap( partons.jeta, partons.jetb ); 

  return p;
}


// ------------------------- ======= ------------------------- ======= -------------------------

bool TT1j::isPossible( const TLorentzVector& total ) const
{
  return ( (!m_me_flag) || ( (total.M() <= 2 * hepstd::beamEnergy) && (std::abs(total.Pz()) < hepstd::beamEnergy) ) );
}

// ------------------------- ======= ------------------------- ======= -------------------------

bool TT1j::getPartonEnergies(double& ea, double& eb, const TLorentzVector& total ) const
{
  if( !m_me_flag )
    return true; 
  
  using hepstd::beamEnergy;
  
  double partial = std::sqrt( (total.Pz()*total.Pz()) + (total.M()*total.M()) );
  ea = (partial + total.Pz()) / 2;
  eb = (partial - total.Pz()) / 2;
  
  if (ea < 0)
    return false;
  if (ea > beamEnergy)
    return false;
  if (eb < 0)
    return false;
  if (eb > beamEnergy)
    return false;
  
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double TT1j::getME( unsigned& idiagram ) const
{
  using namespace hepstd;
  
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
  
  // FIXME for now we just consider gg initial state ( ~ 90% of production anyway )
  
  // calculateDiagram( &sqqbarttbar_, array, me_buff, idiagram );
  // me += me_buff * PDF( ku, kubar );
  // me += me_buff * PDF( kd, kdbar );
  
  //
  // note: d/dbar initial state truly can't be trivially added in 
  //       with different PDF weight b/c s-channel qq->Z->tt 
  //       diagram is different for up/down quarks
  //       however, expect this to be small since Z is far off-shell
  //       so qq->g->tt surely dominates ... (?)
  //
  
  return me;
}
