//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//


#include "matrix/WWLeptonsBaseIntegrand.hh"
#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/Fortran.hh"
#include "matrix/METResolutionTF.hh"
#include "matrix/PhaseSpace.hh"

#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>

#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>

// ------------------------- ======= ------------------------- ======= -------------------------
WWLeptonsBaseIntegrand::WWLeptonsBaseIntegrand( const std::string& name, 
						Integrator* integrator,
						TransferFunction* mytf, 
						IntegrationStrategy strategy ) :
  FeynIntegrand<WW0jFeynDiagram>( name, integrator, static_cast<int>( abs(strategy) ), 1 ), 
  m_me_ready_flag( false ),
  m_strategy( strategy ), 
  m_tf( mytf ), 
  m_num_ps_nonnull( 0 ),
  m_ps_counter( 0 ),
  m_phase_space_flag( false )
{ 
  configure();
  
  psdata::ofile = NULL;
  psdata::otree = NULL;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool WWLeptonsBaseIntegrand::setDynamicLimits(  )
{
  double xlim[MAX_DIMENSIONS];
  
  if( m_strategy == kNEUTRINO_5D || m_strategy == kWMASS_5D ||
      m_strategy == kNEUTRINO_6D || m_strategy == kWMASS_6D )
  {    
    double parameters[2] = { measured.recoil.Pt(),
			     measured.recoil.Phi() };
    xlim[0] = 0;
    xlim[1] = 0;
    // set limits for system recoil pT
    //
    // selecting 'agressive' usage of recoil information, or just using 
    // theory spectrum, is done by choosing a transfer function ... the underlying code here does not change
    //
    if( m_tf->limits( "pT", measured.recoil.Pt(), parameters, xlim ) ) 
    {
      setIntegrationLimits( ISPT, xlim[0], xlim[1] ); 
      std::cout << "DEBUG limits for variable [" << ISPT << "] : " << xlim[0] << "," << xlim[1] << std::endl;
    }
    else
    {
      throw std::runtime_error( "ERROR cannot set integration limits" );
    }
    if( m_strategy == kNEUTRINO_6D || m_strategy == kWMASS_6D )
    {
      xlim[0] = 0;
      xlim[1] = 0;
      // set limits for system recoil azimuthal angle 
      if( m_tf->limits( "azimuth", hepstd::phiMPiPi( measured.recoil.Phi() ), parameters, xlim ) ) 
      {
	setIntegrationLimits( ISDPHI, xlim[0], xlim[1] ); 
	std::cout << "DEBUG limits for variable [" << ISDPHI << "] : " << measured.recoil.Phi() << " +/- "
		  << xlim[0] << "," << xlim[1] << std::endl;
      }
      else
      {
	throw std::runtime_error( "ERROR cannot set integration limits" );
      }
    }
    return true;
  }
  return false;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool WWLeptonsBaseIntegrand::setKinematics( double parameters[] )
{
  using namespace hepstd;
  
  double E = 0;  
  
  double qwp_sqr = 0;
  double qwm_sqr = 0;
  
  double xa = 0;
  double xb = 0;
  
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
      E = norm( -(partons.lp + partons.lm).Px() - partons.nur.Px(),
			-(partons.lp + partons.lm).Py() - partons.nur.Py(),
			parameters[IPZB] );
      partons.nul.SetPxPyPzE( -(partons.lp + partons.lm).Px() - partons.nur.Px(),
			      -(partons.lp + partons.lm).Py() - partons.nur.Py(),
			      parameters[IPZB], E );
      break;
    case kNEUTRINO_5D : 
      partons.lp = measured.lp;
      partons.lm = measured.lm;      
      partons.recoil.SetPx( parameters[ISPT] * cos(measured.recoil.Phi()) );
      partons.recoil.SetPy( parameters[ISPT] * sin(measured.recoil.Phi()) ); 
      partons.recoil.SetPz( 0 );      
      partons.recoil.SetE( partons.recoil.Pt() );
      E = norm( parameters[IPTA], parameters[IPZA] );
      partons.nur.SetPxPyPzE( parameters[IPTA] * cos( parameters[IPHIA] ),
			      parameters[IPTA] * sin( parameters[IPHIA] ), 
			      parameters[IPZA], E );
      
      E = norm( -(partons.lp + partons.lm + partons.recoil).Px() - partons.nur.Px(),
			-(partons.lp + partons.lm + partons.recoil).Py() - partons.nur.Py(), 
			parameters[IPZB] );
      partons.nul.SetPxPyPzE( -(partons.lp + partons.lm + partons.recoil).Px() - partons.nur.Px(),
			      -(partons.lp + partons.lm + partons.recoil).Py() - partons.nur.Py(), 
			      parameters[IPZB], E );
      applyRecoilBoost( partons );
      break;
    case kNEUTRINO_6D : 
      partons.lp = measured.lp;
      partons.lm = measured.lm;      
      partons.recoil.SetPx( parameters[ISPT] * cos(measured.recoil.Phi() + parameters[ISDPHI]) );
      partons.recoil.SetPy( parameters[ISPT] * sin(measured.recoil.Phi() + parameters[ISDPHI]) ); 
      partons.recoil.SetPz( 0 );      
      partons.recoil.SetE( partons.recoil.Pt() );
      E = norm( parameters[IPTA], parameters[IPZA] );
      partons.nur.SetPxPyPzE( parameters[IPTA] * cos( parameters[IPHIA] ),
			      parameters[IPTA] * sin( parameters[IPHIA] ), 
			      parameters[IPZA], E );
      
      E = norm( -(partons.lp + partons.lm + partons.recoil).Px() - partons.nur.Px(),
			-(partons.lp + partons.lm + partons.recoil).Py() - partons.nur.Py(), 
			parameters[IPZB] );
      partons.nul.SetPxPyPzE( -(partons.lp + partons.lm + partons.recoil).Px() - partons.nur.Px(),
			      -(partons.lp + partons.lm + partons.recoil).Py() - partons.nur.Py(), 
			      parameters[IPZB], E );
      applyRecoilBoost( partons );
      break;
    case kWMASS_4D : 
      partons.lp = measured.lp;
      partons.lm = measured.lm;
      qwqwxaxb<WW0jFeynDiagram>( partons, 0, 0, 
				 qwp_sqr, qwm_sqr,
				 xa, xb,
				 m_kine_solutions );
      break;
    case kWMASS_5D : 
      partons.lp = measured.lp;
      partons.lm = measured.lm;
      partons.recoil.SetPx( parameters[ISPT] * cos(measured.recoil.Phi()) );
      partons.recoil.SetPy( parameters[ISPT] * sin(measured.recoil.Phi()) ); 
      partons.recoil.SetPz( 0 );
      partons.recoil.SetE( partons.recoil.Pt() );
      SHOW_DEBUG( std::cout << "recoil: " << partons.recoil.Px() << ", " << partons.recoil.Py() << std::endl );
      qwqwxaxb<WW0jFeynDiagram>( partons, 
				 0, 0, 
				 qwp_sqr, qwm_sqr,
				 xa, xb, 
				 m_kine_solutions, partons.recoil );
      break;
    case kWMASS_6D : 
      partons.lp = measured.lp;
      partons.lm = measured.lm;
      partons.recoil.SetPx( parameters[ISPT] * cos(measured.recoil.Phi() + parameters[ISDPHI]) );
      partons.recoil.SetPy( parameters[ISPT] * sin(measured.recoil.Phi() + parameters[ISDPHI]) ); 
      partons.recoil.SetPz( 0 );
      partons.recoil.SetE( partons.recoil.Pt() );
      SHOW_DEBUG( std::cout << "recoil: " << partons.recoil.Px() << ", " << partons.recoil.Py() << std::endl );
      qwqwxaxb<WW0jFeynDiagram>( partons, 
				 0, 0, 
				 qwp_sqr, qwm_sqr,
				 xa, xb, 
				 m_kine_solutions, partons.recoil );
      break;
  }
  
  if( static_cast<int>( m_strategy ) < 0 )
  {
    for( unsigned int isol = 0; isol < m_kine_solutions.size(); ++isol )
    {
      if( m_kine_solutions[isol].weight < 0 )
      	continue;      
      SHOW_DEBUG( std::cout << "weight: " << m_kine_solutions[isol].weight << std::endl );      
      if( m_strategy == kWMASS_6D || m_strategy == kWMASS_5D ) // need to boost from LAB to CM frame after solving for neutrino momenta ... 
      {
	m_kine_solutions[isol].recoil = partons.recoil;
	applyRecoilBoost( m_kine_solutions[isol] );
	SHOW_DEBUG( std::cout << "ll pT (boosted): " << (m_kine_solutions[isol].lp + m_kine_solutions[isol].lm).Pt() << std::endl );
      }
      m_kine_solutions[isol].weight *= fabs( pow( (qwp_sqr - pow(wMass,2)), 2 ) + pow( wMass * wWidth, 2 ) ) / ( ( wMass * wWidth ) );
      m_kine_solutions[isol].weight *= fabs( pow( (qwm_sqr - pow(wMass,2)), 2 ) + pow( wMass * wWidth, 2 ) ) / ( ( wMass * wWidth ) );
    }
  }
  
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool WWLeptonsBaseIntegrand::initialize()
{
  if( !FeynIntegrand<WW0jFeynDiagram>::initialize( ) )
    return false;
  
  hepstd::prepareCommonBlocks( mass(), hepstd::tMass, width() );
  
  setDynamicLimits( );
  
  if( m_phase_space_flag )
  {
    bookPhaseSpaceNtuple();
  }

  m_num_ps_nonnull = 0;

  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool WWLeptonsBaseIntegrand::finalize()
{
  using namespace psdata;
  
  if( m_phase_space_flag )
  {
    TDirectory* pwd = gDirectory;
    ofile->cd();
    
    otree->Write();
    ofile->Close();
    
    pwd->cd();
  }
  
  if( !FeynIntegrand<WW0jFeynDiagram>::finalize( ) )
    return false;
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WWLeptonsBaseIntegrand::eventScale( double& sa, double& sb ) const
{
  double scale = partons.sHat();
  if (scale < 0)
    sa = sb = 0;
  else
    sa = sb = std::sqrt( scale );
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WWLeptonsBaseIntegrand::totalTF( ) const 
{
  double tf = 1.0;
  double parameters[2] = { measured.recoil.Pt(),
			   measured.recoil.Phi() };
  if( m_strategy == kNEUTRINO_5D || m_strategy == kWMASS_5D ||
      m_strategy == kNEUTRINO_6D || m_strategy == kWMASS_6D )
  {
    tf *= (*m_tf)( "pT", partons.recoil.Pt(), parameters );
    if( m_strategy == kNEUTRINO_6D || m_strategy == kWMASS_6D ) 
    {
      tf *= (*m_tf)( "angle", partons.recoil.Phi(), parameters );
    }
  }
  return tf;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WWLeptonsBaseIntegrand::phaseSpace() const 
{
  double ps = partons.phaseSpace();
  
  if( static_cast<int>( m_strategy ) < 0 )
  {
    return 1;
  }
  else // integration over neutrino momenta
  {
    if( m_strategy == kNEUTRINO_4D || 
	m_strategy == kNEUTRINO_5D || 
	m_strategy == kNEUTRINO_6D ) 
    {
      ps *= partons.nur.Pt(); 
    }
    return ps;
  }   
  return 1.0;
}

// ------------------------- ======= ------------------------- ======= -------------------------

bool WWLeptonsBaseIntegrand::isPossible( const TLorentzVector& total ) const
{
  return ( (!m_me_ready_flag && static_cast<int>( m_strategy ) < 0) || ( FeynIntegrand<WW0jFeynDiagram>::isPossible( total ) ) );
}

// ------------------------- ======= ------------------------- ======= -------------------------

bool WWLeptonsBaseIntegrand::getPartonEnergies(double& ea, double& eb, const TLorentzVector& total ) const
{
  return ( (!m_me_ready_flag && static_cast<int>( m_strategy ) < 0) || ( FeynIntegrand<WW0jFeynDiagram>::getPartonEnergies( ea, eb, total ) ) );
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WWLeptonsBaseIntegrand::matrixElement() const
{
  double array[6][4];
  double me = 0.0;
  getArray( partons, array );
  if( static_cast<int>( m_strategy ) > 0 )
  {
    SHOW_DEBUG( partons.print() );
    me = getME( array );
    
    if( m_phase_space_flag )
      dumpPhaseSpacePoint( partons, me );
    
    return me;
  }
  else
  {
    double me_buff = 0.0;
    double wgt = 0.0;
    m_me_ready_flag = true;
    for( unsigned int isol = 0; isol < m_kine_solutions.size(); ++isol )
    {      
      if( m_kine_solutions[isol].weight < 0 )
	continue;
      
      if( !isPossible( m_kine_solutions[isol].total() ) )
	continue;
      
      SHOW_DEBUG( m_kine_solutions[isol].print() );
      
      getArray( m_kine_solutions[isol], array );
      
      partons.pa = m_kine_solutions[isol].pa; // this is needed for querying the PDF tables 
      partons.pb = m_kine_solutions[isol].pb;
      
      me_buff = getME( array ); 
      wgt = m_kine_solutions[isol].weight;
      
      if( std::isnan( me_buff ) || std::isnan( wgt ) )
	continue;
      
      if( m_phase_space_flag )
	dumpPhaseSpacePoint( m_kine_solutions[isol], me_buff );
      
      me_buff *= wgt;
      
      me += me_buff;
    }
    m_me_ready_flag = false;
    
    return me;
  }
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WWLeptonsBaseIntegrand::applyRecoilBoost( WW0jFeynDiagram& lpartons, int sign ) 
{  
  int s = sign >= 0 ? 1 : -1;

  if( lpartons.recoil.Pt() < TINY )
    return;

  lpartons.recoil.SetXYZM( s * lpartons.recoil.Px(), 
			   s * lpartons.recoil.Py(),
			   ( lpartons.lp + lpartons.lm + lpartons.nul + lpartons.nur ).Pz(), 
			   ( lpartons.lp + lpartons.lm + lpartons.nul + lpartons.nur ).M() ); 
  
  double boost[3] = { lpartons.recoil.BoostVector().X(),
		      lpartons.recoil.BoostVector().Y(), 0. };
  
  lpartons.lp.Boost( boost );
  lpartons.lm.Boost( boost );
  lpartons.nul.Boost( boost );
  lpartons.nur.Boost( boost );

  boost[0] = -boost[0];
  boost[1] = -boost[1];

  lpartons.pa.Boost( boost );
  lpartons.pb.Boost( boost );
    
  if( (lpartons.lp + lpartons.lm + lpartons.nul + lpartons.nur).Pt() > TINY ||
      (lpartons.pa + lpartons.pb).Pt() > TINY )
  {
    std::cerr << "ERROR system pT is not balanced: "
	      << (lpartons.lp + lpartons.lm + lpartons.nul + lpartons.nur).Pt() << " | " 
	      << (lpartons.pa + lpartons.pb).Pt() << std::endl;
    throw std::runtime_error( "invalid PS" );
  }
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WWLeptonsBaseIntegrand::dumpPhaseSpacePoint( const WW0jFeynDiagram& p, double m ) const 
{
  using namespace psdata;
  
  iter = this->attempts();
  
  pvx	= p.nul.Px();
  pvy	= p.nul.Py();
  pvz	= p.nul.Pz();
  plx	= p.lm.Px();
  ply	= p.lm.Py();
  plz	= p.lm.Pz();
  pvbx	= p.nur.Px();
  pvby	= p.nur.Py();
  pvbz	= p.nur.Pz();
  plbx	= p.lp.Px();
  plby	= p.lp.Py();
  plbz	= p.lp.Pz();
  pwmx	= (p.lm + p.nur).Px();
  pwmy	= (p.lm + p.nur).Py();
  pwmz	= (p.lm + p.nur).Pz();
  mwm	= (p.lm + p.nur).M();
  pwpx	= (p.lp + p.nul).Px();
  pwpy	= (p.lp + p.nul).Py();
  pwpz	= (p.lp + p.nul).Pz();
  mwp	= (p.lp + p.nul).M();   
  myme	= m;
  mywgt = p.weight;  
  
  TDirectory* pwd = gDirectory;
  ofile->cd();
  
  otree->Fill();
  
  pwd->cd();
}

// ------------------------- ======= ------------------------- ======= -------------------------

void WWLeptonsBaseIntegrand::bookPhaseSpaceNtuple() 
{  
  m_ps_counter += 1;
  
  using namespace psdata;
  
  char buffer[80];
  sprintf( buffer, "ps_%04d.root", m_ps_counter );
  
  TDirectory* pwd = gDirectory;
  
  ofile = TFile::Open( buffer, "recreate" );
  ofile->cd();
  otree = new TTree( "ps", "ps" );
  otree->Branch( "iter"	, &iter , "iter/I"	);
  otree->Branch( "pvx"	, &pvx  , "pvx/D"	);
  otree->Branch( "pvy"	, &pvy  , "pvy/D"	);
  otree->Branch( "pvz"	, &pvz  , "pvz/D"	);
  otree->Branch( "pvbx"	, &pvbx , "pvbx/D"	);
  otree->Branch( "pvby"	, &pvby , "pvby/D"	);
  otree->Branch( "pvbz"	, &pvbz , "pvbz/D"	);
  otree->Branch( "plx"	, &plx  , "plx/D"	);
  otree->Branch( "ply"	, &ply  , "ply/D"	);
  otree->Branch( "plz"	, &plz  , "plz/D"	);
  otree->Branch( "plbx"	, &plbx , "plbx/D"	);
  otree->Branch( "plby"	, &plby , "plby/D"	);
  otree->Branch( "plbz"	, &plbz , "plbz/D"	);
  otree->Branch( "pwmx"	, &pwmx , "pwmx/D"	);
  otree->Branch( "pwmy"	, &pwmy , "pwmy/D"	);
  otree->Branch( "pwmz"	, &pwmz , "pwmz/D"	);
  otree->Branch( "mwm"	, &mwm  , "mwm/D"	);
  otree->Branch( "pwpx"	, &pwpx , "pwpx/D"	);
  otree->Branch( "pwpy"	, &pwpy , "pwpy/D"	);
  otree->Branch( "pwpz"	, &pwpz , "pwpz/D"	);
  otree->Branch( "mwp"	, &mwp  , "mwp/D"	);
  otree->Branch( "me"	, &myme , "me/D"	);
  otree->Branch( "wt"	, &mywgt, "wt/D"	);
  
  pwd->cd();    
}
