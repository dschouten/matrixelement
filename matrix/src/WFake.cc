//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/WFake.hh"
#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/Fortran.hh"

#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>

#include <TF1.h>

extern "C"
{
  void* sgubarwmj_(double[][4], double*, int*); // << W + j ME's
  void* sgdwmj_(double[][4], double*, int*);
  void* sdubarwmj_(double[][4], double*, int*);

  void* sgdbarwpj_(double[][4], double*, int*);
  void* sguwpj_(double[][4], double*, int*);
  void* sudbarwpj_(double[][4], double*, int*);

  void* sqqwpg_(double[][4], double*, int*); // << W + gamma ME's
  void* sqqwmg_(double[][4], double*, int*);
}

// ------------------------- ======= ------------------------- ======= -------------------------
WFake::WFake( Integrator* integrator, 
	IntegrationStrategy strategy,
	TransferFunction* faketf ) :
  FeynIntegrand<WFakeFeynDiagram>( "WFake", integrator, static_cast<int>( abs(strategy) ), 1 ), 
  m_me_flag( false ),
  m_wg_mode( false ),
  m_wj_mode( false ),
  m_strategy( strategy ), 
  m_fake_tf( faketf ),
  m_fake_response( 0x0 ),
  m_w_mass( hepstd::wMass ) 
{ 
  std::cout << "DEBUG create WFake ME" << std::endl;
  configure();
}

// ------------------------- ======= ------------------------- ======= -------------------------
WFake::~WFake() 
{
  if( m_fake_response != 0x0 ) delete m_fake_response;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool WFake::setDynamicLimits(  )
{
  double xlim[2] = { 0, 0 };
  double xpar[1] = { measured.fake.Eta() };

  if( m_wg_mode )
  {
    return true; // TF for Wg has fixed boundaries
  }

  if( m_wj_mode && m_strategy == k2D )
  {
    if( m_fake_tf->limits( "e", measured.fake.E(), xpar, xlim ) )
    {
      setIntegrationLimits( IFAKE, xlim[0], xlim[1] );
      SHOW_DEBUG( std::cout << "\tfake E limits: " << xlim[0] << " - " << xlim[1] << std::endl );
    }
    else
    {
      throw std::runtime_error( "ERROR cannot set integration limits" );
    }
    return true;
  }

  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool WFake::setKinematics( double parameters[] )
{
  m_w_mass     = std::sqrt( pow(hepstd::wMass,2) + hepstd::wMass * hepstd::wWidth * tan(parameters[IW]) ); // std::sqrt( parameters[IW] ); // 
  m_fake_param = 0;

  partons.real.SetPtEtaPhiM( partons.real.Pt(), partons.real.Eta(), partons.real.Phi(), 0 );
  
  if( m_strategy == k1D && m_wj_mode )
  {
    if( m_fake_response != 0x0 )
    {
      double E = measured.fake.E(); // << use the measured energy for (hypothesized) fake lepton
      E = E / (m_fake_response->Eval( E / cosh(measured.fake.Eta()) )); // << divide by simple fake response function (non-TF)

      partons.fake.SetPtEtaPhiM( E / cosh(measured.fake.Eta()),
				 measured.fake.Eta(), 
				 measured.fake.Phi(), 
				 0 );
    } 
  }
  
  if( m_strategy == k2D )
  {
    m_fake_param = parameters[IFAKE];

    double E = parameters[IFAKE]; // << use the sampled energy for (hypothesized) fake lepton
    
    if( m_wg_mode ) // << integration is over y = E[fake] / E[gamma] 
    {
      E = ( parameters[IFAKE] > TINY ? measured.fake.E() / parameters[IFAKE] : measured.fake.E() / TINY );
    }
    
    if( m_wj_mode )
    {
      if( m_fake_response != 0x0 )
      {
	E = E / (m_fake_response->Eval( E / cosh(measured.fake.Eta()) )); // << divide by simple response function
      } 
    }
    
    partons.fake.SetPtEtaPhiM( E / cosh(measured.fake.Eta()), 
			       measured.fake.Eta(), 
			       measured.fake.Phi(), measured.fake.M() );
  }  

  partons.nu.SetPx( -(partons.fake.Px() + partons.real.Px()) );
  partons.nu.SetPy( -(partons.fake.Py() + partons.real.Py()) );
  partons.nu.SetPz( 0 );
  partons.nu.SetE( partons.nu.Pt() );
  
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool WFake::initialize()
{
  if( !FeynIntegrand<WFakeFeynDiagram>::initialize( ) )
    return false;
  hepstd::prepareCommonBlocks( );  
  setDynamicLimits( );
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool WFake::finalize()
{
  if( !FeynIntegrand<WFakeFeynDiagram>::finalize( ) )
    return false;
  return true;
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WFake::eventScale( double& sa, double& sb ) const
{
  sa = sb = std::sqrt( hepstd::wMass );
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WFake::totalTF( ) const 
{
  double parameters[2];

  if( !m_me_flag )
  {
    return 1.0;
  }  
  
  double tf = 1.0;

  if( m_strategy == k2D && m_fake_tf != 0x0 )
  {
    parameters[0] = measured.fake.E();
    parameters[1] = measured.fake.Eta();
    tf *= (*m_fake_tf)( "e", m_fake_param, parameters );
  }

  return tf;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WFake::phaseSpace( ) const 
{
  if( !m_me_flag )
  {
    return 1.0;
  }  

  double ps = partons.phaseSpace();
  double j = fabs( partons.real.E() * partons.nu.Z() / partons.nu.E() - partons.real.Z() ); // 
  if( j < TINY ) 
  {
    return 0;
  }
  ps /= j; // transform from pZ to q**2
  ps *= ( ( pow( (pow(m_w_mass,2) - pow(hepstd::wMass,2)), 2 ) + 
  	    pow( hepstd::wMass * hepstd::wWidth, 2 ) ) / ( hepstd::wMass * hepstd::wWidth ) ); // transform from q**2 to t
  if( m_wg_mode && m_strategy == k2D ) // << for Wg, transform variables from E[gamma] to y = E[fake] / E[gamma]
  {
    ps *= measured.fake.E() / pow( m_fake_param, 2 );
  }

  return ps;
}

// ------------------------- ======= ------------------------- ======= -------------------------

bool WFake::isPossible( const TLorentzVector& total ) const
{
  return ( (!m_me_flag) || ( FeynIntegrand<WFakeFeynDiagram>::isPossible( total ) ) );
}

// ------------------------- ======= ------------------------- ======= -------------------------

bool WFake::getPartonEnergies(double& ea, double& eb, const TLorentzVector& total ) const
{
  return ( (!m_me_flag) || ( FeynIntegrand<WFakeFeynDiagram>::getPartonEnergies( ea, eb, total ) ) );
}

// ------------------------- ======= ------------------------- ======= -------------------------

void WFake::getArray( const WFakeFeynDiagram& t, double array[][4] ) const
{  
  using namespace hepstd;
  
  tlvToArray( t.pa,    array[0] );
  tlvToArray( t.pb,    array[1] );
  tlvToArray( t.real,  array[2] );
  tlvToArray( t.nu,    array[3] );
  tlvToArray( t.fake,  array[4] );
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WFake::matrixElement( ) const
{
  double array[5][4];
  double zsol_buff[2];

  unsigned zsol_n;

  solveZ( m_w_mass, zsol_n, zsol_buff );
  
  double me = 0.0;
  double ea = 0;
  double eb = 0;

  TLorentzVector total;
  
  m_me_flag = true;
  for( unsigned int isol = 0; isol < zsol_n; ++isol )
  {
    partons.nu.SetPx( -(partons.fake.Px() + partons.real.Px()) );
    partons.nu.SetPy( -(partons.fake.Py() + partons.real.Py()) );
    partons.nu.SetPz( zsol_buff[isol] );
    partons.nu.SetE( partons.nu.Rho() );

    if( fabs((partons.nu + partons.real).M() - m_w_mass)/m_w_mass > 0.1e-2 ||
	( (m_strategy == k1D || m_strategy == k2D) && fabs((partons.nu + partons.real + partons.fake).Pt())/partons.real.Pt() > 0.1e-2 ) )
    {
      std::cout << "ERROR " << isol << " : " << m_w_mass << " != " << (partons.nu + partons.real).M() << std::endl;
      throw std::runtime_error( "bad solution for neutrino" );
    }
    
    partons.total( total );

    if( !isPossible( total ) )
      continue;

    getPartonEnergies(ea, eb, total);
    partons.pa.SetPxPyPzE(0,0, ea,ea);
    partons.pb.SetPxPyPzE(0,0,-eb,eb);

    if( partons.pa.E() / hepstd::beamEnergy > 1 ||
	partons.pb.E() / hepstd::beamEnergy > 1 )
      continue;

    SHOW_DEBUG( partons.print() );

    getArray( partons, array );
    me += getME( array ) * phaseSpace() * totalTF();
  }
  m_me_flag = false;
  
  setDoClearHelicityCombinations( false );
  return me;
}

// ------------------------- ======= ------------------------- ======= -------------------------
double WFake::getME( double array[][4] ) const
{
  using namespace hepstd;
  
  unsigned idiagram = 0;
  double me = 0.0;
  double me_buff = 0.0;
  
  if( partons.realcharge < 0 )
  {
    if( m_wj_mode )
    {
      calculateDiagram( &sgubarwmj_, array, me_buff, idiagram );
      me += me_buff * PDF( kgluon, kubar );
      me += me_buff * PDF( kgluon, kcbar );
      
      calculateDiagram( &sgdwmj_, array, me_buff, idiagram );
      me += me_buff * PDF( kgluon, kd );
      me += me_buff * PDF( kgluon, ks );
      
      calculateDiagram( &sdubarwmj_, array, me_buff, idiagram );
      me += me_buff * PDF( kd, kubar );
      me += me_buff * PDF( ks, kcbar );
    }
    if( m_wg_mode )
    {
      calculateDiagram( &sqqwmg_, array, me_buff, idiagram );
      me += me_buff * PDF( kd, kubar );
      me += me_buff * PDF( ks, kcbar );
    }
  }
  else
  {
    if( m_wj_mode )
    {
      calculateDiagram( &sgdbarwpj_, array, me_buff, idiagram );
      me += me_buff * PDF( kgluon, kdbar );
      me += me_buff * PDF( kgluon, ksbar );
      
      calculateDiagram( &sguwpj_, array, me_buff, idiagram );
      me += me_buff * PDF( kgluon, ku );
      me += me_buff * PDF( kgluon, kc );
      
      calculateDiagram( &sudbarwpj_, array, me_buff, idiagram );
      me += me_buff * PDF( ku, kdbar );
      me += me_buff * PDF( kc, ksbar );
    }
    if( m_wg_mode )
    {
      calculateDiagram( &sqqwpg_, array, me_buff, idiagram );
      me += me_buff * PDF( ku, kdbar );
      me += me_buff * PDF( kc, ksbar );
    }
  }
  
  SHOW_DEBUG( std::cout << "matrix element: " << me << std::endl );
  
  return me;
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WFake::solveZ( double q, unsigned& zsol_n, double zsol_buff[] ) const
{
  zsol_n = 0;
  
  //
  // solve for pZ of neutrino given q**2 = (l + nu)**2, with q**2, lx, ly, lz and nux, nuy all known
  //
  
  double ea, pxa, pya, pza, ma;
  double     pxb, pyb;
  double disc;  
  
  ea  = partons.real.E();
  pxa = partons.real.Px();
  pya = partons.real.Py();
  pza = partons.real.Pz();
  ma  = partons.real.M();
  
  pxb = partons.nu.Px();
  pyb = partons.nu.Py();
  
  disc = pow(-4*pow(ma,2)*pza + 8*pxa*pxb*pza + 8*pya*pyb*pza +
	     4*pza*pow(q,2),2) -
    4*(-4*pow(ea,2) + 4*pow(pza,2))*
    (pow(ma,4) - 4*pow(ma,2)*pxa*pxb - 4*pow(ea,2)*pow(pxb,2) +
     4*pow(pxa,2)*pow(pxb,2) - 4*pow(ma,2)*pya*pyb +
     8*pxa*pxb*pya*pyb - 4*pow(ea,2)*pow(pyb,2) +
     4*pow(pya,2)*pow(pyb,2) - 2*pow(ma,2)*pow(q,2) +
     4*pxa*pxb*pow(q,2) + 4*pya*pyb*pow(q,2) + pow(q,4));
    
  if( disc >= 0 )
  {
    zsol_buff[0] =   (4*pow(ma,2)*pza - 8*pxa*pxb*pza - 8*pya*pyb*pza - 4*pza*pow(q,2) -
			sqrt(disc))/
      (2.*(-4*pow(ea,2) + 4*pow(pza,2)));
    
    zsol_buff[1] =   (4*pow(ma,2)*pza - 8*pxa*pxb*pza - 8*pya*pyb*pza - 4*pza*pow(q,2) +
			sqrt(disc))/
      (2.*(-4*pow(ea,2) + 4*pow(pza,2)));
    zsol_n = 2;
  }
}

// ------------------------- ======= ------------------------- ======= -------------------------
void WFake::setFakeResponse( const std::string& f )
{
  if( m_fake_response != 0x0 )
  {
    delete m_fake_response;
  }
  if( f.size() > 0 )
  {
    m_fake_response = new TF1( "FakeResponse", f.c_str() );
  }
  else
  {
    m_fake_response = 0x0;
  }
}

// ------------------------- ======= ------------------------- ======= -------------------------
std::string WFake::getFakeResponse( ) const
{ 
  return m_fake_response != 0x0 ? m_fake_response->GetExpFormula().Data() : ""; 
}

