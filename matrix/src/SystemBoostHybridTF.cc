//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/SystemBoostHybridTF.hh"

#include <iostream>
#include <string>
#include <math.h>
#include <stdexcept>

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TSpline.h>
#include <TROOT.h>
#include <TDirectory.h>

/**
 * 
 * Read in filename and histogram name that contains the desired pT spectrum, and 
 * function that parametrizes phi resolution as function of estimated (measured!) pT
 *
 * notes/caveats: 
 *    - the functions (of 'x') must be specified without spaces in the .tf file
 *    - the integration limits for phi are reported in terms of +/- dphi instead 
 *      of absolute coordinates to avoid unnecessarily complex code to deal
 *      with phi-wrapping
 *    - the pT TF is evaluated from a (smoothed) spline of the histogram over
 *      the range contained by its x-axis
 *    - the TF's for pT and phi are accessed with names "pT" and "azimuth"
 *
 **/

// ------------------------- ======= ------------------------- ======= -------------------------
SystemBoostHybridTF::SystemBoostHybridTF( const std::string& defn ) : 
  TransferFunction( defn ),
  m_theory_fit( NULL ),
  m_theory_spl( NULL ),
  m_theory_dat( NULL ),
  m_sigma_phi_fun( 0x0 ),
  m_delta_phi_fun( 0x0 )
{
  bool has_fit;

  std::string file = str_parameter( "file_name" );
  std::string hstr = str_parameter( "hist_name" );
  std::string func = str_parameter( "hist_func", has_fit );

  TDirectory* cwd = gDirectory;

  TFile* f = TFile::Open( file.c_str(), "read" );

  gROOT->cd( );

  m_theory_dat = dynamic_cast<TH1*>( f->Get( hstr.c_str() )->Clone( "THINPUT" ) );
  m_theory_spl = new TSpline3( m_theory_dat );   
  
  f->Close();

  // allow for limit's spec'd in .tf file

  double low  = m_theory_dat->GetBinLowEdge(1);
  double high = m_theory_dat->GetBinLowEdge( m_theory_dat->GetNbinsX() + 1 );

  bool flag;

  m_a = parameter( "min_cutoff", flag );
  if( !has_fit ) m_a = flag ? std::max( low, m_a ) : low;

  m_b = parameter( "max_cutoff", flag );
  if( !has_fit ) m_b = flag ? std::min( high, std::max( m_a, m_b ) ) : high;

  std::string expr = str_parameter( "offset_phi", flag );
  if( flag ) 
  {
    m_delta_phi_fun = new TF1( "OFFSET_PHI", expr.c_str() );
  }
  
  m_sigma_phi_fun = new TF1( "SIGMA_PHI", str_parameter( "sigma_phi" ).c_str() );
  
  if( has_fit )
  {
    m_theory_fit = new TF1( "HPTFIT", func.c_str(), m_a, m_b );
  }

  cwd->cd();

  // m_tf_fun = new TF1( "TFPHIFUN", 
  // 		      ( std::string( "1.0 / ( sqrt(2 * TMath::Pi()) * [0] )" ) + "*"
  // 			std::string( "exp( -((x)*(x)) / (2.0 * [0] * [0]) ) * (1 - [1]) + [1] / (TMath::Pi())" ) ).c_str(), 0, M_PI );
  // m_tf_fun->SetParName( 0, "sigma" );
  // m_tf_fun->SetParName( 1, "delta" );
}

// ------------------------- ======= ------------------------- ======= -------------------------
SystemBoostHybridTF::~SystemBoostHybridTF() 
{ 
  if( m_theory_dat != NULL ) delete m_theory_dat; 
  if( m_theory_fit != NULL ) delete m_theory_fit; 
  if( m_theory_spl != NULL ) delete m_theory_spl; 
  if( m_sigma_phi_fun != NULL ) delete m_sigma_phi_fun;
  if( m_delta_phi_fun != NULL ) delete m_delta_phi_fun;
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool SystemBoostHybridTF::limits( const char* name, const double& val, const double par[], double vlim[] ) 
{
  double sigma;  
  switch( name[0] )
  {
    case 'p':
      if( m_theory_dat != NULL && m_theory_fit != NULL )
      {
	vlim[0] = m_a;
	vlim[1] = m_b; // integrate over the range of the input pT histogram
      }
      else
      {
	vlim[0] = vlim[1] = 0;
      }
      break;
    case 'a':
      sigma = fabs(m_sigma_phi_fun->Eval( par[0] )); 
      vlim[0] = -M_PI; // -std::min( M_PI/2., 2.0 * sigma );
      vlim[1] =  M_PI; //  std::min( M_PI/2., 2.0 * sigma );
      break;
    default:
      std::cout << "WARNING variable [" << name << "] not known." << std::endl;
      break;
  }
  return true;
}
  
// ------------------------- ======= ------------------------- ======= -------------------------
double SystemBoostHybridTF::operator()( const char* name, const double& x, const double par[] ) 
{
  //
  // integrate around +/- 3 \sigma using a Gaussian smearing of the measured
  // (event by event) recoil \phi, but use theory spectrum for pT
  //

  double sigma, delta, dx;

  switch( name[0] )
  {
    case 'p':
      return m_theory_fit != NULL ?
	std::max( m_theory_fit->Eval( x ), 0.0 ) :
	std::max( m_theory_spl->Eval( x ), 0.0 );
      break;
    case 'a':
      sigma = std::min( fabs(m_sigma_phi_fun->Eval( par[0] )), M_PI );
      delta = fabs(m_delta_phi_fun->Eval( par[0] ));
      dx    = fabs(hepstd::phiMPiPi(x-par[1]));
      return phi_kernel( dx, sigma, delta ) / phi_integral( sigma, delta );    
      break;
    default:
      std::cout << "WARNING variable [" << name << "] not known." << std::endl;
      break;
  }
  return 0.0;
}
