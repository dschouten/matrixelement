//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/SystemBoostAverageTF.hh"

#include <iostream>
#include <string>
#include <math.h>
#include <stdexcept>

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TSpline.h>
#include <TROOT.h>

/**
 * 
 * Read in filename and histogram name that contains the desired pT spectrum
 *
 * notes/caveats: 
 *    - the phi TF is uniform over the range [-pi, pi]
 *    - the pT TF is evaluated from a (smoothed) spline of the histogram over
 *      the range contained by its x-axis
 *    - the TF's for pT and phi are accessed with names "pT" and "azimuth"
 *
 **/

// ------------------------- ======= ------------------------- ======= -------------------------
SystemBoostAverageTF::SystemBoostAverageTF( ) :
  TransferFunction( ),
  m_theory_fit( NULL ),
  m_theory_spl( NULL ),
  m_theory_dat( NULL )
{ }

// ------------------------- ======= ------------------------- ======= -------------------------
SystemBoostAverageTF::SystemBoostAverageTF( const std::string& defn ) : 
  TransferFunction( defn ),
  m_theory_fit( NULL ),
  m_theory_spl( NULL ),
  m_theory_dat( NULL )
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
  m_a = flag ? std::max( low, m_a ) : low;

  m_b = parameter( "max_cutoff", flag );
  m_b = flag ? std::min( high, std::max( m_a, m_b ) ) : high;

  if( has_fit )
  {
    m_theory_fit = new TF1( "HPTFIT", func.c_str(), m_a, m_b );
  }

  cwd->cd();
}

// ------------------------- ======= ------------------------- ======= -------------------------
SystemBoostAverageTF::~SystemBoostAverageTF() 
{ 
  if( m_theory_dat != NULL ) delete m_theory_dat; 
  if( m_theory_fit != NULL ) delete m_theory_fit; 
}

// ------------------------- ======= ------------------------- ======= -------------------------
bool SystemBoostAverageTF::limits( const char* name, const double& val, const double par[], double vlim[] ) 
{
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
      vlim[0] = -M_PI;
      vlim[1] =  M_PI;
      break;
    default:
      std::cout << "WARNING variable [" << name << "] not known." << std::endl;
      break;
  }
  return true;
}
  
// ------------------------- ======= ------------------------- ======= -------------------------
double SystemBoostAverageTF::operator()( const char* name, const double& x, const double par[] ) 
{

  //
  // ignore the measured recoil, and simply weight by the theory spectrum
  // for the sampled recoil
  //

  switch( name[0] )
  {
    case 'p':
      return m_theory_fit != NULL ?
	std::max( m_theory_fit->Eval( x ), 0.0 ) :
	std::max( m_theory_spl->Eval( x ), 0.0 );
      break;
    case 'a':
      return 1.0;
      break;
    default:
      std::cout << "WARNING variable [" << name << "] not known." << std::endl;
      break;
  }
  return 0.0;
}
