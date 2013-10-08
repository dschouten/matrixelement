//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/PhotonConversionTF.hh"

#include <iostream>
#include <stdexcept>
#include <limits>
#include <cmath>

#include <TH2.h>
#include <TFile.h>
#include <TAxis.h>

using namespace std;

PhotonConversionTF::PhotonConversionTF() : 
  TransferFunction(), m_y_fit_params( 0 )
{

}

PhotonConversionTF::PhotonConversionTF(const std::string & defn) : 
  TransferFunction( defn ), m_y_fit_params( 0 )
{
  bool IS_VERSION_0 = true;
  
  // char par_name_buffer[64];
  
  // int npar = (int) std::floor( parameter( "degree", IS_VERSION_0 ) + 0.5 );

  // // std::cout << "INFO initializing PhotonConversionTF with " 
  // // 	    << npar << "-degree fit for y" << std::endl;

  // bool READOK = true;
  // double par = 0;
  // for( unsigned i = 0; i <= (unsigned)npar; ++i )
  // {
  //   sprintf( par_name_buffer, "a%d", i );
  //   par = parameter( par_name_buffer, READOK );
  //   if( !READOK )
  //   {
  //     break;
  //   }
  //   m_y_fit_params.push_back( par );
  // } 

  // if( !IS_VERSION_0 && npar > 0 )
  // {
  //   std::string file = str_parameter( "file" );
  //   std::string hist = str_parameter( "histogram" );
    
  //   m_pt_cutoff = parameter( "pt_cut" );
    
  //   TFile* f = TFile::Open( file.c_str(), "read" );
  //   if( f == NULL or f->IsZombie() )
  //   {
  //     throw runtime_error( "bad TF file for PhotonConversionTF" );
  //   }
  //   TH2* h = dynamic_cast<TH2*>( f->Get( hist.c_str() ) );
  //   if( h == NULL )
  //   {
  //     throw runtime_error( "no data found for PhotonConversionTF" );
  //   }
    
  //   TAxis* xaxis = h->GetXaxis();
  //   TAxis* yaxis = h->GetYaxis();
    
  //   double xdn, xup, ydn, yup;
  //   data1D_t* tmp_ptr;
    
  //   m_ymax = -std::numeric_limits<double>::max();
  //   m_ymin =  std::numeric_limits<double>::max();
  //   m_xmax = -std::numeric_limits<double>::max();
  //   m_xmin =  std::numeric_limits<double>::max();
    
  //   for( unsigned ibx = 1; ibx < (unsigned)h->GetNbinsX() + 1; ++ibx )
  //   {
  //     xdn = xaxis->GetBinLowEdge( ibx );
  //     xup = xdn + xaxis->GetBinWidth( ibx );
      
  //     m_eff_map.insert( DataBin(xdn,xup), data1D_t() );
      
  //     for( unsigned iby = 1; iby < (unsigned)h->GetNbinsY() + 1; ++iby )
  //     {
  // 	ydn = yaxis->GetBinLowEdge( iby );
  // 	yup = ydn + yaxis->GetBinWidth( iby );
	
  // 	// std::cout << xdn << "," << xup << "," << ydn << "," << yup << std::endl;
	
  // 	tmp_ptr = &(m_eff_map[DataBin(xdn,xup)]);
  // 	tmp_ptr->insert( DataBin(ydn,yup), h->GetBinContent( ibx, iby ) );
	
  // 	if( ibx == 1 )
  // 	{
  // 	  if( yup > m_ymax ) m_ymax = yup;
  // 	  if( ydn < m_ymin ) m_ymin = ydn;
  // 	}
	
  // 	// std::cout << "test (" << ibx << "," << iby << "): "
  // 	// 		<< m_eff_map[(xdn + xup)/2][(ydn+yup)/2] << " =? " 
  // 	// 		<< h->GetBinContent( ibx, iby ) << std::endl;
	
  //     }
  //     if( xup > m_xmax ) m_xmax = xup;
  //     if( xdn < m_xmin ) m_xmin = xdn;
  //   }
  //   f->Close();
  // }
}

PhotonConversionTF::~PhotonConversionTF()
{
  
}

bool PhotonConversionTF::limits(const char * name, const double& val, const double par[], double vlim[])
{
  vlim[0] = 0;
  vlim[1] = 1; // << integration is simply over [0,1]
  return true;
}

/*
 * TF is p(y) m(y), where m(y) is (1 - o) where o is probability to observe subleading e+/e- from conversion
 */
double PhotonConversionTF::operator()(const char * name, const double& y, const double parameters[])
{
  double tf = 1;
  
  double e_fake  = parameters[0];   // energy of observed 'fake' electron
  double e_gamm  = e_fake / y;      // energy of gamma in gamma -> e+,e-
  double e_miss  = e_gamm - e_fake; // energy taken away by unobserved electron
  double pt_lead = 0;
  double pt_subl = 0;
  double eta     = parameters[1];

  double yy = y < 0.5 ? 1 - y : y;
  
  if( m_y_fit_params.size() == 0 )
  {    
    tf = 9.0/3.5 * (1 - 4.0/3.0*yy + 4.0/3.0*pow(yy,2));

    // if( y < 0.5 )
    // {
    //   pt_lead = e_miss / cosh( eta );
    //   pt_subl = e_fake / cosh( eta );
    // }
    // else
    // {
    //   pt_subl = e_miss / cosh( eta );
    //   pt_lead = e_fake / cosh( eta );
    // }
    
    // if( e_miss / cosh(eta) < m_pt_cutoff )
    //   return 1;

    return tf; // ignore id efficiency for sub-leading lepton
    
    // pt_subl = std::min( pt_subl, m_xmax );
    // pt_subl = std::max( pt_subl, m_xmin );
    
    // pt_lead = std::min( pt_lead, m_ymax );
    // pt_lead = std::max( pt_lead, m_ymin );
    
    // return tf * ( 1 - m_eff_map[pt_lead][pt_subl] ); // @FIXME interpolation here?
  }
  // else
  // {
  //   tf = m_y_fit_params[0];
  //   for( unsigned i = 1; i < m_y_fit_params.size(); ++i )
  //   {
  //     tf += m_y_fit_params[i] * std::pow(y, static_cast<int>(i));
  //   }
  //   return tf;
  // }
  return 0;
}
