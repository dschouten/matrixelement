//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/PDFGrid.hh"

// ------------------------- ======= ------------------------- ======= -------------------------

#ifdef USECTEQ
extern "C"
{
  void* setctq6_(int&);
  double ctq6pdf_(int&, double&, double&);
  void* setct10_(int&);
  double ct10pdf_(int&, double&, double&);
}

//   1  CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl
//   2  CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl
//   3  CTEQ6L   Leading Order           0.118**   326** 226    cteq6l.tbl
//   4  CTEQ6L1  Leading Order           0.130**   215** 165    cteq6l1.tbl
// 100  CT10     central NLO             0.118                  ct10.00.pds  

namespace 
{
  int CT10 = 100;
  int CTEQ6M  = 1;
  int CTEQ6D  = 2;
  int CTEQ6L  = 3;
  int CTEQ6L1 = 4;

  typedef double (*pdf_fun)(int&, double&, double&);

  pdf_fun ctpdf;
}
#else
#include "lhapdf/LHAPDF.h"
#endif

// ------------------------- ======= ------------------------- ======= -------------------------

PDFGrid* PDFGrid::m_instance = NULL;

PDFGrid* PDFGrid::Load( const std::string& pdfn, unsigned nx, unsigned nQ, unsigned int isub )
{
  if( !m_instance )
  {
    m_instance = new PDFGrid( pdfn, nx, nQ, isub );
    return m_instance;
  }
  if( m_instance->get_pdfn() != pdfn ||
      m_instance->get_isub() != isub ||
      m_instance->get_nx() != nx || m_instance->get_nQ() != nQ )
  {
    delete m_instance;
    m_instance = new PDFGrid( pdfn, nx, nQ, isub );
    return m_instance;
  }
  return m_instance;
}

// ------------------------- ======= ------------------------- ======= -------------------------

PDFGrid::PDFGrid( const std::string& pdfn, unsigned nx, unsigned nQ, unsigned int isub ) :
  m_nx( nx ),
  m_nQ( nQ ),
  m_pdfn( pdfn ),
  m_isub( isub )
{
#ifndef USECTEQ
  LHAPDF::initPDFSet( pdfn.c_str(), LHAPDF::LHGRID, isubset );
  m_lxmin = std::log( LHAPDF::getXmin(isubset) );
  m_lxmax = 0;
  m_lQmin = std::log( std::sqrt( LHAPDF::getQ2min(isubset) ) );
  m_lQmax = std::log( std::sqrt( LHAPDF::getQ2max(isubset) ) );
#else
  std::cout << "\tinitializing " << pdfn << " PDF table" << std::endl;
  if( pdfn == "CTEQ6L" || pdfn == "cteq6l" ) {
    setctq6_( CTEQ6L ); ctpdf = ctq6pdf_; }
  if( pdfn == "CTEQ6D" || pdfn == "cteq6d" ) {
    setctq6_( CTEQ6D ); ctpdf = ctq6pdf_; }
  if( pdfn == "CTEQ6M" || pdfn == "cteq6m" ) {
    setctq6_( CTEQ6M ); ctpdf = ctq6pdf_; }
  if( pdfn == "CTEQ6L1" || pdfn == "cteq6l1" ) {
    setctq6_( CTEQ6L1 ); ctpdf = ctq6pdf_; }
  if( pdfn == "CT10" || pdfn == "ct10" ) {
    setct10_( CT10 ); ctpdf = ct10pdf_; }
  m_lxmin = std::log( 1.0e-6 );
  m_lxmax = 0;
  m_lQmin = std::log( 1.3e00 );
  m_lQmax = std::log( 1.0e04 );

  
#endif

  m_dlQ = (m_lQmax - m_lQmin) / m_nQ;
  m_dlx = (m_lxmax - m_lxmin) / m_nx;
 
  init( hepstd::kd, &m_d_data );
  init( hepstd::ku, &m_u_data );
  init( hepstd::ks, &m_s_data );
  init( hepstd::kc, &m_c_data );

  init( hepstd::kdbar, &m_dbar_data );
  init( hepstd::kubar, &m_ubar_data );
  init( hepstd::ksbar, &m_sbar_data );
  init( hepstd::kcbar, &m_cbar_data );
  
  init( hepstd::kgluon, &m_g_data );
}

void PDFGrid::init( hepstd::PartonType t, grid_t* data )
{
#ifdef USECTEQ
  return;
#else

  int id = static_cast<int>( t );
  double x;
  double Q;
  for( unsigned int iQ = 0; iQ < m_nQ; ++iQ )
  {
    double lQ = m_lQmin + iQ * m_dlQ;
    row_t buffer;
    for( unsigned int ix = 0; ix < m_nx; ++ix )
    {
      double lx = m_lxmin + ix * m_dlx;
#ifndef USECTEQ
      buffer.push_back( LHAPDF::xfx(exp(lx), exp(lQ), static_cast<int>( t )) / (exp(lx)) );
#else
      x = exp( lx );
      Q = exp( lQ );
      buffer.push_back( ctpdf(id, x, Q) );
#endif
    }
    data->push_back( buffer );
  }
  m_grid_map[static_cast<int>(t)] = data;
  std::cout << "\t\tloaded PDF for " << id << " with " 
            << data->size() << "x" << data->at(0).size() << " grid" << std::endl;
#endif
}

unsigned int PDFGrid::get_ix( double lx ) const
{
  if( lx <= m_lxmin ) return 0;
  if( lx >= m_lxmax ) return m_nx-1; 
  unsigned int ix = static_cast<unsigned int>( (lx - m_lxmin) / m_dlx );
  return ix < m_nx ? ix : m_nx-1;
}

unsigned int PDFGrid::get_iQ( double lQ ) const
{
  if( lQ <= m_lQmin ) return 0;
  if( lQ >= m_lQmax ) return m_nQ-1;
  unsigned int iQ = static_cast<unsigned int>( (lQ - m_lQmin) / m_dlQ );
  return iQ < m_nQ ? iQ : m_nQ-1;
}

double PDFGrid::lookup( unsigned int ix, unsigned int iQ, const grid_t* data ) const
{
  return (*data)[iQ][ix];
}  

double PDFGrid::operator()( hepstd::PartonType t, double x, double Q ) const
{
  return this->operator()( static_cast<int>(t), x, Q );
}

double PDFGrid::operator()( int t, double x, double Q ) const
{
#ifdef USECTEQ

  return ctpdf(t, x, Q);

#else
  register double lx = std::log(x);
  register double lQ = std::log(Q);

  register unsigned int ix = get_ix( lx );
  register unsigned int iQ = get_iQ( lQ );

  register double c11x = m_lxmin + m_dlx * ix;
  register double c11q = m_lQmin + m_dlQ * iQ;

  register double c22x = m_lxmin + m_dlx * (ix+1);
  register double c22q = m_lQmin + m_dlQ * (iQ+1);

  register double norm = m_dlx*m_dlQ;
  
  const grid_t* data = m_grid_map.find( t )->second;

  return ( (*data)[iQ][ix]     / (norm) * (c22x - lx)*(c22q - lQ) +
	   (*data)[iQ][ix+1]   / (norm) * (lx - c11x)*(c22q - lQ) + 
	   (*data)[iQ+1][ix]   / (norm) * (c22x - lx)*(lQ - c11q) + 
	   (*data)[iQ+1][ix+1] / (norm) * (lx - c11x)*(lQ - c11q) );
#endif
}
