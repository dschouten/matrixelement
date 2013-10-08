
#include "matrix/Constants.hh"
#include "matrix/Utility.hh"

#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <iostream>

/*
 * a unified interface to CTEQ and LHAPDF PDF libraries and datasets
 */

#define NUMXSAMPLES 4000
#define NUMQSAMPLES 100


class PDFGrid
{
public:

  typedef std::vector< std::vector<double> > grid_t;
  typedef std::vector<double> row_t;

  struct coord_t 
  {
  public:
    coord_t( double mlx, double mlQ )
    {
      lx = mlx;
      lQ = mlQ;
    }
    double lx;
    double lQ;
  };

  static PDFGrid* Load( const std::string& pdfn, unsigned nx = NUMXSAMPLES, unsigned nQ = NUMQSAMPLES, unsigned int isubset = 0 );

  double operator()( hepstd::PartonType t, double x, double Q ) const;
  double operator()( int t, double x, double Q ) const;

  double u( double x, double Q ) const { return (*this)( hepstd::ku, x, Q ); }
  double d( double x, double Q ) const { return (*this)( hepstd::kd, x, Q ); }
  double s( double x, double Q ) const { return (*this)( hepstd::ks, x, Q ); }
  double c( double x, double Q ) const { return (*this)( hepstd::kc, x, Q ); }

  double ubar( double x, double Q ) const { return (*this)( hepstd::kubar, x, Q ); }
  double dbar( double x, double Q ) const { return (*this)( hepstd::kdbar, x, Q ); }
  double sbar( double x, double Q ) const { return (*this)( hepstd::ksbar, x, Q ); }
  double cbar( double x, double Q ) const { return (*this)( hepstd::kcbar, x, Q ); }

  double g( double x, double Q ) const { return (*this)( hepstd::kgluon, x, Q ); }

  double get_dx() const { return exp( m_dlx ); }
  double get_dQ() const { return exp( m_dlQ ); }

  double get_xmin() const { return exp( m_lxmin ); }
  double get_xmax() const { return exp( m_lxmax ); }
  double get_Qmin() const { return exp( m_lQmin ); }
  double get_Qmax() const { return exp( m_lQmax ); }

  unsigned int get_nx( ) const { return m_nx; }
  unsigned int get_nQ( ) const { return m_nQ; }

  std::string get_pdfn( ) const { return m_pdfn; }
  unsigned int get_isub( ) const { return m_isub; }
  
private:
  static PDFGrid* m_instance;

  PDFGrid( const std::string& pdfn, unsigned nx, unsigned nQ, unsigned int isub );

  PDFGrid( ) { }
  ~PDFGrid( ) { }

  PDFGrid( const PDFGrid& );
  PDFGrid& operator=(const PDFGrid& );

  unsigned int m_nx;
  unsigned int m_nQ;

  double m_dlx;
  double m_dlQ;

  double m_lxmax;
  double m_lxmin;
  double m_lQmax;
  double m_lQmin;

  grid_t m_u_data;
  grid_t m_d_data;
  grid_t m_s_data;
  grid_t m_c_data;

  grid_t m_ubar_data;
  grid_t m_dbar_data;
  grid_t m_sbar_data;
  grid_t m_cbar_data;

  grid_t m_g_data;

  std::map<int, const grid_t*> m_grid_map;
  
  std::string m_pdfn;
  unsigned int m_isub;

  unsigned int get_ix( double x ) const;
  unsigned int get_iQ( double Q ) const;

  void init( hepstd::PartonType t, grid_t* data );

  double lookup( unsigned int ix, unsigned int iQ, const grid_t* ) const;

};
