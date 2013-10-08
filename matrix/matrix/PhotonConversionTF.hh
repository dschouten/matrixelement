//
// author: Doug Schouten <doug dot schouten at triumf dot ca>
//


#ifndef PH_CONV_TF
#define PH_CONV_TF

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <stdexcept>

#include "matrix/TransferFunction.hh"

class DataBin
{
public:
  DataBin( double x );
  DataBin( double xlo, double xhi );
  DataBin( const DataBin& );
  DataBin& operator=( const DataBin& );
  
  ~DataBin() {} 

  bool operator< ( const DataBin& other ) const { return (m_xhi < other.hi() && m_xlo < other.lo()); }
  bool operator> ( const DataBin& other ) const { return (m_xlo > other.lo() && m_xhi > other.hi()); }
  bool operator==( const DataBin& other ) const { return (m_xhi <= other.hi() && m_xlo >= other.lo()); }

  double lo() const { return m_xlo; }
  double hi() const { return m_xhi; }

private:
  double m_xlo;
  double m_xhi;
};

inline DataBin::DataBin( double x )
{
  m_xlo = x;
  m_xhi = x;
}

inline DataBin::DataBin( double xlo, double xhi )
{
  if( xlo > xhi )
  {
    throw std::runtime_error( "bins must be in order" );
  }
  m_xlo = xlo;
  m_xhi = xhi;
}

inline DataBin::DataBin( const DataBin& cl )
{
  m_xlo = cl.m_xlo;
  m_xhi = cl.m_xhi;
}

inline DataBin& DataBin::operator=( const DataBin& cl )
{
  m_xlo = cl.m_xlo;
  m_xhi = cl.m_xhi;
  return (*this);
}

//
// Class for storing binned data
//    - uses std::map for fast lookhi O(log n) instead of O(n)
//    - doesn't assume fixed bin widths
//

template<typename T>
class BinnedDataMap : public std::map<DataBin,T>
{
public:
  BinnedDataMap( ) : std::map<DataBin,T>() { }
  ~BinnedDataMap( ) { }

  const T& operator[]( double x ) const;
  T& operator[]( double x );

  const T& operator[]( DataBin x ) const;
  T& operator[]( DataBin x );

  void insert( DataBin k, T v ) 
  {
    // std::cout << "insert @ " << k.lo() << "," << k.hi() << std::endl;
    std::map<DataBin,T>::insert( std::pair<DataBin,T>( k, v ) );
  }
};

#ifndef __CINT__
template<typename T>
inline const T& BinnedDataMap<T>::operator[]( double x ) const
{ 
  typename std::map<DataBin,T>::const_iterator itr = this->find( DataBin(x) );
  if( itr == std::map<DataBin,T>::end() )
    throw std::runtime_error( "bin not found" );
  return itr->second; 
}

template<typename T>
inline T& BinnedDataMap<T>::operator[]( double x )
{ 
  typename std::map<DataBin,T>::iterator itr = this->find( DataBin(x) );
  if( itr == std::map<DataBin,T>::end() )
    throw std::runtime_error( "bin not found" );
  return itr->second; 
}

template<typename T>
inline const T& BinnedDataMap<T>::operator[]( DataBin b ) const
{ 
  typename std::map<DataBin,T>::const_iterator itr = this->find( b );
  if( itr == std::map<DataBin,T>::end() )
    throw std::runtime_error( "bin not found" );
  return itr->second; 
}

template<typename T>
inline T& BinnedDataMap<T>::operator[]( DataBin b )
{ 
  typename std::map<DataBin,T>::iterator itr = this->find( b );
  if( itr == std::map<DataBin,T>::end() )
    throw std::runtime_error( "bin not found" );
  return itr->second; 
}
#endif

typedef BinnedDataMap<double>   data1D_t;
typedef BinnedDataMap<data1D_t> data2D_t;
typedef BinnedDataMap<data2D_t> data3D_t;

//////////////////////////////////////////////////////////
//
// Photon pair production TF 
//
// - calculate the probability for a prompt photon with
//  energy k to produce electron/positron with energy E
// - multiply by probability for positron/electron with 
//   energy (k - E) to be missed (curl away in B-field)
//
////////////////////////////////////////////////////////// 

class PhotonConversionTF : public TransferFunction
{
public:
  
  PhotonConversionTF(const std::string & defn);
  
  virtual ~PhotonConversionTF();

  virtual bool limits(const char * name, const double & val, const double par[], double vlim[]);
  virtual double operator()(const char * name, const double & hypothesized_tau_energy, const double measured_lepton_energies[]);
  
private:
  PhotonConversionTF();
  
  data2D_t m_eff_map; // << efficiency map for e+/e- pair

  double m_ymax;
  double m_xmax;
  double m_ymin;
  double m_xmin;

  double m_pt_cutoff;

  std::vector<double> m_y_fit_params;
};

#endif
