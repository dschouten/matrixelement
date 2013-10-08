//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef METRESOLUTION_TF
#define METRESOLUTION_TF

#include "matrix/TransferFunction.hh"

//////////////////////////////////////////////////////////
//
// Transfer function for missing ET 
//
// - provides probability distribution for px, py given 
//   measured 
//
// Version 1.0
// - resolution based solely on scalar sum ET scaling
// 
//////////////////////////////////////////////////////////

class METResolutionTF : public TransferFunction
{
public:
  METResolutionTF( const std::string& defn );
  virtual ~METResolutionTF() { }

  virtual bool limits( const char*, const double&, const double[], double[] );
  virtual double operator()( const char*, const double&, const double[] );

private:
  double m_sigma_x_scale;
  double m_sigma_y_scale;
};

#endif

  
  
