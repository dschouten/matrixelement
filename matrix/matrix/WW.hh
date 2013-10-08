
//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//         
//

#ifndef WW_HH
#define WW_HH

#include "matrix/WWLeptonsBaseIntegrand.hh"
#include "matrix/TransferFunction.hh"

//////////////////////////////////////////////////////////////
//
// The Feynman integrand: ME(p) pdf(x,Q) dp for WW scattering
// diagram in which the W+- decay leptonically
//
// This is for the 0-jet case
//
//////////////////////////////////////////////////////////////

class WW : public WWLeptonsBaseIntegrand
{
public:
  WW( Integrator*, int, TransferFunction* );
  WW( Integrator* );

  virtual ~WW() { }
  
  virtual bool initialize();
    
  void setUseStrangePDF( bool flag = true ) { m_doStrange = flag; }
protected:
  virtual double getME( double[][4] ) const;
  virtual void getArray( const WW0jFeynDiagram&, double[][4] ) const;
  virtual void eventScale( double&, double& ) const;

private:
  bool m_doStrange;
};

#endif
