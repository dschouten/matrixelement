//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef WW1J_HH
#define WW1J_HH

#include "matrix/WW1jLeptonsBaseIntegrand.hh"

//////////////////////////////////////////////////////////////
//
// The Feynman integrand: ME(p) pdf(x,Q) dp for H->WW + jets
// diagram in which the W+- decay leptonically
//
// This is for the 1-jet case, in which gg fusion is dominant
// diagram. Use the HEFT model in Madgraph
//
//////////////////////////////////////////////////////////////

class WW1j : public WW1jLeptonsBaseIntegrand
{
public:
  WW1j( Integrator*, int, TransferFunction* );

  virtual ~WW1j() { }

  virtual bool initialize();

  void setUseStrangePDF( bool flag = true ) { m_doStrange = flag; }

private:
  bool m_doStrange; ///< flag to do strange quark PDF

protected:
  virtual void eventScale( double&, double& ) const;
  virtual double getME( double[][4] ) const;
};

#endif
