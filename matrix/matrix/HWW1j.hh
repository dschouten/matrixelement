//
// author: Doug Schouten <doug dot schouten at triumf dot ca>
//

#ifndef HWW1J_HH
#define HWW1J_HH

#include "matrix/WW1jLeptonsBaseIntegrand.hh"

//////////////////////////////////////////////////////////////
//
// The Feynman integrand: ME(p) pdf(x,Q) dp for gH->WW + j
// diagram in which the W+- decay leptonically
//
//////////////////////////////////////////////////////////////

class HWW1j : public WW1jLeptonsBaseIntegrand
{
public:
  HWW1j( Integrator*, double, int, TransferFunction*, bool );

  virtual ~HWW1j() { }

  virtual bool initialize();

  virtual bool setKinematics( double[] );

  double mass( )  const { return m_mass;  }
  double width( ) const { return m_width; }

  void setMass( double mass );
  void setWidth( double width );

  int getIH()		const	{ return IH;	}
  int getIY()		const	{ return IY;	}

  void setUseStrangePDF( bool flag = true ) { m_doStrange = flag; }

protected:
  virtual bool isPossible( const TLorentzVector& t ) const;
  virtual double getME( double[][4] ) const;
  virtual void getArray( const WW1jFeynDiagram&, double[][4] ) const;

  void reconfigure( );

  int IY; //< rapidity of H 
  int IH; //< ... and Q^{2} 

private:
  double m_mass; ///< mass of the Higgs
  double m_width;

  bool m_doStrange;
  bool m_useNW;
};

#endif
