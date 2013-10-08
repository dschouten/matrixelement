//
// author: Doug Schouten <doug dot schouten at triumf dot ca>
//

#ifndef HWW_HH
#define HWW_HH

#include "matrix/WWLeptonsBaseIntegrand.hh"
#include "matrix/TransferFunction.hh"

//////////////////////////////////////////////////////////////
//
// the integrand: ME(p) pdf(x,Q) dp for H->WW
// diagram in which the W+- decay leptonically
//
//////////////////////////////////////////////////////////////

class HWW : public WWLeptonsBaseIntegrand
{
public:

  HWW( Integrator*, double, int, TransferFunction*, bool useNarrowWidth = false );

  virtual ~HWW() { }

  virtual void getPeaks( dvectorlist& answer, 
			 const double bounds[] );

  virtual bool setKinematics( double[] );

  virtual double mass( )  const { return m_mass;  }
  virtual double width( ) const { return m_width; }

  void setMass( double mass );
  void setWidth( double width );

  void setUseSM( bool flag = true ) { m_useSM = flag; if( flag ) { m_useHEFT = !flag; m_useRS = !flag; m_useJHU = !flag; } }
  bool useSM( ) const { return m_useSM; }

  void setUseRS( bool flag = true ) { m_useRS = flag; if( flag ) { m_useSM = !flag; m_useHEFT = !flag; m_useJHU = !flag; } }
  bool useRS( ) const { return m_useRS; }

  void setUseHEFT( bool flag = true ) { m_useHEFT = flag; if( flag ) { m_useSM = !flag; m_useRS = !flag; m_useJHU = !flag; } }
  bool useHEFT( ) const { return m_useHEFT; }

  void setUseJHU( bool flag = true, int J = 2, int CP = +1 );
  bool useJHU( ) const { return m_useJHU; }

  bool useNarrowWidth( ) const { return m_useNW; } 

  void setUseBoxPS( bool flag = false ) { m_useBoxPS = flag; }
  bool useBoxPS( ) const { return m_useBoxPS; }

  int getIH()	const { return IH; }
  int getIY()	const { return IY; }

  void setJHUFractionQQ( double frac = 0.0 ) { m_jhufracqq = frac <= 1 ? frac : 0.0; }

protected:
  virtual double getME( double[][4] ) const;
  virtual void getArray( const WW0jFeynDiagram&, double[][4] ) const;

  void reconfigure( );

  int IY; //< rapidity of H 
  int IH; //< ... and Q^{2} 
  
private:
  double m_mass;
  double m_width;

  int m_jhuspin; //< J=0,2 for JHU ME
  int m_jhucp;   //< CP=-1,+1 for JHU ME
  double m_jhufracqq; //< fraction of qqbar production for J=2 JHU ME

  bool m_useSM; //< flag to use SM bb > H kludge instead of HEFT
  bool m_useRS; //< flag to use RS p gg > y instead of HEFT
  bool m_useHEFT; //< flag to use HEFT model 
  bool m_useJHU; //< flag to use JHU S=0,2 
  bool m_useNW; //< flag to use narrow width approximation for H mass integral
  bool m_useBoxPS; //< flag to flatten out BW by transformation q -> tan(t)
};

// ------------------------- ======= ------------------------- ======= -------------------------
inline void HWW::setUseJHU( bool flag, int J, int CP ) 
{ 
  m_useJHU = flag;
  m_jhuspin = J; 
  m_jhucp = CP; 
  if( flag ) 
  { 
    m_useSM = !flag; 
    m_useRS = !flag; 
    m_useHEFT = !flag; 
  } 
}

#endif
