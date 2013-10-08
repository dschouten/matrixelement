//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          Bernd Stelzer (stelzer at sfu dot ca)
//          David Hall    (david dot hall at physics dot ox dot ac dot uk)
//          Philip Chang  (philip dot chang at cern dot ch
//

#ifndef WFAKEINTEGRAND_HH
#define WFAKEINTEGRAND_HH

#include "matrix/WFakeFeynDiagram.hh"
#include "matrix/FeynIntegrand.hh"
#include "matrix/TransferFunction.hh"

#include <limits>
#include <string>

class TF1;

// ///////////////////////////////////////////////////////////////////////////////
//
// The Feynman integrand: ME(p) pdf(x,Q) dp for W(->l,nu) + X(->fake electron)
//
// ///////////////////////////////////////////////////////////////////////////////

class WFake : public FeynIntegrand<WFakeFeynDiagram>
{
public:
  
  //! IntegrationStrategy - control integrand configuration, and also == dimensions of integration
  enum IntegrationStrategy 
    {
      k1D = 1, // < integrate over W q**2
      k2D = 2  // < also integrate over fake - lepton momentum
    };
  
  WFake( Integrator*,  
	 IntegrationStrategy strategy = k1D,
	 TransferFunction* faketf = NULL );
  
  virtual ~WFake();
  
  virtual bool setDynamicLimits( );
  virtual bool setKinematics( double[] );
  
  virtual bool initialize( );
  virtual bool finalize( );

  void setWgammaMode( ) { m_wj_mode = false; m_wg_mode = true; }
  void setWjetMode( )   { m_wj_mode = true; m_wg_mode = false; }

  bool isWjetMode( ) const { return m_wj_mode; }
  bool isWgammaMode( ) const { return m_wg_mode; }
  
  int  getIW()    const { return IW;    }
  int  getIFAKE() const { return IFAKE; }
  
  //  ////////////////////////////////////////////////////////////
  //
  // these kinematic functions must be redefined 
  // because they depend on neutrino momenta solutions
  //
  
  virtual bool isPossible( const TLorentzVector& ) const;
  virtual bool getPartonEnergies( double&, double&, const TLorentzVector& ) const;
  
  //  ////////////////////////////////////////////////////////////

  void setFakeResponse( const std::string& f );
  std::string getFakeResponse( ) const;
  
protected:
  mutable bool m_me_flag; ///< flag indicates whether ME calculation is ready

  bool m_wg_mode;
  bool m_wj_mode;
  
  virtual double matrixElement()        const;
  virtual void eventScale( double& sa, 
			   double& sb ) const;
  virtual double totalTF( )             const;
  virtual double phaseSpace( )          const;
  
  double getME( double[][4] )   const;
  
  IntegrationStrategy m_strategy; ///< integration strategy to use
  
  TransferFunction* m_fake_tf; ///< transfer function for fake  

  TF1* m_fake_response; ///< optionally use a simple response function (not a TF!)
  
  int IFAKE;
  int IW;
  
  // solve for z momenta of neutrinos given Q values for W bosons
  void solveZ( double q, unsigned& zsol_n, double zsol_buff[] ) const;

  double m_w_mass;
  double m_fake_param;
  
  void configure( );
  
  void getArray( const WFakeFeynDiagram&, double[][4] ) const;
};

inline void WFake::configure( )
{
  IW    =  0;
  IFAKE = -1;
  if( m_strategy == k2D )
  {
    IFAKE = 1;
  }
}

#endif
