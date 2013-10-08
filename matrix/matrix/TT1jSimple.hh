//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef TT1JSIMPLE_HH
#define TT1JSIMPLE_HH

#include "matrix/FeynIntegrand.hh"
#include "matrix/WW2jFeynDiagram.hh"
#include "matrix/TransferFunction.hh"

// ////////////////////////////////////////////////////////////
//
// The Feynman integrand: ME(p) pdf(x,Q) dp for 
// tT with t->l,nu,b
// 
// This is for the single jet case (either second jet pT < X GeV or 
// outside of detector acceptance)
//
// ////////////////////////////////////////////////////////////

class TT1jSimple : public FeynIntegrand<WW2jFeynDiagram>
{
public:
  
  //! IntegrationStrategy - control integrand configuration, and also == dimensions of integration
  enum IntegrationStrategy 
    {
      kJET_7D = 7, ///< integrate over missing jet pT,\eta,\phi, px,py,pz for \nu and pz for \bar{\nu}
      kJET_8D = 8 ///< integrate over missing jet pT,\eta,\phi, observed jet E, px,py,pz for \nu and pz for \bar{\nu}
    };
    
  TT1jSimple( Integrator*, int strategy );
  
  TT1jSimple( Integrator*, int strategy, double mass );

  virtual ~TT1jSimple( ) { }

  //  ////////////////////////////////////////////////////////////
  
  virtual bool setDynamicLimits( );
  virtual bool setKinematics( double[] );

  virtual bool initialize( );
  virtual bool finalize( );

  int getIPTJB( )  const { return IPTJB;  } 
  int getIPHIJB( ) const { return IPHIJB; }
  int getIETAJB( ) const { return IETAJB; }
  int getIEJA( )   const { return IEJA;	  }
  int getIPTA( )   const { return IPTA;   }
  int getIPHIA( )  const { return IPHIA;  }
  int getIPZA( )   const { return IPZA;   }
  int getIPZB( )   const { return IPZB;   }

  void setMass( double mass );

  void setCutoffScale( double scale );
  
  void setGluonOnly( bool flag = true ) { m_gg_only = flag; }

  void setTransferFunctions( TransferFunction* inv_jet_tf,
			     TransferFunction* obs_jet_tf )
  {
    m_inv_jet_tf	= inv_jet_tf;
    m_obs_jet_tf	= obs_jet_tf;
  }

  void setObservedJetTF( TransferFunction* obs_jet_tf )		{ m_obs_jet_tf	= obs_jet_tf;	}
  void setInvisibleJetTF( TransferFunction* inv_jet_tf )	{ m_inv_jet_tf	= inv_jet_tf;	}

protected:
  virtual double phaseSpace() const;

  virtual void eventScale( double& sa, 
			   double& sb ) const;
  virtual double totalTF()              const;

  virtual double matrixElement()  const;

private:

  double getME( unsigned& idiagram ) const;

  void configure( );

  int IPTJB; ///< integration variable index for jet momenta
  int IPHIJB;
  int IETAJB;
  int IEJA; ///< integration variable index for jet pT
  int IPTA; ///< integration variable index for \nu x,y,z
  int IPHIA;
  int IPZA;
  int IPZB; 

  IntegrationStrategy m_strategy; ///< integration strategy to use

  mutable TransferFunction* m_obs_jet_tf; ///< transfer function for observed jet
  mutable TransferFunction* m_inv_jet_tf; ///< transfer function for invisible jet

  double m_mass; ///< top quark mass

  double m_cutoffScale; ///< low pT cutoff for integration of invisible jet

  bool m_gg_only;
};

// ------------------------- ======= ------------------------- ======= -------------------------
inline void TT1jSimple::configure( )
{

  std::cout << "\t mT = " << m_mass << " width = " << hepstd::topDecayWidth( m_mass ) << std::endl;
  std::cout << "\t integration strategy: " << m_strategy << std::endl;
  
  IPTJB = IETAJB = IPHIJB = -1;
  IEJA  = -1;
  IPTA = IPHIA = IPZA  = -1;
  IPZB = -1;

  switch( m_strategy )
  {
    case kJET_7D:
      IPTJB = 0;
      IETAJB = 1;
      IPHIJB = 2;
      IPTA = 3;
      IPHIA = 4;
      IPZA = 5;
      IPZB = 6;
      break;
    case kJET_8D:
      IPTJB = 0;
      IETAJB = 1;
      IPHIJB = 2;
      IEJA = 3;
      IPTA = 4;
      IPHIA = 5;
      IPZA = 6;
      IPZB = 7;
      break;
  }
}

// ------------------------- ======= ------------------------- ======= -------------------------
inline void TT1jSimple::setMass( double mass )
{
  m_mass = mass;
}

// ------------------------- ======= ------------------------- ======= -------------------------
inline void TT1jSimple::setCutoffScale( double scale )
{
  m_cutoffScale = scale;
}

#endif
