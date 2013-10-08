//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef TT1J_HH
#define TT1J_HH

#include "matrix/FeynIntegrand.hh"
#include "matrix/WW2jFeynDiagram.hh"
#include "matrix/TransferFunction.hh"

#include "ttanalysis/tTSolver.h"
#include "mtt/Solver.h"

// ////////////////////////////////////////////////////////////
//
// The Feynman integrand: ME(p) pdf(x,Q) dp for 
// tT with t->l,nu,b
// 
// This is for the single jet case (either second jet pT < X GeV or 
// outside of detector acceptance)
//
// This is a complicated calculation:
//  - redefine MET given 'missing' jet (provided by integrator)
//  - solve for neutrino momenta
//  - loop over jet <-> top quark combinations
//  
// Also specify low pT cutoff scale for the jet integration
//
// ////////////////////////////////////////////////////////////

class TT1j : public FeynIntegrand<WW2jFeynDiagram>
{
public:
  
  //! IntegrationStrategy - control integrand configuration, and also == dimensions of integration
  enum IntegrationStrategy 
    {
      kJET_3D = 3, ///< integrate over an unobserved b-jet pT, rapidity and azimuth
      kJET_4D = 4, ///< use TF for observed jet E resolution, integrate over unobserved jet
      kJET_WMASS_6D = -6 ///< integrate over W virtualities, use TF for jet E resolution, integrate over missing jet
    }; 
  
  TT1j( Integrator*, int strategy );

  virtual ~TT1j( ) { }

  //  ////////////////////////////////////////////////////////////
  //
  // these kinematic functions must be redefined 
  // because they depend on neutrino momenta solutions
  //

  virtual bool isPossible( const TLorentzVector& ) const;
  virtual bool getPartonEnergies( double&, double&, const TLorentzVector& ) const;

  //  ////////////////////////////////////////////////////////////
  
  virtual bool setDynamicLimits( );
  virtual bool setKinematics( double[] );

  virtual bool initialize( );
  virtual bool finalize( );

  int getIPTB( )  const	{ return IPTB;	} 
  int getIPHIB( ) const { return IPHIB; }
  int getIETAB( ) const { return IETAB; }
  int getIEA( )   const	{ return IEA;	}
  int getIWP( )   const { return IWP; }
  int getIWM( )   const { return IWM; }

  void setMass( double mass );

  void setCutoffScale( double scale );

  void setTransferFunctions( TransferFunction* met_tf,
			     TransferFunction* inv_jet_tf,
			     TransferFunction* obs_jet_tf )
  {
    m_met_tf		= met_tf	;
    m_inv_jet_tf	= inv_jet_tf	;
    m_obs_jet_tf	= obs_jet_tf	;
  }
  void setMissingETTF( TransferFunction* met_tf )		{ m_met_tf	= met_tf;	}
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

  int IPTB; ///< integration variable index for jet momenta
  int IPHIB;
  int IETAB;
  int IEA; ///< integration variable index for jet pT
  int IWP; ///< integration variable index for W virtualities
  int IWM; 

  IntegrationStrategy m_strategy; ///< integration strategy to use

  mutable TransferFunction* m_obs_jet_tf; ///< transfer function for observed jet
  mutable TransferFunction* m_inv_jet_tf; ///< transfer function for invisible jet
  mutable TransferFunction* m_met_tf; ///< transfer function for MET

  mutable tTSolver m_solver; ///< use Sonennschein, Phys. Rev. D 73 054015 (2006)

  mutable unsigned m_nsolutions;

  bool solveTTBar( unsigned count = 0 ) const;

  double m_mass; ///< top quark mass

  double m_cutoffScale; ///< low pT cutoff for integration of invisible jet

private:
  mutable std::vector<double> m_pnux;
  mutable std::vector<double> m_pnuy;
  mutable std::vector<double> m_pnuz;
  mutable std::vector<double> m_pnubx;
  mutable std::vector<double> m_pnuby;
  mutable std::vector<double> m_pnubz; 
  mutable std::vector<double> m_cd_diff;

  mutable bool m_me_flag; ///< set this flag while performing ME calculation

  double m_wp_virt;
  double m_wm_virt;

  mutable int m_icomb;
};

// ------------------------- ======= ------------------------- ======= -------------------------
inline void TT1j::configure( )
{

  std::cout << "\t mT = " << m_mass << " width = " << hepstd::topDecayWidth( m_mass ) << std::endl;
  std::cout << "\t integration strategy: " << m_strategy << std::endl;
  
  IPTB = IETAB = IPHIB = -1;
  IEA = -1;

  switch( m_strategy )
  {
    case kJET_3D:
      IPTB  = 0;
      IETAB = 1;
      IPHIB = 2;
      break;
    case kJET_4D:
      IPTB  = 0;
      IETAB = 1;
      IPHIB = 2;
      IEA  = 3;
      break;
    case kJET_WMASS_6D:
      IPTB  = 0;
      IETAB = 1;
      IPHIB = 2;
      IEA  = 3;
      IWP  = 4;
      IWM  = 5;
      break;
  }
}

// ------------------------- ======= ------------------------- ======= -------------------------
inline void TT1j::setMass( double mass )
{
  m_mass = mass;
}

// ------------------------- ======= ------------------------- ======= -------------------------
inline void TT1j::setCutoffScale( double scale )
{
  m_cutoffScale = scale;
}

#endif
