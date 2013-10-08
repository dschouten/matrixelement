//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef WWLEPTONSBASEINTEGRAND_HH
#define WWLEPTONSBASEINTEGRAND_HH

#include "matrix/WW0jFeynDiagram.hh"
#include "matrix/FeynIntegrand.hh"
#include "matrix/TransferFunction.hh"

#include <limits>
#include <fstream>

#ifndef __CINT__
#include "matrix/WWFinalStateSolver.hh"
#endif

// ////////////////////////////////////////////////////////////
//
// The Feynman integrand: ME(p) pdf(x,Q) dp for WW
// diagrams in which the W+- decay leptonically
//
// This is for the 0-jet case
//
// Calculation modes: 
//
//    NEUTRINO_2D - 2d integration: pz for \nu and \bar{\nu},
//                    assuming (unphysical) case where px
//                    and py for both are known
//
//    NEUTRINO_4D - 4d integration: pT, phi, pz for \nu and
//                    pz for \bar{\nu}; the MET px and py are
//                    assumed to be distributed over delta fn's
//
//    NEUTRINO_6D - 6d integration over px, py and pz for \nu
//                    and pz for \bar{\nu}, and using TF 
//                    for p_{T}^{H}, \phi^{H}
//    WMASS_4/6D  - same, except two degrees of freedom are 
//                   set using W mass integration
//
//
// NOTE the W-mass coordinates integration is not yet validated
//     
//
// \nu = B
// \bar{\nu} = A
//
// ////////////////////////////////////////////////////////////

class WWLeptonsBaseIntegrand : public FeynIntegrand<WW0jFeynDiagram>
{
public:
  
  //! IntegrationStrategy - control integrand configuration, and also == dimensions of integration
  enum IntegrationStrategy 
    {
      kNEUTRINO_2D =  2,
      kNEUTRINO_4D =  4,
      kNEUTRINO_5D =  5, // include integration over recoil pT
      kNEUTRINO_6D =  6, // include integration over recoil pT, phi
      kWMASS_4D    = -4,
      kWMASS_5D    = -5, // include integration over recoil pT
      kWMASS_6D    = -6  // include integration over recoil pT, phi
    };
  
  WWLeptonsBaseIntegrand( const std::string& name,
			  Integrator*,
			  TransferFunction* = NULL, 
			  IntegrationStrategy strategy = kNEUTRINO_4D );
  
  virtual ~WWLeptonsBaseIntegrand() { }
  
  virtual bool setDynamicLimits( );
  virtual bool setKinematics( double[] );
  
  virtual double mass( )  const { return -1; } 
  virtual double width( ) const { return -1; }
  
  virtual bool initialize( );
  virtual bool finalize( );

  virtual void applyRecoilBoost( WW0jFeynDiagram&, int sign = 1 );
  
  void setStorePhaseSpace( bool flag = true ) { m_phase_space_flag = flag; }

  TransferFunction* getTF( ) { return m_tf; }  

  int getIXA()          const   { return IXA;    }
  int getIXB()          const   { return IXB;    }

  int getIT()           const   { return IT;     }
  int getIY()           const   { return IY;     }

  int getIWP()		const	{ return IWP;	 }
  int getIWM()		const	{ return IWM;	 }

  int getISPT()	        const	{ return ISPT;   } 
  int getISDPHI()	const	{ return ISDPHI; } 

  int getIPTA()         const   { return IPTA;   }
  int getIPHIA()        const   { return IPHIA;  } 
  int getIPZA()		const	{ return IPZA;	 }
  int getIPZB()		const	{ return IPZB;	 }

  
  //  ////////////////////////////////////////////////////////////
  //
  // these kinematic functions must be redefined 
  // because they depend on neutrino momenta solutions
  //
  
  virtual bool isPossible( const TLorentzVector& ) const;
  virtual bool getPartonEnergies( double&, double&, const TLorentzVector& ) const;
  
  //  ////////////////////////////////////////////////////////////
  
protected:
  mutable bool m_me_ready_flag; ///< flag indicates whether ME calculation is ready

  virtual double matrixElement()        const;
  virtual void eventScale( double& sa, 
			   double& sb ) const;
  virtual double totalTF( )             const;
  virtual double phaseSpace( )          const;

  virtual double getME( double[][4] )   const = 0;
  
  IntegrationStrategy m_strategy; ///< integration strategy to use
  
  TransferFunction* m_tf; ///< transfer function for MET
  
  int IPZA; ///< integration variable index for neutrino momentum components
  int IPZB;
  int IPTA; 
  int IPHIA;
  
  int IWP; ///< W+ and W- 
  int IWM; 
  
  int ISPT; ///< system transverse recoil p_{T} and \phi
  int ISDPHI;

  int ISX; ///< system boost x,y parameters
  int ISY;

  int IXA; ///< x fraction for parton A
  int IXB;

  int IT; ///< change vars xa,xb -> tau = xa*xb, y = 1/2 log( xa / xb )
  int IY; 
  
  virtual void configure( );

  virtual void getArray( const WW0jFeynDiagram&, double[][4] ) const = 0;

  mutable long m_num_ps_nonnull; ///< allow to count phase space points with ME != 0

private:
  unsigned m_ps_counter; ///< count phase space integrations, not really trustworthy (based on # of initialize() calls)
  bool m_phase_space_flag; ///< flag to dump out the full phase space sampled in the integration

  void bookPhaseSpaceNtuple( );
  void dumpPhaseSpacePoint( const WW0jFeynDiagram&, double ) const;

protected:  
  std::vector<WW0jFeynDiagram> m_kine_solutions; ///< solutions from coupled kinematic constraints
  
};

inline void WWLeptonsBaseIntegrand::configure( )
{
  /////////////////////////////
  // \bar{\nu} = A
  // \nu       = B
  /////////////////////////////
  
  IPZA = IPZB   = -1;
  ISPT = ISDPHI = -1;
  IPTA = IPHIA  = -1;
  IWP  = IWM    = -1;
  IXA  = IXB    = -1;
  IT   = IY     = -1;
  
  switch( m_strategy ) {
    case kNEUTRINO_2D :
      IPZA	= 0;
      IPZB	= 1;
      break;
    case kNEUTRINO_4D :
      IPTA	= 0;
      IPHIA	= 1;
      IPZA	= 2;
      IPZB	= 3;
      break;
    case kNEUTRINO_5D :
      IPTA	= 0;
      IPHIA	= 1;
      IPZA	= 2;
      IPZB	= 3;
      ISPT	= 4;
      break;
    case kNEUTRINO_6D :
      IPTA	= 0;
      IPHIA	= 1;
      IPZA	= 2;
      IPZB	= 3;
      ISPT	= 4;
      ISDPHI	= 5;
      break;
    case kWMASS_4D :
      IWP	= 0;
      IWM	= 1;
      IT        = 2;
      IY        = 3;
      break;
    case kWMASS_5D :
      IWP	= 0;
      IWM	= 1;
      IT        = 2;
      IY        = 3;
      ISPT	= 4;
      break;
    case kWMASS_6D :
      IWP	= 0;
      IWM	= 1;
      IT        = 2;
      IY        = 3;
      ISPT	= 4;
      ISDPHI	= 5;
      break;
  }
}
  
#endif
