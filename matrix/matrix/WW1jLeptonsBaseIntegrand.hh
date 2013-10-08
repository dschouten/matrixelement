//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef WW1JLEPTONSBASEINTEGRAND_HH
#define WW1JLEPTONSBASEINTEGRAND_HH

#include "matrix/WW1jFeynDiagram.hh"
#include "matrix/FeynIntegrand.hh"
#include "matrix/TransferFunction.hh"

#ifndef __CINT__
#include "matrix/WWFinalStateSolver.hh"
#endif

// ////////////////////////////////////////////////////////////
//
// The Feynman integrand: ME(p) pdf(x,Q) dp for WW + jet
// diagrams in which the W+- decay leptonically
//
// This is for the 1-jet case
//
// Calculation modes: 
//
//    NEUTRINO_2D - 2d integration: pz for \f$ \nu\f$  and \f$ \bar{\nu}\f$ ,
//                    assuming (unphysical) case where px
//                    and py for both are known
//
//    NEUTRINO_4D - 4d integration: pT, phi, pz for \f$ \nu\f$  and
//                    pz for \f$ \bar{\nu}\f$ ; the MET px and py are
//                    assumed to be distributed over delta fn's
//
//    NEUTRINO_5D - 5d integration: pT, phi, pz for \f$ \nu\f$  and
//                    pz for \f$ \bar{\nu}\f$ , pT for jet; the MET px and py are
//                    assumed to be distributed over delta fn's, use
//                    transfer function for jet pT
//
// \f$ \bar{\nu}\f$  = A
// \f$ \nu\f$        = B
//
// ////////////////////////////////////////////////////////////

class WW1jLeptonsBaseIntegrand : public FeynIntegrand<WW1jFeynDiagram>
{
public:

  //! IntegrationStrategy - control integrand configuration, and also == dimensions of integration
  enum IntegrationStrategy 
    {
      kNEUTRINO_2D =  2, ///< the transverse components for \f$ \nu\f$  and \f$ \bar{\nu}\f$  are known
      kNEUTRINO_4D =  4, ///< the x,y MET components and jet E are known
      kNEUTRINO_5D =  5, ///< the x,y MET components are known, integrate over jet E resolution
      kWMASS_4D    = -4, ///< MET x,y are \f$ \delta\f$ -functions, integrate W virtuality and longitudinal \f$ \nu\f$  
      kWMASS_5D    = -5, ///< MET x,y are \f$ \delta\f$ -functions, integrate W virtuality, jet E, and longitudinal \f$ \nu\f$  
    };
  
  WW1jLeptonsBaseIntegrand( const std::string& name,
			    Integrator*,
			    TransferFunction* = NULL,
			    IntegrationStrategy strategy = kNEUTRINO_4D );

  virtual ~WW1jLeptonsBaseIntegrand() { }

  virtual bool setDynamicLimits( );
  virtual bool setKinematics( double[] );
  
  virtual bool initialize( );
  virtual bool finalize( );

  int getIPXA() const { return IPXA; }
  int getIPYA() const { return IPYA; }

  int getIPTA() const { return IPTA; }
  int getIPHIA() const { return IPHIA; }

  int getIPZA() const { return IPZA; }
  int getIPZB() const { return IPZB; }

  int getIWP() const { return IWP; }
  int getIWM() const { return IWM; }

  int getIT() const { return IT; }
  int getIY() const { return IY; }

  int getIE() const { return IE; }  
  
  //  ////////////////////////////////////////////////////////////
  //
  // these kinematic functions must be redefined 
  // because they depend on neutrino momenta solutions
  //
  
  virtual bool isPossible( const TLorentzVector& ) const;
  virtual bool getPartonEnergies( double&, double&, const TLorentzVector& ) const;  

protected:

  virtual double matrixElement()        const;
  virtual double phaseSpace()           const;

  virtual void eventScale( double& sa, 
			   double& sb ) const;
  virtual double totalTF( )             const;

  virtual double getME( double[][4] )   const = 0;

protected:
  IntegrationStrategy m_strategy; ///< integration strategy to use

  TransferFunction* m_tf_jet; ///< transfer function for jet  

  mutable bool m_me_flag; ///< flag indicates whether ME calculation is ready

  mutable std::vector<WW1jFeynDiagram> m_kine_solutions; ///< solutions from coupled kinematic constraints

  int IPXA; //< integration variable index for neutrino transverse momentum components
  int IPYA;

  int IPTA; //< integration variable index for neutrino transverse momentum components
  int IPHIA;

  int IPZA; //< integration variable index for neutrino longitudinal momentum components
  int IPZB;
  
  int IWP; //< integration variable index for W+ and W- mass
  int IWM; 

  int IT; //< integration variable index for tau, rapidity
  int IY; 

  int IE; //< integration variable index for jet energy

  void configure( );

  void getArray( const WW1jFeynDiagram&, double[][4] ) const;
};

// ------------------------- ======= ------------------------- ======= -------------------------
inline void WW1jLeptonsBaseIntegrand::configure( )
{
  IPXA = IPYA  = -1;
  IPZA = IPZB  = -1;
  IPTA = IPHIA = -1;
  IT   = IY    = -1;
  IE   = -1;

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
      IE	= 4;
      break;
    case kWMASS_4D :
      IWP	= 0;
      IWM	= 1;
      IT	= 2;
      IY	= 3;
      break;
    case kWMASS_5D :
      IWP	= 0;
      IWM	= 1;
      IT	= 2;
      IY	= 3;
      IE	= 4;
      break;
  }
}
#endif
