//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef DY2j_HH
#define DY2j_HH

#include "matrix/DY2jFeynDiagram.hh"
#include "matrix/FeynIntegrand.hh"

class TauTF;
class JetEnergyResolutionTF;

// ////////////////////////////////////////////////////////////
//
// the integrand: ME(p) pdf(x,Q) dp for Z+jj
// diagrams in which the Z/gamma decays leptonically
//
// ////////////////////////////////////////////////////////////

class DY2j : public FeynIntegrand<DY2jFeynDiagram>
{
public:
  
  //! IntegrationStrategy - control integrand configuration, and also dimensions of integration
  enum IntegrationStrategy 
    {
      kNONE = 0, // dummy configuration ... for no-integration testing of ME
      k2D = 2,   // integrate over jet resolutions
      k4D = 4    // for Z->tt use TF's also for ratios of tauon energy to muon/electron energy
    };
  
  DY2j( Integrator* integrator, int strategy = kNONE );
  
  virtual ~DY2j() { delete m_tau_tf; }
  
  virtual bool setKinematics( double[] );
  
  virtual bool initialize( );
  virtual bool finalize( );

  int getIXP()      const   { return IXP; }
  int getIXM()      const   { return IXM; }

  void setJetTF( JetEnergyResolutionTF* tf ) { m_jet_tf = tf; }
  
protected:
  
  virtual double matrixElement()        const;

  virtual void eventScale( double& sa, 
			   double& sb ) const;
  virtual double totalTF( )             const;
  virtual double phaseSpace( )          const;
  
  IntegrationStrategy m_strategy; ///< integration strategy to use

  int IXP;  ///< and for the tau energy fractions 
  int IXM;
  
  int IJETA;
  int IJETB;

  double m_xp, m_xm;
  
  mutable TauTF* m_tau_tf;
  mutable JetEnergyResolutionTF* m_jet_tf;

  void configure( );  
};

#endif
