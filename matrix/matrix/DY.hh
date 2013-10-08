//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//          
//

#ifndef DY_HH
#define DY_HH

#include "matrix/DYFeynDiagram.hh"
#include "matrix/FeynIntegrand.hh"
#include "matrix/TransferFunction.hh"
#include "matrix/TauTF.hh"

// ////////////////////////////////////////////////////////////
//
// the integrand: ME(p) pdf(x,Q) dp for Z+j
// diagrams in which the Z/gamma decays leptonically
//
// ////////////////////////////////////////////////////////////

class DY : public FeynIntegrand<DYFeynDiagram>
{
public:
  
  //! IntegrationStrategy - control integrand configuration, and also dimensions of integration
  enum IntegrationStrategy 
    {
      kNONE = 0,      // dummy configuration ... for no-integration testing of ME
      kJET_1D = 1,    // integrate over rapidity of the jet
      kJET_TT_3D = 3  // for Z->tt use TF's also for ratios of tauon energy to muon/electron energy
    };
  
  DY( Integrator* integrator, int strategy = kNONE );
  
  virtual ~DY() { }
  
  virtual bool setKinematics( double[] );
  
  virtual bool initialize( );
  virtual bool finalize( );

  int getIETA()     const   { return IETA; }
  int getIXP()      const   { return IXP; }
  int getIXM()      const   { return IXM; }
  
protected:
  
  virtual double matrixElement()        const;

  virtual void eventScale( double& sa, 
			   double& sb ) const;
  virtual double totalTF( )             const;
  virtual double phaseSpace( )          const;
  
  IntegrationStrategy m_strategy; ///< integration strategy to use
  
  int IETA; ///< index of the integration variable for jet rapidity
  int IXP;  ///< and for the tau energy fractions 
  int IXM;
  
  double m_xp, m_xm;
  
  mutable TauTF m_tauTF;

  void configure( );  
};

#endif
