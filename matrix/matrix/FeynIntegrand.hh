//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//         
//

#ifndef FEYNINTEGRAND_HH
#define FEYNINTEGRAND_HH

#include "matrix/Constants.hh"
#include "matrix/Utility.hh"
#include "matrix/PDFGrid.hh"

#include "integrator/Integrator.hh"

#include <string>
#include <utility>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <algorithm>

#include <TMath.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>


using std::string;
using std::vector;

typedef vector< vector<double> > dvectorlist;

// ////////////////////////////////////////////////////////////////////////////
//
// interface template for a matrix element (amplitude) integrand
// 
// ////////////////////////////////////////////////////////////////////////////

template<typename FEYNMANDIAGRAM>
class FeynIntegrand
{
public: 
  FeynIntegrand( const string& name,
		 Integrator* integrator, 
		 unsigned int nVars, 
		 unsigned int nLoop );
  
  virtual ~FeynIntegrand() {}

  static const unsigned ITOT = 0; // << the components of the calculation
  static const unsigned IME  = 1;
  static const unsigned IPDF = 2;
  static const unsigned IPS  = 3;
  static const unsigned ITF  = 4;
  static const unsigned IVOL = 5;
  static const unsigned NFAC = 5;
  
  // performs amplitude calculation @ provided integration point
  double operator()( const double parameters[] );
  
  // ////////////////////////////////////////////////
  //
  // API - sub-classes will define
  // these methods, used by the main
  // calculation in FeynIntegrans()
  //
  // ////////////////////////////////////////////////

  // set integration limits for variables (implements physical constraints)
  virtual bool setDynamicLimits( ) { return true; } 

  // set kinematic quantities in event
  virtual bool setKinematics( double parameters[] ) = 0;

  // calculate locations of poles in amplitude
  virtual void getPeaks( dvectorlist& answer,
			 const double bounds[] ) { }

  // ////////////////////////////////////////////////
  //
  // accessors
  //
  // ////////////////////////////////////////////////

  string name()             const { return m_name;        }
  long ncalls()             const { return m_ncalls;      }
  Integrator* integrator()  const { return m_integrator;  }
  unsigned int dimension()  const { return m_nvars;       }
  unsigned int attempts()   const { return m_nattempts;   }
  unsigned int maximum()    const { return m_maxattempts; }
  int ihelicity()           const { return m_ihel;        }
  bool limits( unsigned, double&, double& ) const;

  // prepare for integration
  virtual bool initialize();

  // increment the number of integration attempts
  virtual bool increment();

  // wrap up the integration
  virtual bool finalize();

  // ////////////////////////////////////////////////
  //
  // setters
  //
  // ////////////////////////////////////////////////
  
  void setIntegrator( Integrator* ptr ) { if( ptr != NULL ) m_integrator = ptr; }
  void setIntegrationLimits( unsigned int, double, double ); 
  void setName( const std::string& name ) { m_name = name; }
  void setUseOnlyPDF( bool flag = true ) { m_pdfonly = flag; }
  void setNoPDF( bool flag = true ) { m_nopdf = flag; }
  void setPDF( const std::string& pdf ) { m_pdfset = pdf; }
  void setIHelicity( int ihel = -1 ) { m_ihel = ihel; }

  void setDoClearHelicityCombinations( bool flag = false ) const { m_reset = flag ? 1 : 0; }

protected:
  
  // /////////////////////////////////////////////////
  // 
  // the following methods are common to all
  // ME calculations 
  //
  // /////////////////////////////////////////////////

  // transform integration parameters to phase space (matrix element) variables
  void transformParameters( double par[] ) const;

  // perform quick check on kinematics at phase space point
  virtual bool isPossible( const TLorentzVector& ) const;

  // calculate energies of incoming partons (2->N scattering)
  virtual bool getPartonEnergies( double&, 
				  double&, 
				  const TLorentzVector& ) const;

  // ////////////////////////////////////////////////
  // 
  // these methods depend on the Feynman diagram
  // used in the calculation; they encapsulate 
  // the full amplitude: ME pdf(x,Q) TF dp
  //
  // ////////////////////////////////////////////////

  // calculate the matrix element, ME
  virtual double matrixElement() const = 0;

  // calculate the total transfer function TF, for measured -> parton
  virtual double totalTF() const = 0; 

  // calculate the phase space factor, dp
  virtual double phaseSpace() const { return partons.phaseSpace( ); }

  // calculate the scale at which to sample the parton PDF's
  virtual void eventScale( double&, double& ) const = 0;

  // calculate the PDF factor, pdf(x, Q)
  virtual double PDF( int, int ) const;

  // /////////////////////////////////////////////////

#ifndef __CINT__ // ugh, CINT ... WTF?
  virtual void calculateDiagram( void* (*diagram)(double[][4], double*, int*), double[][4], double&, unsigned& ) const;
#else
  virtual void calculateDiagram( void*, double[][4], double&, unsigned& ) const;
#endif 

  // this setter allows descendants to modify integration strategies (and correspondingly
  // # of integration dimensions) on the fly in the c'tor 
  void setNumIntegrationDimensions( unsigned int ndim );
  
private:
  PDFGrid* m_pdfgrid; ///< the PDF data grid
  mutable Integrator*  m_integrator;  ///< the integrator that will be used for this integrand

  string m_name;  ///< name of amplitude calculation
  string m_pdfset; ///< PDF set to use
  unsigned int m_maxattempts; ///< maximum number of attempts for integration
  unsigned int  m_nattempts;  ///< number of integration attempts for this amplitude
  unsigned int m_nvars;  ///< number of free parameters to integrate over
  float  m_volume;  ///< integration volume
  unsigned int  m_nbounds;  ///< number of integration limits
  long m_ncalls;  ///< # of times the integrand function has been evaluated
  int m_ihel; ///< select a particular helicity configuration; 
  
  // ... bounds type stores the pair {a, b - a} for an interval (a, b)
  typedef std::pair<double, double> bounds_t;
  vector<bounds_t> m_bounds;

  mutable int m_reset;

protected:
  bool m_nointegration; ///< just calculate the matrix element w/o integration 
  bool m_pdfonly; ///< flag to disable calculating |ME|, just calculate PDF
  bool m_nopdf; //< don't use PDF factor
  // //////////////////////////////////////////////////////////////////////////
  // store the factors of the integrand separately; the first element will 
  // be the product of all the factors; the individual diagrams are in 
  // elements [IVOL+1:]
  // //////////////////////////////////////////////////////////////////////////
  mutable std::vector<double> m_integrand_factors; ///< factors that comprise the integrand 

public:
  static FEYNMANDIAGRAM partons;   ///< partons used in evaluation of the integrand
  static FEYNMANDIAGRAM measured;  ///< actual measured objects in the event

  static FeynIntegrand<FEYNMANDIAGRAM>* glbl_integrand;

  // ------------------------- ======= ------------------------- ======= -------------------------

  static void peakFinder( const int* nDim, 
			  const double bounds[],
			  int* nPoints,
			  double peaks[] ) 
  {
    dvectorlist answer;
    glbl_integrand->getPeaks( answer, bounds );
    if( answer.size() > static_cast<unsigned>(*nPoints) )
    {
      throw std::runtime_error("too many peaks returned in peakFinder");
    }
    *nPoints = answer.size();
    unsigned counter = 0;
    for(unsigned int ii = 0; ii < answer.size(); ++ii)
    {
      for(unsigned int jj = 0; jj < answer[ii].size(); ++jj)
      {
	peaks[counter++] = answer[ii][jj];
      }
    }
  }

  // ------------------------- ======= ------------------------- ======= -------------------------
  
  static int fcn_wrapper (const int*, const double params[],
			  const int*, double results[], void *userdata) 
  {
    if( glbl_integrand != NULL )
    {
      results[0] = (*glbl_integrand)( params );
      if( std::isnan(results[0]) )
	results[0]=0;
      return true;
    }
    std::cout << "ERROR integrand pointer is 0x0!" << std::endl;
    return false;
  }
  
  // ------------------------- ======= ------------------------- ======= -------------------------
  
  TStopwatch tPrelim, tSetKine, tQuarks, tME, tPhase, tTF, tPDF, tTotal; ///< timers for benchmarking  
};

#ifndef __CINT__
#include "matrix/FeynIntegrand.icc"
template<typename FEYNMANDIAGRAM>
inline bool FeynIntegrand<FEYNMANDIAGRAM>::limits( unsigned il, double& xl, double& xh ) const
{
  if( il < dimension() )
  {
    xl = m_bounds[il].first;
    xh = xl + m_bounds[il].second;
    return true;
  }
  return false;
}
#endif

#endif
