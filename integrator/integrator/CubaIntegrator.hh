//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#ifndef CUBAINTEGRATOR_HH
#define CUBAINTEGRATOR_HH

#include <iostream>
#include <vector>

#include "integrator/Integrator.hh"

namespace CubaIntegrators
{
  class Vegas : public Integrator
  {
  public:
    Vegas();
    virtual ~Vegas() {}

    virtual void doIntegral( double returnVal[], double error[], int* fail,
			     int* neval, double prob[] ) const;

    void setParams( int start, int increase )
    {m_start = start; m_increase = increase;}

  private:
    int m_start;
    int m_increase;
  };


  class Suave : public Integrator
  {
  public:
    Suave();
    virtual ~Suave() {}

    virtual void doIntegral( double returnVal[], double error[], int* fail,
			     int* neval, double prob[] ) const;

    void setNew( int input ) {m_new = input;}
    void setFlatness( double input ) {m_flatness = input;}

  private:
    int m_new;
    double m_flatness;
  };


  class Divonne : public Integrator
  {
  public:
    Divonne();
    virtual ~Divonne() { delete m_points; }

    virtual void doIntegral( double returnVal[], double error[], int* fail,
			     int* neval, double prob[] ) const;

    void setPartitioningRule( int input ) {m_key_par = input;}
    void setIntegrationRule( int input ) {m_key_int = input;}
    void setRefinementRule( int input ) {m_key_ref = input;}
    void setMaxPass( int input ) {m_maxpass = input;}
    void setBorder( double input ) {m_border = input;}
    void setChiSqr( double input ) {m_chisqr = input;}
    void setMinDev( double input ) {m_minDev = input;}
    void setMaxPeaks( int input ) {m_maxPeaks = input;}

    void setPeakFinder( void (*input)( const int*, const double[], int*, double[] ) ) { m_peakFinder = input; }

    void addPoint( std::vector<double> point );
  private:
    int m_key_par;
    int m_key_int;
    int m_key_ref;
    int m_maxpass;
    double m_border;
    double m_chisqr;
    double m_minDev;

    typedef std::vector<double> VecDouble;
    std::vector<VecDouble>* m_points;

    int m_maxPeaks;

    void (*m_peakFinder)( const int*, const double[], int*, double[] );     
  };


  class Cuhre : public Integrator
  {
  public:
    Cuhre();
    virtual ~Cuhre() {}

    virtual void doIntegral( double returnVal[], double error[], int* fail,
			     int* neval, double prob[] ) const;

    void setCubatureRule( int input ) {m_rule = input;}

  private:
    int m_rule;
  };
}

#endif
