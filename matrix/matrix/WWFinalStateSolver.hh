#ifndef WWFINALSTATESOLVER_HH
#define WWFINALSTATESOLVER_HH

#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <limits>
#include <vector>
#include <iostream>
#include <cmath>

#include <TLorentzVector.h>

#include "matrix/Constants.hh"

struct KineStore 
{
  TLorentzVector nul; 
  TLorentzVector nur; 
  TLorentzVector lp;  
  TLorentzVector lm;  
  TLorentzVector wp;   
  TLorentzVector wm;   
  TLorentzVector pa;  
  TLorentzVector pb;  
  double weight;

  double phaseSpace() { return 1.0; }
};

template<typename TKIN>
int qwqwqhgy( const TKIN& inputs, double  qX, double  qY,
	      double qh_sqr, double qwp_sqr, double qwm_sqr, double gy,
              std::vector<TKIN>& sol, TLorentzVector branches = TLorentzVector() );

template<typename TKIN>
int qwqwxaxb( const TKIN& inputs, double  qX, double  qY,
	      double qwp_sqr, double qwm_sqr, double xa, double xb,
              std::vector<TKIN>& sol, TLorentzVector branches = TLorentzVector() );

// template<typename TKIN>
// int qwqwpzpz( const TKIN& inputs, double  qX, double  qY,
// 	      double qwp_sqr, double qwm_sqr, double nuZ, double nbZ,
// 	      std::vector<TKIN>& sol, TLorentzVector branches = TLorentzVector() );

#include "matrix/WWFinalStateSolver.icc"

// void quad(double a, double b1, double c, double *sr, double *si, double *lr, double *li);
// void fxshfr(int l2, int *nz);
// void quadit(double *uu, double *vv,int *nz);
// void realit(double sss, int *nz, int *iflag);
// void calcsc(int *type);
// void nextk(int *type);
// void newest(int type, double *uu, double *vv);
// void quadsd(int n, double *u, double *v, double *p, double *q, double *a, double *b);
// int rpoly(double *op, int degree, double *zeror, double *zeroi, int info[] );

#endif
