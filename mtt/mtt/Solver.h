#ifndef SOLVER_H
#define SOLVER_H

/* A translation for J. Gunion's reconstruction code from Fortran to C++ */

// arXiv 0809.4487

class Solver
{
public:      

  // inputs

  double ptot[4],pmiss[4];
  double p3o[4],p4o[4],p5o[4],p6o[4];
  double m3o,m4o,m5o,m6o;
  double *p3,*p4,*p5,*p6;
  double m3,m4,m5,m6,m531,m642,m31,m42,m1,m2;

  bool solveone; 

  int nevt;

  Solver();   

  void setMasses(double m531a,
		 double m642a,
		 double m31a,
		 double m42a,
		 double m1a,
		 double m2a);

  void setMomenta(double *p3a,
		  double *p4a,
		  double *p5a,
		  double *p6a,
		  double *p7a);
  bool solve();
  bool allcombi();
  bool halfcombi();
  void printl();    
  void swap();
  void order(int,int,int,int);
  void set1d();
    
  // outputs

  double p1[4][4],p2[4][4],p531[4][4],p642[4][4],p31[4][4],p42[4][4],pg1[4][4],pg2[4][4];
  double p1n[16], p2n[16];
  double totaljac[4];
  bool solved;
  int nsolutions;
  int n, i, j, row, col;
};    

#endif

      
