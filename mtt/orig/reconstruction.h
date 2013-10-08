#include <complex>
#include <fstream>
#ifndef RECONSTRUCT_H
#define RECONSTRUCT_H
using namespace std;
typedef complex<double> dcmplx;
inline double dot(double* p1,double* p2)
{  
   return p1[0]*p2[0]-p1[1]*p2[1]-p1[2]*p2[2]-p1[3]*p2[3]; 
};
inline dcmplx dot(dcmplx *p1, dcmplx *p2)
{
   return p1[0]*p2[0]-p1[1]*p2[1]-p1[2]*p2[2]-p1[3]*p2[3]; 
};
class event
{
  
  private:

    double ptot[4],pmiss[4];
    double p3o[4],p4o[4],p5o[4],p6o[4];
    int kf3,kf4,kf5,kf6;
   
    double m531,m642,m31,m42,m1,m2;
    dcmplx p1cpl[4][4],p2cpl[4][4],p531cpl[4][4],p642cpl[4][4],
            p31cpl[4][4],p42cpl[4][4],pg1cpl[4][4],pg2cpl[4][4];
    dcmplx totaljac[4];
    void realcut();
    void solvecpl();

  public:    
     event();
     int nevt; 
    void order(int,int,int,int);
    // inputs  
    void setmasses(double m531a,
          double m642a,double m31a,double m42a,double m1a,double m2a);
    void setmomenta(double *p3a,double *p4a,double *p5a,double *p6a,double *p7a);
    void setkf(int kf3a,int kf4a,int kf5a,int kf6a);
    bool solve();
    bool solve_allcombi();
    bool solve_onecombi();
    void print();
    void printcpl();
    void printfile(ofstream &, int);
    bool solveone; 
    //TRUE=pass after finding one solution, FALSE=find all solutions
    
    int debug0;
    //outputs
     double *p3,*p4,*p5,*p6;
    double p1[4][4],p2[4][4],p531[4][4],p642[4][4],p31[4][4],p42[4][4],pg1[4][4],pg2[4][4];
    double m2res[4];
    double mtcl;
    bool solved;
    int nsolutions;    
    
};    
#endif

      
