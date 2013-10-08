/* A translation for J. Gunion's reconstruction code from Fortran to C++ */

#include "mtt/Solver.h"

#include <iostream>
#include <math.h>

// #define M_PI  3.1415926535897932

#define M_1_3 0.333333333333333
#define M_EPS 1.0e-29

using namespace std;

inline double dot( double* p1,
		   double* p2 )
{  
  return ( p1[0]*p2[0] - 
	   p1[1]*p2[1] - 
	   p1[2]*p2[2] - 
	   p1[3]*p2[3] );
};

void cubic( double *a, double *x, int &l )
{
  double u[3];
  double w,p,q,dis;  
  if (a[3]!=0.) 
  {     
    w = a[2]/a[3]*M_1_3;
    p = a[1]/a[3]*M_1_3-w*w;
    p = p*p*p;
    q = -0.5*(2.*w*w*w-(a[1]*w-a[0])/a[3]);
    dis = q*q+p;	
    if (dis<0.)
    {   
      double phi = acos(min(1.,max(-1.,q/sqrt(-p))));
      p = 2.*sqrt(cbrt(-p));
      //       P = 2.D+0*(-P)**(5.D-1*M_1_3)
      for(int i = 1;i<4;i++) u[i-1] = p*cos((phi+2.*i*M_PI)*M_1_3)-w;
      x[0] = min(u[0],min(u[1],u[2]));
      x[1] = max(min(u[0],u[1]),max(min(u[0],u[2]),min(u[1],u[2])));
      x[2] = max(u[0],min(u[1],u[2]));
      l = 3;
    }
    else
    {
	
      //         only one real solution!
      dis = sqrt(dis);
      x[0] = cbrt(q+dis)+cbrt(q-dis)-w;
      l = 1;
    }
  }  
  else if (a[2]!=0.)
  {
    //        quadratic problem
    p = 0.5*a[1]/a[2];
    dis = p*p-a[0]/a[2];
    if (dis>=0.) 
    {
      x[0] = -p-sqrt(dis); 
      x[1] = -p+sqrt(dis);
      l = 2;
    }
    else l = 0;
  }
  else if (a[1]!=0.) 
  {

    //       linear equation
    x[0] = -a[0]/a[1];
    l = 1;

  }
  else l = 0;
  for (int i = 0;i<l;i++)
    x[i] = x[i]-(a[0]+x[i]*(a[1]+x[i]*(a[2]+x[i]*a[3])))/(a[1]+x[i]*(2.*a[2]+x[i]*3.*a[3]));
};

void quartic( double *dd, double *sol, double *soli, int & nsol )
{
  double aa[4],z[3];
  int i;
  nsol = 0;
  double a,b,c,d,e;
  a = dd[4];
  b = dd[3];
  c = dd[2];
  d = dd[1];
  e = dd[0];
  if (fabs(dd[4])<M_EPS)
  { 
    cout<<"ERROR: NOT A QUARTIC EQUATION";
    return;
  };   
  double p,q,r;
  p = (-3.*b*b + 8.*a*c)/(8.*a*a);
  q = (b*b*b - 4.*a*b*c + 8.*d*a*a)/(8.*a*a*a);
  r = (-3.*b*b*b*b + 16.*a*b*b*c - 64.*a*a*b*d + 
       256.*a*a*a*e)/(256.*a*a*a*a);

  //       solve cubic resolvent
  aa[3] =  8.;
  aa[2] = -4.*p; 
  aa[1] = -8.*r;
  aa[0] =  4.*p*r - q*q;
  int ncube;
  cubic(aa,z,ncube);
  double zsol;      
  zsol = -1.e+99;
  for(i=0;i<ncube;i++)zsol = max(zsol,z[i]);
  z[0] = zsol;
  double xk2,xk;
  xk2 = 2. * z[0] - p;
  xk  = sqrt(xk2);
  //-----------------------------------------------
  double xl2,xl; 
  if (fabs(xk)<M_EPS) 
  {
    xl2 = z[0]*z[0] - r;
    if (xl2<0.) 
    {
      return;
    }  
    xl  = sqrt(xl2);
  }
  else
  {
    xl = q/(2. * xk);
  };
      
  //-----------------------------------------------
  double sqp,sqm; 
  sqp = xk2 - 4.*(z[0] + xl);
  sqm = xk2 - 4.*(z[0] - xl);

  for(i=0;i<4;i++)soli[i]=0.;

  if (sqp>=0. && sqm>=0.) 
  {
    sol[0] = 0.5*( xk + sqrt(sqp));
    sol[1] = 0.5*( xk - sqrt(sqp));
    sol[2] = 0.5*(-xk + sqrt(sqm));
    sol[3] = 0.5*(-xk - sqrt(sqm));
    nsol = 4;
  }
  else if(sqp>=0.&&sqm<0.) 
  {
    sol[0] =  0.5*(xk + sqrt(sqp));
    sol[1] =  0.5*(xk - sqrt(sqp));
    sol[2] = -0.5*xk; 
    sol[3] = -0.5*xk; 
    soli[2] =  sqrt(-0.25 * sqm);
    soli[3] = -sqrt(-0.25 * sqm);
    nsol = 2;
  }
  else if  (sqp<0.&& sqm>=0.) 
  {
    sol[0] = 0.5*(-xk + sqrt(sqm));
    sol[1] = 0.5*(-xk - sqrt(sqm));
    sol[2] =  0.5*xk; 
    sol[3] =  0.5*xk; 
    soli[2] =  sqrt(-0.25 * sqp);
    soli[3] = -sqrt(-0.25 * sqp);
    nsol = 2;
  }
  else if  (sqp<0.&& sqm<0.) 
  {
    sol[0] = -0.5*xk; 
    sol[1] = -0.5*xk; 
    soli[0] =  sqrt(-0.25 * sqm);
    soli[1] = -sqrt(-0.25 * sqm);
    sol[2] =  0.5*xk; 
    sol[3] =  0.5*xk; 
    soli[2] =  sqrt(-0.25 * sqp);
    soli[3] = -sqrt(-0.25 * sqp);
    nsol = 0;
  } 
  for(i=0;i<4;i++)sol[i] = sol[i] - b/(4.*a);
};  

Solver::Solver()
{
  solved = false;
  solveone = false;
  nsolutions = 0;
};

bool Solver::solve()   
{
  double dd[5],sol[4],soli[4];
  double m3sq,m4sq,m5sq,m6sq,m531sq,m642sq,m31sq,m42sq,m1sq,m2sq;
  int i;
  solved=true;
  nsolutions=0;
  m3sq=m3*m3;
  m4sq=m4*m4;
  m5sq=m5*m5;
  m6sq=m6*m6;
  m531sq=m531*m531;
  m642sq=m642*m642;
  m31sq=m31*m31;
  m42sq=m42*m42;
  m1sq=m1*m1;
  m2sq=m2*m2;
   
  double pvisx,pvisy,pvisz,evis;
  pvisz = 0;
  pvisx=-pmiss[1];
  pvisy=-pmiss[2];
  // pvisz=p3[3]+p4[3]+p5[3]+p6[3];
  // pvis=p3[0]+p4[0]+p5[0]+p6[0];
  double p3dp5,p5dp3,p4dp6,p6dp4;
  p3dp5=dot(p3,p5);
  p5dp3=p3dp5;
  p4dp6=dot(p4,p6);
  p6dp4=p4dp6;
   
  double del2b,del3b,del31,del531;
  del2b=m3sq+m1sq-m31sq+m42sq-m2sq-m4sq;
  del3b=m31sq+m5sq-m531sq+m642sq-m42sq-m6sq+2.*p3dp5-2.*p4dp6;
  del31=m31sq-m1sq-m3sq;
  del531=m531sq-m31sq-m5sq-2.0*p3dp5;
   
  double e3,p3x,p3y,p3z,e4,
    p4x,p4y,p4z,e5,
    p5x,p5y,p5z,e6,
    p6x,p6y,p6z;
  e3=p3[0];
  p3x=p3[1];
  p3y=p3[2];
  p3z=p3[3];

  e4=p4[0];
  p4x=p4[1];
  p4y=p4[2];
  p4z=p4[3];

  e5=p5[0];
  p5x=p5[1];
  p5y=p5[2];
  p5z=p5[3];

  e6=p6[0];
  p6x=p6[1];
  p6y=p6[2];
  p6z=p6[3];

  double detval,cxe1,cxe2,
    cye1,cye2,cze1,
    cze2,czte1,czte2,
    cx,cy,cz,czt;

  detval=
    -4*p3z*p4z*p5y*p6x + 4*p3y*p4z*p5z*p6x + 
    4*p3z*p4z*p5x*p6y - 4*p3x*p4z*p5z*p6y - 
    4*p3z*p4y*p5x*p6z + 4*p3z*p4x*p5y*p6z - 
    4*p3y*p4x*p5z*p6z + 4*p3x*p4y*p5z*p6z;
  cxe1=(-((-4*e5*p3z*p4z*p6y + 4*e3*p4z*p5z*p6y + 
	   4*e5*p3z*p4y*p6z - 4*e3*p4y*p5z*p6z)/
          (-4*p3z*p4z*p5y*p6x + 4*p3y*p4z*p5z*p6x + 
	   4*p3z*p4z*p5x*p6y - 4*p3x*p4z*p5z*p6y - 
	   4*p3z*p4y*p5x*p6z + 4*p3z*p4x*p5y*p6z - 
	   4*p3y*p4x*p5z*p6z + 4*p3x*p4y*p5z*p6z)));
  cxe2=(-((-4*e6*p3z*p4z*p5y + 4*e6*p3y*p4z*p5z + 
	   4*e4*p3z*p5y*p6z - 4*e4*p3y*p5z*p6z)/
          (-4*p3z*p4z*p5y*p6x + 4*p3y*p4z*p5z*p6x + 
	   4*p3z*p4z*p5x*p6y - 4*p3x*p4z*p5z*p6y - 
	   4*p3z*p4y*p5x*p6z + 4*p3z*p4x*p5y*p6z - 
	   4*p3y*p4x*p5z*p6z + 4*p3x*p4y*p5z*p6z)));
  cye1=(((e5*p3z - e3*p5z)*(p4z*p6x - p4x*p6z))/
        (p5z*(-(p3y*p4z*p6x) + p3x*p4z*p6y + p3y*p4x*p6z - 
	      p3x*p4y*p6z) + 
	 p3z*(p4z*p5y*p6x - p4z*p5x*p6y + p4y*p5x*p6z - 
	      p4x*p5y*p6z)));
  cye2=(-(((p3z*p5x - p3x*p5z)*(e6*p4z - e4*p6z))/
          (p5z*(p3y*p4z*p6x - p3x*p4z*p6y - p3y*p4x*p6z + 
		p3x*p4y*p6z) + 
	   p3z*(-(p4z*p5y*p6x) + p4z*p5x*p6y - p4y*p5x*p6z + 
		p4x*p5y*p6z))));
  cze1=((e5*(p3y*p4z*p6x - p3x*p4z*p6y - p3y*p4x*p6z + 
             p3x*p4y*p6z) + 
	 e3*(-(p4z*p5y*p6x) + p4z*p5x*p6y - p4y*p5x*p6z + 
             p4x*p5y*p6z))/
        (p5z*(p3y*p4z*p6x - p3x*p4z*p6y - p3y*p4x*p6z + 
	      p3x*p4y*p6z) + 
	 p3z*(-(p4z*p5y*p6x) + p4z*p5x*p6y - p4y*p5x*p6z + 
	      p4x*p5y*p6z)));
  cze2=(((p3y*p5x - p3x*p5y)*(e6*p4z - e4*p6z))/
        (p5z*(p3y*p4z*p6x - p3x*p4z*p6y - p3y*p4x*p6z + 
	      p3x*p4y*p6z) + 
	 p3z*(-(p4z*p5y*p6x) + p4z*p5x*p6y - p4y*p5x*p6z + 
	      p4x*p5y*p6z)));
  czte1=((e5*(-(p3z*p4y*p6x) + p3y*p4z*p6x + p3z*p4x*p6y - 
	      p3x*p4z*p6y - p3y*p4x*p6z + p3x*p4y*p6z) + 
          e3*(-(p4z*p5y*p6x) + p4y*p5z*p6x + p4z*p5x*p6y - 
	      p4x*p5z*p6y - p4y*p5x*p6z + p4x*p5y*p6z))/
	 (p5z*(p3y*p4z*p6x - p3x*p4z*p6y - p3y*p4x*p6z + 
	       p3x*p4y*p6z) + 
          p3z*(-(p4z*p5y*p6x) + p4z*p5x*p6y - p4y*p5x*p6z + 
	       p4x*p5y*p6z)));
  czte2=((e6*(-(p3z*p4y*p5x) + p3y*p4z*p5x + p3z*p4x*p5y - 
	      p3x*p4z*p5y - p3y*p4x*p5z + p3x*p4y*p5z) + 
          e4*(-(p3z*p5y*p6x) + p3y*p5z*p6x + p3z*p5x*p6y - 
	      p3x*p5z*p6y - p3y*p5x*p6z + p3x*p5y*p6z))/
	 (p5z*(p3y*p4z*p6x - p3x*p4z*p6y - p3y*p4x*p6z + 
	       p3x*p4y*p6z) + 
          p3z*(-(p4z*p5y*p6x) + p4z*p5x*p6y - p4y*p5x*p6z + 
	       p4x*p5y*p6z)));
   
  cx=((del3b*p4z*(p3z*p5y - p3y*p5z) - 
       del31*p4z*p5z*p6y - del2b*p3z*p5y*p6z - 
       del31*p3z*p5y*p6z + del2b*p3y*p5z*p6z + 
       del31*p3y*p5z*p6z + del31*p4y*p5z*p6z + 
       del531*(-(p3y*p4z*p5z) + 
	       p3z*(p4z*(p5y + p6y) - p4y*p6z)) - 
       2*p3z*p4z*p5y*p6x*pvisx + 2*p3y*p4z*p5z*p6x*pvisx + 
       2*p3z*p4x*p5y*p6z*pvisx - 2*p3y*p4x*p5z*p6z*pvisx - 
       2*p3z*p4z*p5y*p6y*pvisy + 2*p3y*p4z*p5z*p6y*pvisy + 
       2*p3z*p4y*p5y*p6z*pvisy - 2*p3y*p4y*p5z*p6z*pvisy)/
      (2*(p5z*(-(p3y*p4z*p6x) + p3x*p4z*p6y + p3y*p4x*p6z - 
               p3x*p4y*p6z) + 
	  p3z*(p4z*p5y*p6x - p4z*p5x*p6y + p4y*p5x*p6z - 
               p4x*p5y*p6z))));
   
  cy= ((del3b*p4z*(p3z*p5x - p3x*p5z) - 
	del31*p4z*p5z*p6x - del2b*p3z*p5x*p6z - 
	del31*p3z*p5x*p6z + del2b*p3x*p5z*p6z + 
	del31*p3x*p5z*p6z + del31*p4x*p5z*p6z + 
	del531*(-(p3x*p4z*p5z) + 
		p3z*(p4z*(p5x + p6x) - p4x*p6z)) - 
	2*p3z*p4z*p5x*p6x*pvisx + 2*p3x*p4z*p5z*p6x*pvisx + 
	2*p3z*p4x*p5x*p6z*pvisx - 2*p3x*p4x*p5z*p6z*pvisx - 
	2*p3z*p4z*p5x*p6y*pvisy + 2*p3x*p4z*p5z*p6y*pvisy + 
	2*p3z*p4y*p5x*p6z*pvisy - 2*p3x*p4y*p5z*p6z*pvisy)/
       (2.*(p5z*(p3y*p4z*p6x - p3x*p4z*p6y - p3y*p4x*p6z + 
		 p3x*p4y*p6z) + 
            p3z*(-(p4z*p5y*p6x) + p4z*p5x*p6y - p4y*p5x*p6z + 
		 p4x*p5y*p6z))));


  cz= ((del3b*p4z*(p3y*p5x - p3x*p5y) - 
	del31*p4z*p5y*p6x + del31*p4z*p5x*p6y - 
	del2b*p3y*p5x*p6z - del31*p3y*p5x*p6z - 
	del31*p4y*p5x*p6z + del2b*p3x*p5y*p6z + 
	del31*p3x*p5y*p6z + del31*p4x*p5y*p6z + 
	del531*(p3y*(p4z*(p5x + p6x) - p4x*p6z) - 
		p3x*(p4z*(p5y + p6y) - p4y*p6z)) - 
	2*p3y*p4z*p5x*p6x*pvisx + 2*p3x*p4z*p5y*p6x*pvisx + 
	2*p3y*p4x*p5x*p6z*pvisx - 2*p3x*p4x*p5y*p6z*pvisx - 
	2*p3y*p4z*p5x*p6y*pvisy + 2*p3x*p4z*p5y*p6y*pvisy + 
	2*p3y*p4y*p5x*p6z*pvisy - 2*p3x*p4y*p5y*p6z*pvisy)/
       (2.*(p5z*(-(p3y*p4z*p6x) + p3x*p4z*p6y + p3y*p4x*p6z - 
		 p3x*p4y*p6z) + 
            p3z*(p4z*p5y*p6x - p4z*p5x*p6y + p4y*p5x*p6z - 
		 p4x*p5y*p6z))));
  czt= (del3b*(p3z*p4y*p5x - p3y*p4z*p5x - p3z*p4x*p5y + 
	       p3x*p4z*p5y + p3y*p4x*p5z - p3x*p4y*p5z) + 
	del2b*p3z*p5y*p6x + del31*p3z*p5y*p6x + 
	del31*p4z*p5y*p6x - del2b*p3y*p5z*p6x - 
	del31*p3y*p5z*p6x - del31*p4y*p5z*p6x - 
	del2b*p3z*p5x*p6y - del31*p3z*p5x*p6y - 
	del31*p4z*p5x*p6y + del2b*p3x*p5z*p6y + 
	del31*p3x*p5z*p6y + del31*p4x*p5z*p6y + 
	del2b*p3y*p5x*p6z + del31*p3y*p5x*p6z + 
	del31*p4y*p5x*p6z - del2b*p3x*p5y*p6z - 
	del31*p3x*p5y*p6z - del31*p4x*p5y*p6z + 
	del531*(p3z*(p4y*(p5x + p6x) - p4x*(p5y + p6y)) + 
		p3y*(-(p4z*(p5x + p6x)) + p4x*(p5z + p6z)) + 
		p3x*(p4z*(p5y + p6y) - p4y*(p5z + p6z))) - 
	2*p3z*p4y*p5x*p6x*pvisx + 2*p3y*p4z*p5x*p6x*pvisx - 
	2*p3x*p4z*p5y*p6x*pvisx + 2*p3x*p4y*p5z*p6x*pvisx + 
	2*p3z*p4x*p5x*p6y*pvisx - 2*p3x*p4x*p5z*p6y*pvisx - 
	2*p3y*p4x*p5x*p6z*pvisx + 2*p3x*p4x*p5y*p6z*pvisx - 
	2*p3z*p4y*p5y*p6x*pvisy + 2*p3y*p4y*p5z*p6x*pvisy + 
	2*p3y*p4z*p5x*p6y*pvisy + 2*p3z*p4x*p5y*p6y*pvisy - 
	2*p3x*p4z*p5y*p6y*pvisy - 2*p3y*p4x*p5z*p6y*pvisy - 
	2*p3y*p4y*p5x*p6z*pvisy + 2*p3x*p4y*p5y*p6z*pvisy - 
	2*p3z*p4z*p5y*p6x*pvisz + 2*p3y*p4z*p5z*p6x*pvisz + 
	2*p3z*p4z*p5x*p6y*pvisz - 2*p3x*p4z*p5z*p6y*pvisz - 
	2*p3z*p4y*p5x*p6z*pvisz + 2*p3z*p4x*p5y*p6z*pvisz - 
	2*p3y*p4x*p5z*p6z*pvisz + 2*p3x*p4y*p5z*p6z*pvisz)/
    (2*(p5z*(p3y*p4z*p6x - p3x*p4z*p6y - p3y*p4x*p6z + 
	     p3x*p4y*p6z) + 
	p3z*(-(p4z*p5y*p6x) + p4z*p5x*p6y - p4y*p5x*p6z + 
	     p4x*p5y*p6z)));
  double a11,b11,a22,b22,a12,b12,af,bf,cf,df,ef,a1,b1,a2,b2,a,b;
  a11=-1 + cxe1*cxe1 + cye1*cye1 + cze1*cze1;
  b11=cxe1*cxe1 + cye1*cye1 + cze1*cze1 - 2*cze1*czte1 + czte1*czte1;
  a22=cxe2*cxe2 + cye2*cye2 + cze2*cze2;
  b22=-1 + cxe2*cxe2 + cye2*cye2 + cze2*cze2 - 
    2*cze2*czte2 + czte2*czte2;
  a12=2*cxe1*cxe2 + 2*cye1*cye2 + 2*cze1*cze2;
  b12=    2*cxe1*cxe2 + 2*cye1*cye2 + 2*cze1*cze2 - 
    2*cze2*czte1 - 2*cze1*czte2 + 2*czte1*czte2;
  a1=2*cx*cxe1 + 2*cy*cye1 + 2*cz*cze1;
  b1=   2*cx*cxe1 + 2*cy*cye1 + 2*cz*cze1 - 2*cze1*czt - 
    2*cz*czte1 + 2*czt*czte1 + 2*cxe1*pvisx + 
    2*cye1*pvisy + 2*cze1*pvisz - 2*czte1*pvisz;
  a2=2*cx*cxe2 + 2*cy*cye2 + 2*cz*cze2;
  b2=    2*cx*cxe2 + 2*cy*cye2 + 2*cz*cze2 - 2*cze2*czt - 
    2*cz*czte2 + 2*czt*czte2 + 2*cxe2*pvisx + 
    2*cye2*pvisy + 2*cze2*pvisz - 2*czte2*pvisz;
  a=cx*cx + cy*cy + cz*cz + m1sq;
  b=    cx*cx + cy*cy + cz*cz - 2*cz*czt + czt*czt + 
    m2sq + 2*cx*pvisx + pvisx*pvisx + 2*cy*pvisy + 
    pvisy*pvisy + 2*cz*pvisz - 2*czt*pvisz + pvisz*pvisz;
  af= a11*a22*a22*b11*b11 - a11*a12*a22*b11*b12 + 
    a11*a11*a22*b12*b12 + a11*a12*a12*b11*b22 - 
    2*a11*a11*a22*b11*b22 - a11*a11*a12*b12*b22 + 
    a11*a11*a11*b22*b22;

  bf=    -(a11*a12*a22*b1*b11) + 2*a11*a2*a22*b11*b11 + 
    2*a11*a11*a22*b1*b12 - a11*a12*a2*b11*b12 - 
    a1*a11*a22*b11*b12 + a11*a11*a2*b12*b12 + 
    a11*a12*a12*b11*b2 - 2*a11*a11*a22*b11*b2 - 
    a11*a11*a12*b12*b2 - a11*a11*a12*b1*b22 + 
    2*a1*a11*a12*b11*b22 - 2*a11*a11*a2*b11*b22 - 
    a1*a11*a11*b12*b22 + 2*a11*a11*a11*b2*b22;
    
  cf= a11*a11*a22*b1*b1 + a11*a12*a12*b*b11 - 
    2*a11*a11*a22*b*b11 - a11*a12*a2*b1*b11 - 
    a1*a11*a22*b1*b11 + a11*a2*a2*b11*b11 + 
    2*a*a11*a22*b11*b11 - a11*a11*a12*b*b12 + 
    2*a11*a11*a2*b1*b12 - a*a11*a12*b11*b12 - 
    a1*a11*a2*b11*b12 + a*a11*a11*b12*b12 - 
    a11*a11*a12*b1*b2 + 2*a1*a11*a12*b11*b2 - 
    2*a11*a11*a2*b11*b2 - a1*a11*a11*b12*b2 + 
    a11*a11*a11*b2*b2 + 2*a11*a11*a11*b*b22 - 
    a1*a11*a11*b1*b22 + a1*a1*a11*b11*b22 - 
    2*a*a11*a11*b11*b22;

  df=  -(a11*a11*a12*b*b1) + a11*a11*a2*b1*b1 + 
    2*a1*a11*a12*b*b11 - 2*a11*a11*a2*b*b11 - 
    a*a11*a12*b1*b11 - a1*a11*a2*b1*b11 + 
    2*a*a11*a2*b11*b11 - a1*a11*a11*b*b12 + 
    2*a*a11*a11*b1*b12 - a*a1*a11*b11*b12 + 
    2*a11*a11*a11*b*b2 - a1*a11*a11*b1*b2 + 
    a1*a1*a11*b11*b2 - 2*a*a11*a11*b11*b2;
    
  ef= a11*(a11*b - a*b11)*(a11*b - a*b11) + 
    a1*(a11*b - a*b11)*(-(a11*b1) + a1*b11) + 
    a*(-(a11*b1) + a1*b11)*(-(a11*b1) + a1*b11);
   
  dd[0]=ef;   
  dd[1]=df;
  dd[2]=cf;
  dd[3]=bf;
  dd[4]=af;
  
  int nsolreal,nerr;
  double e1,e2,detjacf;
  quartic(dd,sol,soli,nsolreal);
   
  if (nsolreal==0) 
  {
    nerr=1;
    return false;
  };
   
  for(i=0;i<nsolreal;i++)
  {
    e2=sol[i];
    // substitute e2 solution into general e1 solution in terms of e2
    e1=
      (a11*b - a*b11 - a2*b11*e2 + a11*b2*e2 - 
       a22*b11*e2*e2 + a11*b22*e2*e2)/
      (-(a11*b1) + a1*b11 + a12*b11*e2 - a11*b12*e2);

    // now evaluate jacobians and momenta
    detjacf= -(a2*b1) + a1*b2 - a12*b1*e1 - 2*a2*b11*e1 + a1*b12*e1 + 
      2*a11*b2*e1 - 2*a12*b11*e1*e1 + 2*a11*b12*e1*e1 - 
      2*a22*b1*e2 - a2*b12*e2 + a12*b2*e2 + 2*a1*b22*e2 - 
      4*a22*b11*e1*e2 + 4*a11*b22*e1*e2 - 2*a22*b12*e2*e2 + 
      2*a12*b22*e2*e2;
    double jac=fabs(detjacf*detval);
      
    if(e1<=0.||e2<=0.)
    {
      nerr=3;
      continue;
    };
    // if(e1>2000.||e2>2000.) continue; 
    
    double p1x,p1y,p1z,p2x,p2y,p2z,ptotz,etot;
    p1x = cxe1*e1 + cxe2*e2 + cx;
    p1y = cye1*e1 + cye2*e2 + cy;
    p1z = cze1*e1 + cze2*e2 + cz;
    ptotz = czte1*e1 + czte2*e2 + czt;
    etot= e1+e2+e3+e4+e5+e6;
    pg1[nsolutions][0]=(etot+ptotz)/2.;
    pg2[nsolutions][0]=(etot-ptotz)/2.;
    pg1[nsolutions][1]=0.;
    pg2[nsolutions][1]=0.;
    pg1[nsolutions][2]=0.;
    pg2[nsolutions][2]=0.;
    pg1[nsolutions][3]=pg1[nsolutions][0];
    pg2[nsolutions][3]=-pg2[nsolutions][0];

    p2x = -pvisx - p1x;
    p2y = -pvisy - p1y;
    p2z = ptotz - p1z - pvisz;

    p1[nsolutions][0]=e1;
    p1[nsolutions][1]=p1x;
    p1[nsolutions][2]=p1y;
    p1[nsolutions][3]=p1z;     

    p2[nsolutions][0]=e2;
    p2[nsolutions][1]=p2x;
    p2[nsolutions][2]=p2y;
    p2[nsolutions][3]=p2z;

    for(int n=0;n<=3;n++)
    {
      p531[nsolutions][n]=p1[nsolutions][n]+p3[n]+p5[n];
      p642[nsolutions][n]=p2[nsolutions][n]+p4[n]+p6[n];
      p31[nsolutions][n]=p1[nsolutions][n]+p3[n];
      p42[nsolutions][n]=p2[nsolutions][n]+p4[n];
    };
    totaljac[nsolutions]=jac;
    nsolutions++;  
    if (solveone) return true;
  };
  if (nsolutions>0) return true;
  else return false;
};

bool Solver::allcombi()
{
  order(3,4,5,6);
  if (solve()) return true;
  order(5,4,3,6);
  if (solve()) return true;
  order(3,6,5,4);
  if (solve()) return true;
  order(5,6,3,4);
  if (solve()) return true;
  order(3,4,6,5);
  if (solve()) return true;
  order(6,4,3,5);
  if (solve()) return true;
  order(3,5,6,4);
  if (solve()) return true;
  order(6,5,3,4);
  if (solve()) return true; 
  return false;  
};

bool Solver::halfcombi()
{
  order(3,4,5,6);
  if (solve()) return true;
  order(5,4,3,6);
  if (solve()) return true;
  order(3,6,5,4);
  if (solve()) return true;
  order(5,6,3,4);
  if (solve()) return true;
  return false;  
};

void Solver::setMomenta( double *p3a,
			 double *p4a,
			 double *p5a,
			 double *p6a,
			 double *pmissa )
{
  int i;
  for (i=0;i<4;i++)
  { 
      
    p3o[i]=p3a[i];
    p4o[i]=p4a[i];      
    p5o[i]=p5a[i];      
    p6o[i]=p6a[i];
    pmiss[i]=pmissa[i];
  };
  double m3temp,m4temp,m5temp,m6temp;
  m3temp=p3o[0]*p3o[0]-p3o[1]*p3o[1]-p3o[2]*p3o[2]-p3o[3]*p3o[3];
  m4temp=p4o[0]*p4o[0]-p4o[1]*p4o[1]-p4o[2]*p4o[2]-p4o[3]*p4o[3];
  m5temp=p5o[0]*p5o[0]-p5o[1]*p5o[1]-p5o[2]*p5o[2]-p5o[3]*p5o[3];
  m6temp=p6o[0]*p6o[0]-p6o[1]*p6o[1]-p6o[2]*p6o[2]-p6o[3]*p6o[3];
  if(m3temp<0.) m3o=0; else m3o=sqrt(m3temp);
  if(m4temp<0.) m4o=0; else m4o=sqrt(m4temp);
  if(m5temp<0.) m5o=0; else m5o=sqrt(m5temp);
  if(m6temp<0.) m6o=0; else m6o=sqrt(m6temp);
  order(3,4,5,6);
  // m7=sqrt(p7[0]*p7[0]-p7[1]*p7[1]-p7[2]*p7[2]-p7[3]*p7[3]);
   
};

void Solver::setMasses( double m531a,
			double m642a,
			double m31a,
			double m42a,
			double m1a,
			double m2a )

{
  m531=m531a;
  m642=m642a;
  m31=m31a;
  m42=m42a;
  m1=m1a;
  m2=m2a;
};

void Solver::swap()
{
  double temp;
  for(int i=0;i<4;i++)
  {
    double temp=p3[i];
    p3[i]=p4[i];
    p4[i]=temp;
  };
  temp=m3;
  m3=m4;
  m4=temp;
};  

void Solver::order(int ip3,int ip4, int ip5, int ip6)
{
  if (ip3==3) {p3=p3o;m3=m3o;} 
  else if(ip3==4) {p3=p4o;m3=m4o;}
  else if(ip3==5) {p3=p5o;m3=m5o;}
  else if(ip3==6) {p3=p6o;m3=m6o;};
  if (ip4==3) {p4=p3o;m4=m3o;} 
  else if(ip4==4) {p4=p4o;m4=m4o;}
  else if(ip4==5) {p4=p5o;m4=m5o;}
  else if(ip4==6) {p4=p6o;m4=m6o;};
  if (ip5==3) {p5=p3o;m5=m3o;} 
  else if(ip5==4) {p5=p4o;m5=m4o;}
  else if(ip5==5) {p5=p5o;m5=m5o;}
  else if(ip5==6) {p5=p6o;m5=m6o;};
  if (ip6==3) {p6=p3o;m6=m3o;} 
  else if(ip6==4) {p6=p4o;m6=m4o;}
  else if(ip6==5) {p6=p5o;m6=m5o;}
  else if(ip6==6) {p6=p6o;m6=m6o;};
}

void Solver::printl()
{
  cout<<"inputs:"<<endl;
  cout<<"p3,p4,p5,p6,pmiss"<<endl;
  cout<<p3[0]<<"  \t"<<p3[1]<<"  \t"<<p3[2]<<"  \t"<<p3[3]<<"  \t"<<m3<<endl;
  cout<<p4[0]<<"  \t"<<p4[1]<<"  \t"<<p4[2]<<"  \t"<<p4[3]<<"  \t"<<m4<<endl;
  cout<<p5[0]<<"  \t"<<p5[1]<<"  \t"<<p5[2]<<"  \t"<<p5[3]<<"  \t"<<m5<<endl;
  cout<<p6[0]<<"  \t"<<p6[1]<<"  \t"<<p6[2]<<"  \t"<<p6[3]<<"  \t"<<m6<<endl;
  cout<<pmiss[0]<<"  \t\t"<<pmiss[1]<<"  \t"<<pmiss[2]<<"  \t"<<pmiss[3]<<"  \t"<<endl;
  
  cout<<"solved="<<solved<<endl;
  cout<<"solutions: "<<nsolutions<<endl;
  if (solved) 
  {  
    for(int i=0;i<nsolutions;i++)
    {
      cout<<"solution "<<i<<endl;
      cout<<"p1,p2"<<endl;
      cout<<p1[i][0]<<"  \t"<<p1[i][1]<<"  \t"
	  <<p1[i][2]<<"  \t"<<p1[i][3]<<"  \t"<<m1<<endl;
      cout<<p2[i][0]<<"  \t"<<p2[i][1]<<"  \t"
	  <<p2[i][2]<<"  \t"<<p2[i][3]<<"  \t"<<m2<<endl;
    }
  }      
};

void Solver::set1d()
{ 
  if (solved)
  {
    int m;
    for(m=0; m<16; m++)
    {
      p1n[m]=0;
      p2n[m]=0;
    }
    i = 4;
    j = 0;
    for (col  = 0; col < nsolutions; col++) 
    {   
      for (row = 0; row < i; row++) 
      {
	p1n[j] =p1[col][row];
	p2n[j] =p2[col][row];
	j++;
      }
    }
  }
};
