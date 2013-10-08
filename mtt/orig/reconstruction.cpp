/* A translation for J. Gunion's reconstruction code from Fortran to C++ */
/* Changed to complex solver to keep complex solutions */
/* Author: Yang Bai and Zhenyu Han; Please refer the paper: arXiv:0809.4487 */

#include "reconstruction.h"
#include <iostream>
#include <math.h>
#include <complex>
using namespace std;

 double signnorm(dcmplx d)
{ 
  double sign;
  if(real(d)>0) sign=1;else sign=-1;
  return sign*sqrt(norm(d));
}
void quartic(double *dd,dcmplx x[4])
{
   double aa[4],z[3];
   double a,b,c,d,e;
   a = dd[4];
   b = dd[3];
   c = dd[2];
   d = dd[1];
   e = dd[0];
   if (fabs(dd[4])<10e-30)
   { 
     //cout<<"ERROR: NOT A QUARTIC EQUATION"<<endl;
     return;
   };   
   double p,q,r;
   p = (-3.*b*b + 8.*a*c)/(8.*a*a);
   q = (b*b*b - 4.*a*b*c + 8.*d*a*a)/(8.*a*a*a);
   r = (-3.*b*b*b*b + 16.*a*b*b*c - 64.*a*a*b*d + 
           256.*a*a*a*e)/(256.*a*a*a*a);

   dcmplx P,Q,R,U,y,W,pp,qq,rr;
   P=dcmplx(-p*p/12.-r,0);
   Q=dcmplx(-p*p*p/108.+p*r/3.-q*q/8.);
   pp=dcmplx(p,0);
   qq=dcmplx(q,0);
   rr=dcmplx(r,0);
   if(fabs(q)<1e-30) 
   {
      x[0]=-b/(4.*a)+sqrt(-pp+sqrt(pp*pp-rr*4.)/2.);
      x[1]=-b/(4.*a)+sqrt(-pp-sqrt(pp*pp-rr*4.)/2.);
      x[2]=-b/(4.*a)-sqrt(-pp+sqrt(pp*pp-rr*4.)/2.);
      x[3]=-b/(4.*a)-sqrt(-pp-sqrt(pp*pp-rr*4.)/2.);
   }
  else
   {
      R=Q/2.+sqrt(Q*Q/4.+P*P*P/27.);
      U=pow(R,1./3.);
      //cout<<"R="<<R<<endl;
      //cout<<"U="<<U<<endl;
      if(norm(U)<10.e-30) y=-5./6.*pp-U;
      else y=-5./6.*pp-U+P/(3.*U);
      W=sqrt(pp+2.*y);
      //cout<<"y="<<y<<endl;
      //cout<<"W="<<W<<endl;
      x[0]=-b/(4.*a)+(W+sqrt(-(3.*pp+2.*y+2.*qq/W)))/2.;
      x[1]=-b/(4.*a)+(W-sqrt(-(3.*pp+2.*y+2.*qq/W)))/2.;
      x[2]=-b/(4.*a)+(-W+sqrt(-(3.*pp+2.*y-2.*qq/W)))/2.;
      x[3]=-b/(4.*a)+(-W-sqrt(-(3.*pp+2.*y-2.*qq/W)))/2.;
    }
   
}
event::event()
{
   solved=false;
   solveone=true;
   nsolutions=0;
   debug0=0;
};
bool event::solve()
{
  solvecpl();

  nsolutions=0;
  realcut();  
  
  //count the number of real solutions
  solved=true;
  if (nsolutions>0) return true;
  else return false;
}
void event::solvecpl()   
{
   double dd[5],soli[4];
   dcmplx sol[4];
   double m3sq,m4sq,m5sq,m6sq,m531sq,m642sq,m31sq,m42sq,m1sq,m2sq;
   int i;
  
   nsolutions=0;
   m3sq=dot(p3,p3);
   m4sq=dot(p4,p4);
   m5sq=dot(p5,p5);
   m6sq=dot(p6,p6);
   m531sq=m531*m531;
   m642sq=m642*m642;
   m31sq=m31*m31;
   m42sq=m42*m42;
   m1sq=m1*m1;
   m2sq=m2*m2;
   
   double pvisx,pvisy,pvisz;
   pvisx=-pmiss[1];
   pvisy=-pmiss[2];
   pvisz=p3[3]+p4[3]+p5[3]+p6[3];
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
   
   double e3,p3x,p3y,p3z,e4,p4x,p4y,p4z,e5,p5x,p5y,p5z,e6,p6x,p6y,p6z;
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

   double detval,cxe1,cxe2,cye1,cye2,cze1,cze2,czte1,czte2,cx,cy,cz,czt;
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

//now e1 and e2 solutions
  
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
   //double e1,e2,detjacf;
   dcmplx e1,e2,detjacf;
   dcmplx x[4];
   for(i=0;i<4;i++) x[i]=dcmplx(0,-1000);
   quartic(dd,x);
   for(i=0;i<4;i++)
   {
      e2=x[i];
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
      
      dcmplx p1x,p1y,p1z,p2x,p2y,p2z,ptotz,etot;
      p1x = cxe1*e1 + cxe2*e2 + cx;
      p1y = cye1*e1 + cye2*e2 + cy;
      p1z = cze1*e1 + cze2*e2 + cz;
      ptotz = czte1*e1 + czte2*e2 + czt;
      etot= e1+e2+e3+e4+e5+e6;
      
      pg1cpl[nsolutions][0]=(etot+ptotz)/2.;
      pg2cpl[nsolutions][0]=(etot-ptotz)/2.;
      pg1cpl[nsolutions][1]=0.;
      pg2cpl[nsolutions][1]=0.;
      pg1cpl[nsolutions][2]=0.;
      pg2cpl[nsolutions][2]=0.;
      pg1cpl[nsolutions][3]=pg1[nsolutions][0];
      pg2cpl[nsolutions][3]=-pg2[nsolutions][0];

      p2x = -pvisx - p1x;
      p2y = -pvisy - p1y;
      p2z = ptotz - p1z - pvisz;
      
      p1cpl[nsolutions][0]=e1;
      p1cpl[nsolutions][1]=p1x;
      p1cpl[nsolutions][2]=p1y;
      p1cpl[nsolutions][3]=p1z;
     
      
      p2cpl[nsolutions][0]=e2;
      p2cpl[nsolutions][1]=p2x;
      p2cpl[nsolutions][2]=p2y;
      p2cpl[nsolutions][3]=p2z;

      for(int n=0;n<=3;n++)
      {
         p531cpl[nsolutions][n]=p1cpl[nsolutions][n]+p3[n]+p5[n];
         p642cpl[nsolutions][n]=p2cpl[nsolutions][n]+p4[n]+p6[n];
         p31cpl[nsolutions][n]=p1cpl[nsolutions][n]+p3[n];
         p42cpl[nsolutions][n]=p2cpl[nsolutions][n]+p4[n];
      };

      nsolutions++;  
      if (solveone) return;
   }

};
void event::realcut()
{
  double cmplxcut=0.4;
  nsolutions=0;
   for (int i=0;i<4;i++)
   {  
     
     //if(fabs(imag(p1cpl[i][0]))<cmplxcut*fabs(real(p1cpl[i][0]))&&
     // fabs(imag(p2cpl[i][0]))<cmplxcut*fabs(real(p2cpl[i][0])))
     { 
        //cout<<x[i]<<endl;
       dcmplx ptt[4];
       double pttreal[4];
       //don't count complex conjugates as two solutions
       int dupe=0;
       for(int j=0;j<i;j++)
        if(fabs(sqrt(norm(p1cpl[i][0]))-sqrt(norm(p1cpl[j][0])))<0.0001)
          dupe=1;
       if (dupe==1) continue;
       for(int j=0;j<4;j++)
       {
         
          p1[nsolutions][j]=signnorm(p1cpl[i][j]);
          p2[nsolutions][j]=signnorm(p2cpl[i][j]);
          p31[nsolutions][j]=signnorm(p31cpl[i][j]);
          p42[nsolutions][j]=signnorm(p42cpl[i][j]);
	  p531[nsolutions][j]=signnorm(p531cpl[i][j]);
          p642[nsolutions][j]=signnorm(p642cpl[i][j]);
          ptt[j]=p531cpl[i][j]+p642cpl[i][j]; 
         
       }
       //if(p531[nsolutions][0]<0||p642[nsolutions][0]<0) continue;
       p531[nsolutions][0]=
	 sqrt(pow(p531[nsolutions][1],2)+pow(p531[nsolutions][2],2)+pow(p531[nsolutions][3],2)+m531*m531);
       p642[nsolutions][0]=
	 sqrt(pow(p642[nsolutions][1],2)+pow(p642[nsolutions][2],2)+pow(p642[nsolutions][3],2)+m642*m642);
       for(int j=0;j<4;j++)
        pttreal[j]=p531[nsolutions][j]+p642[nsolutions][j];
       m2res[nsolutions]=signnorm(dot(ptt,ptt));
       //cout<<"dot(ptt,ptt) = "<<dot(ptt,ptt) <<endl;
       //cout<<"dot(pttreal,pttreal) = "<<dot(pttreal,pttreal) <<endl;
       if (fabs(imag(p531cpl[i][0]))>cmplxcut*fabs(real(p531cpl[i][0]))||
	   fabs(imag(p642cpl[i][0]))>cmplxcut*fabs(real(p642cpl[i][0]))) continue;
       if(real(p531cpl[i][0])<0||real(p642cpl[i][0])<0) continue;
       //if(dot(p1[nsolutions],p1[nsolutions])<-10000.
       // ||dot(p2[nsolutions],p2[nsolutions])<-10000.) continue;
       if( p1[nsolutions][0]<0||p2[nsolutions][0]<0) continue;
       nsolutions++;
     }
   }         
}

/*void event::realcut()
{
  double cmplxcut=1;
  nsolutions=0;
   for (int i=0;i<4;i++)
   {  
     
     if(fabs(imag(p1cpl[i][0]))<cmplxcut*fabs(real(p1cpl[i][0]))&&
        fabs(imag(p2cpl[i][0]))<cmplxcut*fabs(real(p2cpl[i][0])))
     { 
        //cout<<x[i]<<endl;
      
       for(int j=0;j<4;j++)
       {
          p1[nsolutions][j]=signnorm(p1cpl[i][j]);
          p2[nsolutions][j]=signnorm(p2cpl[i][j]);
          p31[nsolutions][j]=p1[nsolutions][j]+p3[j];
          p42[nsolutions][j]=p2[nsolutions][j]+p4[j];
          p531[nsolutions][j]=p1[nsolutions][j]+p3[j]+p5[j];
          p642[nsolutions][j]=p2[nsolutions][j]+p4[j]+p6[j];
       }
       if(dot(p1[nsolutions],p1[nsolutions])<-10000.
          ||dot(p2[nsolutions],p2[nsolutions])<-10000.) continue;
       nsolutions++;
     }
   }         
   }*/
bool event::solve_allcombi()
{
   double * ptemp[7]; int kftemp[7];
   ptemp[3]=p3o;ptemp[4]=p4o;ptemp[5]=p5o;ptemp[6]=p6o;
   kftemp[3]=kf3; kftemp[4]=kf4; kftemp[5]=kf5; kftemp[6]=kf6;
   int i,j,k,l;
   for (i=3;i<=6;i++)
   for (j=3;j<=6;j++)
   for (k=3;k<=6;k++)
   for (l=3;l<=6;l++)
     if( i!=j&&i!=l&&j!=k&&k!=l&&kftemp[i]+kftemp[k]==0&&kftemp[j]+kftemp[l]==0&& kftemp[j]!=3 && kftemp[l]!=3)
     {
        order(i,j,k,l);
        if (solve()) return true;
     }
   return false;
};
void event::setmomenta(double *p3a,double *p4a,double *p5a,double *p6a,double *pmissa)

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
   order(3,4,5,6);
};
void event::setkf(int kf3a,int kf4a,int kf5a,int kf6a)
{
    kf3=kf3a;kf4=kf4a;kf5=kf5a;kf6=kf6a;
};


void event::setmasses(double m531a,double m642a,double m31a,double m42a,double m1a,double m2a)

{
   m531=m531a;
   m642=m642a;
   m31=m31a;
   m42=m42a;
   m1=m1a;
   m2=m2a;
};

void event::order(int ip3,int ip4, int ip5, int ip6)
{
   if (ip3==3) p3=p3o; 
   else if(ip3==4) p3=p4o;
   else if(ip3==5) p3=p5o;
   else if(ip3==6) p3=p6o;
   if (ip4==3) p4=p3o;
   else if(ip4==4) p4=p4o;
   else if(ip4==5) p4=p5o;
   else if(ip4==6) p4=p6o;
   if (ip5==3) p5=p3o;
   else if(ip5==4) p5=p4o;
   else if(ip5==5) p5=p5o;
   else if(ip5==6) p5=p6o;
   if (ip6==3) p6=p3o;
   else if(ip6==4) p6=p4o;
   else if(ip6==5) p6=p5o;
   else if(ip6==6) p6=p6o;
     
}

void event::print()
{
   cout<<"inputs:"<<endl;
   cout<<"p3,p4,p5,p6,pmiss"<<endl;
   /*cout<<"double p3[4]={"<<p3[0]<<",\t"<<p3[1]<<",\t"<<p3[2]<<",\t"
       <<p3[3]<<"},\t"<<sqrt(dot(p3,p3))<<endl;
   cout<<"double p4[4]={"<<p4[0]<<",\t"<<p4[1]<<",\t"<<p4[2]<<",\t"
       <<p4[3]<<"},\t"<<sqrt(dot(p4,p4))<<endl;
   cout<<"double p5[4]={"<<p5[0]<<",\t"<<p5[1]<<",\t"<<p5[2]<<",\t"
       <<p5[3]<<"},\t"<<sqrt(dot(p5,p5))<<endl;
   cout<<"double p6[4]={"<<p6[0]<<",\t"<<p6[1]<<",\t"<<p6[2]<<",\t"
   <<p6[3]<<"},\t"<<sqrt(dot(p6,p6))<<endl; */
   cout<<"double p3[4]="<<p3[0]<<",\t"<<p3[1]<<",\t"<<p3[2]<<",\t"
       <<p3[3]<<",\t"<<sqrt(dot(p3,p3))<<endl;
   cout<<"double p4[4]="<<p4[0]<<",\t"<<p4[1]<<",\t"<<p4[2]<<",\t"
       <<p4[3]<<",\t"<<sqrt(dot(p4,p4))<<endl;
   cout<<"double p5[4]="<<p5[0]<<",\t"<<p5[1]<<",\t"<<p5[2]<<",\t"
       <<p5[3]<<",\t"<<sqrt(dot(p5,p5))<<endl;
   cout<<"double p6[4]="<<p6[0]<<",\t"<<p6[1]<<",\t"<<p6[2]<<",\t"
       <<p6[3]<<",\t"<<sqrt(dot(p6,p6))<<endl; 
   cout<<pmiss[0]<<"  \t\t"<<pmiss[1]<<"  \t"
       <<pmiss[2]<<"  \t"<<pmiss[3]<<"  \t"<<endl;
  
   cout<<"solved = "<<solved<<endl;  
   if (solved) 
   {  
      cout<<"Number of real solutions: "<<nsolutions<<endl;
      for(int i=0;i<nsolutions;i++)
      {
        cout<<"solution "<<i<<endl;
        cout<<"p1,p2"<<endl;
        cout<<p1[i][0]<<"  \t"<<p1[i][1]<<"  \t"
            <<p1[i][2]<<"  \t"<<p1[i][3]<<"  \t"<<m1<<endl;
        cout<<p2[i][0]<<"  \t"<<p2[i][1]<<"  \t"
            <<p2[i][2]<<"  \t"<<p2[i][3]<<"  \t"<<m2<<endl;
        double ptt[4];
        for(int j=0;j<4;j++) ptt[j]=p531[i][j]+p642[i][j];
        cout<<"mtt="<<sqrt(dot(ptt,ptt))<<endl;
      }
   }      
   
};
void event::printcpl()
{
  cout<<"Complex solutions: "<<endl;
  for(int i=0;i<4;i++)
  {
    cout<<"solution "<<i<<endl;
    cout<<"p1 = "<<p1cpl[i][0]<<"  \t"<<p1cpl[i][1]<<"  \t"
        <<p1cpl[i][2]<<"  \t"<<p1cpl[i][3]<<endl;
    cout<<"p2 = "<<p2cpl[i][0]<<"  \t"<<p2cpl[i][1]<<"  \t"
        <<p2cpl[i][2]<<"  \t"<<p2cpl[i][3]<<endl;
  }
       
};
void event::printfile(ofstream & outputfile,int i)
{
   
   outputfile<<p1[i][0]<<"  \t"<<p1[i][1]<<"  \t"
          <<p1[i][2]<<"  \t"<<p1[i][3]<<endl;
   outputfile<<p2[i][0]<<"  \t"<<p2[i][1]<<"  \t"
          <<p2[i][2]<<"  \t"<<p2[i][3]<<"  \t"<<endl;
   outputfile<<p3[0]<<"  \t"<<p3[1]<<"  \t"<<p3[2]<<"  \t"<<p3[3]<<"  \t"<<endl;
   outputfile<<p4[0]<<"  \t"<<p4[1]<<"  \t"<<p4[2]<<"  \t"<<p4[3]<<"  \t"<<endl;
   outputfile<<p5[0]<<"  \t"<<p5[1]<<"  \t"<<p5[2]<<"  \t"<<p5[3]<<"  \t"<<endl;
   outputfile<<p6[0]<<"  \t"<<p6[1]<<"  \t"<<p6[2]<<"  \t"<<p6[3]<<"  \t"<<endl;
   
};
