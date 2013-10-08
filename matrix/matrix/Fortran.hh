#ifndef FORTRAN_HH
#define FORTRAN_HH

#ifndef __CINT__

#include "matrix/Complex.hh"
#include "matrix/Constants.hh"
#include "matrix/Utility.hh"

/* CONS(a,b) should return ab, the concatenation
   of its arguments */

#if  __STDC__ || __APOGEE__
#define CONS(a,b) a##b
#else
#define CONS(a,b) a/**/b
#endif

#ifndef FORTRAN
#define FORTRAN(lcname,ucname)  CONS(lcname,_)
#endif

#define C_masses FORTRAN(masses,MASSES)
#define C_widths FORTRAN(widths,WIDTHS)
#define C_to_mirror FORTRAN(to_mirror,TO_MIRROR)
#define C_polarization FORTRAN(to_polarization,TO_POLARIZATION)
#define C_to_matrix FORTRAN(to_matrix,TO_MATRIX)

///////////////////////////////////////////////
//
// WARNING order matters for the common blocks
//
///////////////////////////////////////////////

extern struct _MASSES { // COMMON/MASSES/ MB,MH,MT,MW,MTA,MZ
  double MB; 
  double MH;  
  double MT;
  double MW;
  double MZ;
} C_masses;

extern struct _WIDTHS { // COMMON/WIDTHS/ WW,WT,WZ,WH
  double WW; 
  double WT; 
  double WZ; 
  double WH;
} C_widths;

extern struct _TO_MIRROR { // COMMON/TO_MIRROR/ IMIRROR
  int imirror;       //< should be 1
} C_to_mirror;

extern struct _POLARIZATION { // COMMON/TO_POLARIZATION
  double pol[2];     //< should be {1,1}
} C_polarization;

extern struct _TO_MATRIX { // COMMON/TO_MATRIX/ ISUM_HEL, MULTI_CHANNEL
  int isum_hel;      //< should be 0
  int multi_channel; //< should be 0 (false)
} C_to_matrix;

#define C_smcouplings FORTRAN(smcouplings,SMCOUPLINGS) // 
extern struct _SMCOUPLINGS { // COMMON/SMCOUPLINGS/
  dbl_cpx GC_1;
  dbl_cpx GC_16;
  dbl_cpx GC_2;
  dbl_cpx GC_21;
  dbl_cpx GC_22;
  dbl_cpx GC_23;
  dbl_cpx GC_24;
  dbl_cpx GC_25;
  dbl_cpx GC_3;
  dbl_cpx GC_31;
  dbl_cpx GC_33;
  dbl_cpx GC_4;
  dbl_cpx GC_5;
  dbl_cpx GC_7;
}  C_smcouplings; // ## END ## 

#define C_heftcouplings FORTRAN(heftcouplings,HEFTCOUPLINGS) // 
extern struct _HEFTCOUPLINGS { // COMMON/HEFTCOUPLINGS/
  dbl_cpx GC_17;
  dbl_cpx GC_37;
  dbl_cpx GC_5;
  dbl_cpx GC_6;
  dbl_cpx GC_8;
  dbl_cpx GC_9;
}  C_heftcouplings; // ## END ## 

#define C_rscouplings FORTRAN(rscouplings,RSCOUPLINGS) // 
extern struct _RSCOUPLINGS { // COMMON/RSCOUPLINGS/
  dbl_cpx GC_11;
  dbl_cpx GC_26;
  dbl_cpx GC_68;
}  C_rscouplings; // ## END ##

#define C_tomecomb FORTRAN(tomecomb,TOMECOMB) //
extern struct _TOMECOMB { // COMMON/TOMECOMB/
  double mecomb[NCOMBSTORED];
  int iselectedhel;
}  C_tomecomb; // ## END ##

#endif

namespace hepstd {
  void prepareCommonBlocks( double MH = 125.0, double MT=hepstd::tMass, double WH=0.001 ); 
}

#endif
