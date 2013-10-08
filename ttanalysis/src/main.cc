#include <iostream>
#include <vector>


#include "ttdilepsolve.h"


int main() {

  ttdilepsolve fdilepsol;

  double lp[4], lm[4], b[4], bb[4];
  double ETmiss[2], nu[4], nub[4];

  //lepton l+
  lp[0]=27.8775;  //E 
  lp[1]=14.7347;  //px
  lp[2]=21.9598;  //py
  lp[3]=-8.8204; //pz

  //lepton l-
  lm[0]=168.71; //E
  lm[1]=-15.1858; //px
  lm[2]=-97.1134; //py
  lm[3]=137.118; //pz

  //b quark
  b[0]=105.321; //E
  b[1]=27.9652; //px
  b[2]=94.1283; //py
  b[3]=-37.7801; //pz

  //bbar quark
  bb[0]=76.3095; //E
  bb[1]=33.7778; //px
  bb[2]=35.977; //py
  bb[3]=58.007; //pz

  //neutrino
  nu[0]=64.8062; //E
  nu[1]=-28.2752; //px
  nu[2]=-54.4052; //py
  nu[3]=-20.9865; //pz
  
  //anti-neutrino
  nub[0]=88.8364; //E
  nub[1]=-33.7231; //px
  nub[2]=0.273651; //py
  nub[3]=82.1863; //pz
  
  //Ex miss
  ETmiss[0]=nu[1]+nub[1];
  
  //Ey miss
  ETmiss[1]=nu[2]+nub[2];
  
  
  vector<double> pnux, pnuy, pnuz, pnubx, pnuby, pnubz;
  vector<double> pnuychi2, pnunubzchi2, pnuyzchi2, cd_diff;
  int cubic_single_root_cmplx;

  double mt, mtb;
  mt = mtb = 175.;
  double mWp, mWm;
  mWp = mWm = 80.41;
  double mnu, mnub;
  mnu = mnub = 0.;

  cout << "generated (anti-)neutrino momenta:" << endl;
  cout << "pnux=" << nu[1] << " pnuy=" << nu[2] << " pnuz=" << nu[3]
       << " pnubx=" << nub[1] << " pnuby=" << nub[2] << " pnubz=" << nub[3] << endl;

  cout << "Exmiss=" << ETmiss[0] << " Eymiss=" << ETmiss[1] << endl;

  cout << "Invoking solve procedure now!" << endl;

  //cd_diff gives the value of this coefficient which can become singular (debugging purposes)
  //cubic_single_root_cmplx gives the number of single complex rootes of the cubic
  // due to numerical inaccuracies (for debugging purposes)
  fdilepsol.solve(ETmiss, b, bb, lp, lm, mWp, mWm, mt, mtb, mnu, mnub, 
		  &pnux, &pnuy, &pnuz, &pnubx, &pnuby, &pnubz, 
		  &cd_diff, cubic_single_root_cmplx);

  for (int i=0; i<pnux.size(); ++i) 
    cout << "solution[" << i << "]: pnux=" << pnux[i] << " pnuy=" << pnuy[i] << " pnuz=" << pnuz[i]
	 << " pnubx=" << pnubx[i] << " pnuby=" << pnuby[i] << " pnubz=" << pnubz[i] << endl;



  return 0;

}
