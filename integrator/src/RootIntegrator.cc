//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "integrator/RootIntegrator.hh"

#include <iostream>
#include <stdexcept>
#include <vector>

#include <TMath.h>

using std::vector;

// ------------------------- ======= ------------------------- ======= -------------------------
RootIntegrator::RootIntegrator() :
  Integrator("ROOT integrator")
{

}

// ------------------------- ======= ------------------------- ======= -------------------------
void RootIntegrator::doIntegral(double returnVal[], double relerr[],
                                int* ifail, int* nfnevl, double prob[]) const
{

  cuba::integrand_t intg = getIntegrand();

  int n = getNDimensions();
  int comp = getNComp();

  double resultArr[1] = {0};

  int minpts = getMinEval();
  int maxpts = getMaxEval();
  double eps = getEpsilonRel();

  const double xl2 = 0.358568582800318073;
  const double xl4 = 0.948683298050513796;
  const double xl5 = 0.688247201611685289;
  const double w2  = 980./6561;
  const double w4  = 200./19683;
  const double wp2 = 245./486;
  const double wp4 = 25./729;

  const double wn1[14] = {-0.193872885230909911, -0.555606360818980835,
			  -0.876695625666819078, -1.15714067977442459,
			  -1.39694152314179743, -1.59609815576893754,
			  -1.75461057765584494, -1.87247878880251983,
			  -1.94970278920896201, -1.98628257887517146,
			  -1.98221815780114818, -1.93750952598689219,
			  -1.85215668343240347, -1.72615963013768225};

  const double wn3[14] = {0.0518213686937966768, 0.0314992633236803330,
			  0.0111771579535639891, -0.00914494741655235473,
			  -0.0294670527866686986, -0.0497891581567850424,
			  -0.0701112635269013768, -0.0904333688970177241,
			  -0.110755474267134071, -0.131077579637250419,
			  -0.151399685007366752, -0.171721790377483099,
			  -0.192043895747599447, -0.212366001117715794};

  const double wn5[14] = {0.871183254585174982e-01, 0.435591627292587508e-01,
			  0.217795813646293754e-01, 0.108897906823146873e-01,
			  0.544489534115734364e-02, 0.272244767057867193e-02,
			  0.136122383528933596e-02, 0.680611917644667955e-03,
			  0.340305958822333977e-03, 0.170152979411166995e-03,
			  0.850764897055834977e-04, 0.425382448527917472e-04,
			  0.212691224263958736e-04, 0.106345612131979372e-04};

  const double wpn1[14] = {-1.33196159122085045, -2.29218106995884763,
			   -3.11522633744855959, -3.80109739368998611,
			   -4.34979423868312742, -4.76131687242798352,
			   -5.03566529492455417, -5.17283950617283939,
			   -5.17283950617283939, -5.03566529492455417,
			   -4.76131687242798352, -4.34979423868312742,
			   -3.80109739368998611, -3.11522633744855959};

  const double wpn3[14] = {0.0445816186556927292, -0.0240054869684499309,
			   -0.0925925925925925875, -0.161179698216735251,
			   -0.229766803840877915, -0.298353909465020564,
			   -0.366941015089163228, -0.435528120713305891,
			   -0.504115226337448555, -0.572702331961591218,
			   -0.641289437585733882, -0.709876543209876532,
			   -0.778463648834019195, -0.847050754458161859};

  double result = 0;

  //   std::cerr << __LINE__ << " " << result << std::endl;

  double abserr = 0;
  *ifail = 3;
  *nfnevl = 0;
  relerr[0] = 0;
  if (n < 2 || n > 15)
    throw std::runtime_error("Improper dimension for RootIntegrator");

  int twondm = pow(2, n);
  int ifncls = 0;
  bool ldv = false;
  int irgnst = 2 * n + 3;
  int irlcls = twondm + 2 * n * (n + 1) + 1;
  int isbrgn = irgnst;
  int isbrgs = irgnst;

  if (minpts < 1) 
    minpts = irlcls;
  if (maxpts < minpts)
    maxpts = 10 * minpts;

  // The original agorithm expected a working space array WK of length IWK
  // with IWK Length ( >= (2N + 3) * (1 + MAXPTS/(2**N + 2N(N + 1) + 1))/2).
  // Here, this array is allocated dynamically

  int iwk = irgnst * (1 + maxpts / irlcls) / 2;
  double* wk = new double[iwk + 10];

  double* ctr = new double[n];
  double* wth = new double[n];
  double* z = new double[n];
  double* wthl = new double[n];

  for (int j = 0; j < n; ++j)
  {
    ctr[j] = 0.5;
    wth[j] = 0.5;
  }

  int idvax0 = 0;

  while (*ifail == 3)
  {
    double rgnvol = twondm;

    for (int j = 0; j < n; ++j)
    {
      rgnvol *= wth[j];
      z[j] = ctr[j];
    }

    //   std::cerr << __LINE__ << std::endl;
    intg(&n, z, &comp, resultArr, m_userdata);
    //   std::cerr << resultArr[0] << std::endl;
    double sum_a = resultArr[0]; //evaluate function

    double difmax = 0, sum_b = 0, sum_c = 0, f2, f3, dif;
    int idvaxn = 0;

    for (int j = 0; j < n; ++j)
    {
      z[j] = ctr[j] - xl2 * wth[j];
      //         std::cerr << "z[" << j << "]: " << z[j] << std::endl;
      intg(&n, z, &comp, resultArr, m_userdata);
      f2 = resultArr[0];

      z[j] = ctr[j] + xl2 * wth[j];
      //      std::cerr << __LINE__ << std::endl;
      intg(&n, z, &comp, resultArr, m_userdata);
      f2 += resultArr[0];

      wthl[j] = xl4 * wth[j];
      z[j] = ctr[j] - wthl[j];
      //      std::cerr << __LINE__ << std::endl;
      intg(&n, z, &comp, resultArr, m_userdata);
      f3 = resultArr[0];

      z[j] = ctr[j] + wthl[j];
      //      std::cerr << __LINE__ << std::endl;
      intg(&n, z, &comp, resultArr, m_userdata);
      f3 += resultArr[0];

      sum_b += f2;
      sum_c += f3;
      dif = TMath::Abs(7 * f2 - f3 - 12 * sum_a);
      if (dif >= difmax)
      {
	difmax = dif;
	idvaxn = j + 1;
      }
      z[j] = ctr[j];
    }

    double sum_d = 0;
    for (int j = 1; j < n; ++j)
    {
      int j1 = j-1;
      for (int k = j; k < n; ++k)
      {
	for (int l = 0; l < 2; ++l)
	{
	  wthl[j1] = -wthl[j1];
	  z[j1]    = ctr[j1] + wthl[j1];
	  for (int m = 0; m < 2; ++m)
	  {
	    wthl[k] = -wthl[k];
	    z[k]    = ctr[k] + wthl[k];
	    intg(&n, z, &comp, resultArr, m_userdata);
	    sum_d += resultArr[0];
	  }
	}
	z[k] = ctr[k];
      }
      z[j1] = ctr[j1];
    }

    double sum_e = 0;
    for (int j = 0; j < n; ++j)
    {
      wthl[j] = -xl5 * wth[j];
      z[j] = ctr[j] + wthl[j];
    }

    bool repeat = true;
    while (repeat)
    {
      intg(&n, z, &comp, resultArr, m_userdata);
      //   std::cerr << resultArr[0] << std::endl;
      sum_e += resultArr[0];
      
      repeat = false;
      for (int j = 0; j < n; ++j)
      {
	wthl[j] = -wthl[j];
	z[j] = ctr[j] + wthl[j];
	if (wthl[j] > 0)
	{
	  repeat = true;
	  break;
	}
      }
    }
      
    //   std::cerr <<wn1[n-2] << " " << sum_a << " " << w2 << " " << sum_b+wn3[n-2] << sum_c+w4 << sum_d+wn5[n-2] << sum_e << std::endl;
    //   std::cerr << __LINE__ << " " << result << std::endl;

    double rgncmp  = rgnvol * (wpn1[n - 2] * sum_a + wp2 * sum_b
			       + wpn3[n - 2] * sum_c + wp4 * sum_d);
    double rgnval  = wn1[n - 2] * sum_a + w2 * sum_b + wn3[n - 2] * sum_c
      + w4 * sum_d + wn5[n - 2] * sum_e;
    rgnval *= rgnvol;
    double rgnerr  = TMath::Abs(rgnval - rgncmp);
    result += rgnval;
    //   std::cerr << "rgnval: " << rgnval << " " << result <<std::endl;
    abserr += rgnerr;
    ifncls += irlcls;
    double aresult = TMath::Abs(result);

    // std::cerr << __LINE__ << " " << result << std::endl;
    // exit(1);

    // if (result > 0 && aresult< 1e-100) {
    //   delete [] wk;
    //   ifail = 0;  //function is probably symmetric ==> integral is null: not an error
    //   return result;
    // }

    int isbtmp = 0;

    bool doThisLoop = true;
    if (ldv)
    {
      //L110:
      while (1)
      {
	isbtmp = 2 * isbrgn;
	if (isbtmp > isbrgs)
	{
	  doThisLoop = false;
	  break;
	}
	if (isbtmp < isbrgs)
	{
	  int isbtpp = isbtmp + irgnst;
	  if (wk[isbtmp - 1] < wk[isbtpp - 1]) isbtmp = isbtpp;
	}
	if (rgnerr >= wk[isbtmp - 1])
	{
	  doThisLoop = false;
	  break;
	}

	for (int k = 0; k < irgnst; ++k)
	{
	  wk[isbrgn - k - 1] = wk[isbtmp - k - 1];
	}
	isbrgn = isbtmp;
      }
    }

    if (doThisLoop)
    {
      isbtmp = (isbrgn / (2 * irgnst)) * irgnst;
      while (isbtmp >= irgnst && rgnerr > wk[isbtmp - 1])
      {
	for (int k = 0; k < irgnst; ++k)
	{
	  wk[isbrgn - k - 1] = wk[isbtmp - k - 1];
	}
	isbrgn = isbtmp;
	isbtmp = (isbrgn / (2 * irgnst)) * irgnst;
      }
    }

    wk[isbrgn - 1] = rgnerr;
    wk[isbrgn - 2] = rgnval;
    wk[isbrgn - 3] = static_cast<double>(idvaxn);
    for (int j = 0; j < n; ++j)
    {
      isbtmp = isbrgn - 2 * j - 4;
      wk[isbtmp] = ctr[j];
      wk[isbtmp - 1] = wth[j];
    }
    if (ldv)
    {
      ldv = false;
      ctr[idvax0] += 2 * wth[idvax0];
      isbrgs += irgnst;
      isbrgn = isbrgs;
      continue;
    }
    relerr[0] = abserr/aresult;
    //      if (relerr[0] < 1e-1 && aresult < 1e-20)
    //         *ifail = 0;
    //      if (relerr[0] < 1e-3 && aresult < 1e-10)
    //         *ifail = 0;
    //      if (relerr[0] < 1e-5 && aresult < 1e-5)
    //         *ifail = 0;
    if (isbrgs + irgnst > iwk)
      *ifail = 2;
    if (ifncls + 2 * irlcls > maxpts)
    {
      if (sum_a == 0 && sum_b == 0 && sum_c == 0 && sum_d == 0 && sum_e == 0)
      {
	*ifail = 0;
	result = 0;
      }
      else
      {
	*ifail = 1;
      }
    }

    //      std::cerr << "Rel error: " << relerr[0] << std::endl;
    if (relerr[0] < eps && ifncls >= minpts)
      *ifail = 0;
    if (*ifail == 3)
    {
      ldv = true;
      isbrgn  = irgnst;
      abserr -= wk[isbrgn - 1];
      result -= wk[isbrgn - 2];
      idvax0  = Int_t(wk[isbrgn - 3] - 1);
      for (int j = 0; j < n; ++j) 
      {
	isbtmp = isbrgn - 2 * j - 4;
	ctr[j] = wk[isbtmp];
	wth[j] = wk[isbtmp - 1];
      }
      wth[idvax0]  = 0.5 * wth[idvax0];
      ctr[idvax0] -= wth[idvax0];
    }
  }

  delete [] wk;
  delete [] ctr;
  delete [] wth;
  delete [] z;
  delete [] wthl;

  *nfnevl = ifncls;       //number of function evaluations performed.
  returnVal[0] = result;         //an approximate value of the integral
  relerr[0] *= result;
}

