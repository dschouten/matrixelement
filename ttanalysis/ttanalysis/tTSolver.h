#ifndef tTSolver_class_h
#define tTSolver_class_h

#include <vector>

using namespace std;



class tTSolver {
 public:
  
  tTSolver();
  ~tTSolver() { } // { cout << "tTSolver_class: destructor executing!" << endl; };

  void solve(double* ETmiss, double* b, double* bb, double* lp, double* lm, 
	     double mWp, double mWm, double mt, double mtb,
	     double mnu, double mnub,
	     vector<double> *pnux, vector<double> *pnuy, vector<double> *pnuz, 
	     vector<double> *pnubx, vector<double> *pnuby, vector<double> *pnubz,
	     vector<double> *cd_diff, int& cubic_single_root_cmplx);

  void quartic(vector<double> poly, vector<double> *pnuy, int& cubic_single_root_cmplx);
  void cubic(vector<double> poly, vector<double> *pnuy);
  void quadratic(vector<double> poly, vector<double> *pnuy);
  int algebraic_pz(double* b, double* lp, double mWp, double mt, double mb, double mlp, 
		   double mnu, double pnux, double pnuy, double* pnuz);
  double evalterma(vector<double> *a1, double pnux, double pnuy);
  double evaltermb(vector<double> *a2, double pnux, double pnuy);
  
 private: 
  
  static const double epsilon = 1.e-6; //numerical precision

}; //end class tTSolver  

#endif






#ifndef PI
#define PI fabs(acos(-1.))
#endif




inline
double sqr(double x) {
  return x*x;
}



inline
double sign(double a) {
  return (a < 0) ? -1 : (a > 0) ? 1 : 0;
}



inline
double quad(double x) {
  return x*x*x*x;
}

