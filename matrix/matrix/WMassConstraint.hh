//
// authors: Doug Schouten (doug dot schouten at triumf dot ca)
//         
//

#ifndef WMASSCONSTRAINT_HH
#define WMASSCONSTRAINT_HH

template <typename T>
class WMassConstraint
{
public:

  static const T* _partons;

  static TLorentzVector _nu;
  static TLorentzVector _nub;

  static double _qwp_sqr;
  static double _qwm_sqr;

  static void function(Vec_I_DP &v, Vec_O_DP &f)
  {
    _nu.SetPxPyPzE( v[0], v[1], _partons->nul.Pz(), hepstd::norm( v[0], v[1], _partons->nul.Pz() ) );
    _nub.SetPxPyPzE( v[2], v[3], _partons->nur.Pz(), hepstd::norm( v[2], v[3], _partons->nur.Pz() ) );
      
    f[0] = _qwp_sqr - _partons->lp.M()*_partons->lp.M() - 2*_partons->lp.E()*hepstd::norm( v[0], v[1], _partons->nul.Pz() ) + 
      2*(_partons->lp.Px()*v[0] + _partons->lp.Py()*v[1] + _partons->lp.Pz()*_partons->nul.Pz());
    f[1] = _qwm_sqr - _partons->lm.M()*_partons->lm.M() - 2*_partons->lm.E()*hepstd::norm( v[2], v[3], _partons->nur.Pz() ) + 
      2*(_partons->lm.Px()*v[2] + _partons->lm.Py()*v[3] + _partons->lm.Pz()*_partons->nur.Pz());
    f[2] = _partons->met.Px() - v[0] - v[2];
    f[3] = _partons->met.Py() - v[1] - v[3];
  }

  static void jacobian(Vec_IO_DP &x, Vec_I_DP &fvec, Mat_O_DP &df, void vecfunc(Vec_I_DP &, Vec_O_DP &))
  {
    double nue  = hepstd::norm( x[0], x[1], _partons->nul.Pz() );
    double nube = hepstd::norm( x[2], x[3], _partons->nur.Pz() );
      
    df[0][0] = -2*_partons->lp.E()*x[0] / nue + 2*_partons->lp.Px();
    df[0][1] = -2*_partons->lp.E()*x[1] / nue + 2*_partons->lp.Py();
    df[0][2] = 0;
    df[0][3] = 0;
      
    df[1][0] = 0;
    df[1][1] = 0;
    df[1][2] = -2*_partons->lm.E()*x[2] / nube + 2*_partons->lm.Px();
    df[1][3] = -2*_partons->lm.E()*x[3] / nube + 2*_partons->lm.Py();
      
    df[2][0] = -1;
    df[2][1] = 0;
    df[2][2] = -1;
    df[2][3] = 0;
      
    df[3][0] = 0;
    df[3][1] = -1;
    df[3][2] = 0;
    df[3][3] = -1;
  }
};

template <typename T>
TLorentzVector WMassConstraint<T>::_nu;
template <typename T>
TLorentzVector WMassConstraint<T>::_nub;

template <typename T>
double WMassConstraint<T>::_qwp_sqr;
template <typename T>
double WMassConstraint<T>::_qwm_sqr;

template <typename T>
const T* WMassConstraint<T>::_partons;

#endif 
