
#
# this is a boosting code snippet that
# uses calorimeter MET only
#

CaloMET_x = itree.MET_x / GeV
CaloMET_y = itree.MET_y / GeV

vvZ = 125
vvM = 30

myMET = ROOT.TLorentzVector()

myMET.SetPxPyPzE( CaloMET_x, CaloMET_y, 0.,
                  math.sqrt( pow(CaloMET_x, 2) + pow(CaloMET_y, 2) + pow(vvZ, 2) + pow(vvM, 2) ) )

WWboost = -( measured_lp + measured_lm + myMET )

boost = WWboost.BoostVector()

if applyBoostCorrection:
    measured_lp.Boost( -boost.X(), -boost.Y(), 0 )
    measured_lm.Boost( -boost.X(), -boost.Y(), 0 )
    measured_jet.Boost( -boost.X(), -boost.Y(), 0 )
