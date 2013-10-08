
#
# this is a boosting code snippet that
# uses calorimeter or track MET depending on average |MET|
#

TrackMET_x = -itree.MET_TrackHWW * math.cos(itree.MET_TrackHWW_phi) / GeV
TrackMET_y = -itree.MET_TrackHWW * math.sin(itree.MET_TrackHWW_phi) / GeV

CaloMET_x = itree.MET_x / GeV
CaloMET_y = itree.MET_y / GeV

avgMET_x = 0.5*(TrackMET_x+CaloMET_x)
avgMET_y = 0.5*(TrackMET_y+CaloMET_y)

preferredMET_x = 0
preferredMET_y = 0

if math.sqrt( math.pow(avgMET_x,2) + math.pow(avgMET_y,2)  ) > 35:
    preferredMET_x = CaloMET_x
    preferredMET_y = CaloMET_y
else:
    preferredMET_x = TrackMET_x
    preferredMET_y = TrackMET_y

vvZ = 125
vvM = 30

myMET = ROOT.TLorentzVector()

myMET.SetPxPyPzE( preferredMET_x, preferredMET_y, 0.,
                  math.sqrt( pow(preferredMET_x, 2) + pow(preferredMET_y, 2) + pow(vvZ, 2) + pow(vvM, 2) ) )

WWboost = -( measured_lp + measured_lm + myMET )

boost = WWboost.BoostVector()

if applyBoostCorrection:
    measured_lp.Boost( -boost.X(), -boost.Y(), 0 )
    measured_lm.Boost( -boost.X(), -boost.Y(), 0 )
    measured_jet.Boost( -boost.X(), -boost.Y(), 0 )
