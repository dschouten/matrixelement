
execfile( 'configs/cfg_HWW125.py' )

myMECalculation.mehndl.setName( 'HWWBoosted' + str( int(mH) ) )

# Choose to apply boost correction using recoil estimate
applyBoostCorrection = True
boostCalculator = 'hwwBoostCalo3D.py'
