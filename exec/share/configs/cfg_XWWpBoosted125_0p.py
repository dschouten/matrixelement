import ROOT
import sys
import math

from matrixelement import toolFactory
from matrixelement import getDefaultIntegrator
from matrixelement import getAlternateIntegrator
from matrixelement import getVegasIntegrator
from matrixelement import getMiserIntegrator
from matrixelement import getHiggsWidth

# Set the Higgs mass & width for H->WW matrix elements
mH = 125
wH = getHiggsWidth( mH )

allowNull            = True  # allow PH == 0, otherwise assume non-convergence and try with finer grid
useNarrowWidthApprox = True  # remove integration over mH
useBoxPhaseSpace     = True  # if true, use q -> atan(t) transformation to sample BW
useSMKludge          = False # flag to use the SM kludge 
useRSModel           = False # flag to use RS graviton model (spin 2)
useJHUModel          = True  # flag to use JHU s=0/2 ME 
selectHelicity       = False # flag to select particular helicity combination (has no effect if not implemented in ME)

iHel = 0  # integrate only using this helicity combination (index starts @ 1) 
spin = 0  # specify spin of resonance in pp -> X -> WW(lv,lv) !! JHU only !!

cuts = None

nGam = 10 # integration range around BW peaks in units of natural width

# SystemBoostTF/recoil.tf, SystemBoostAverageTF/recoilavg_h125.tf, SystemBoostHybrid/recoilhyb_h125.tf
recoilTF = getattr( ROOT, 'SystemBoostAverageTF' )( 'recoilavg_h125.tf' ) 

# Configure default/backup/failsafe integrators
myMEIntegrators = { 'default'  : getVegasIntegrator( 0.015, 7.5e4, 1.5e5 ),
                    'backup'   : getVegasIntegrator( 0.025, 1.5e5, 3.0e5 ),
                    'failsafe' : getVegasIntegrator( 0.050, 3.0e5, 6.0e5 ) }

# Configure the ME calculators
myMECalculation = toolFactory( 'HWW' )( integrator = myMEIntegrators['default'],
                                        args = ( mH, ROOT.HWW.kWMASS_6D, recoilTF, useNarrowWidthApprox ) )

# Set integration limits (use TF to set dynamic limits for jets)
if not useNarrowWidthApprox:
    myMECalculation.mehndl.setWidth( wH )

mehndl = myMECalculation.mehndl

mehndl.setUseBoxPS( useBoxPhaseSpace )
mehndl.setUseSM( useSMKludge )
mehndl.setUseRS( useRSModel )

# Gets integration limits for transformed Q^{2} variable
def solveq( m, w, qlow, qhigh ):
    print 'INFO calculating integration limits [%0.2f,%0.2f] for particle with m,w = (%0.1f, %0.1f)'%( qlow, qhigh, m, w )
    tlow  = math.atan( (pow(qlow,2) - pow(m,2)) / (m * w) )
    thigh = math.atan( (pow(qhigh,2) - pow(m,2)) / (m * w) )
    return (tlow, thigh)
    
mW = ROOT.hepstd.wMass
wW = ROOT.hepstd.wWidth

a,b = 0,0
    
if useBoxPhaseSpace:
    # for low mass Higgs, one W is usually off-shell
    if mH <= (mW + mW): a,b = solveq( mW, wW, 0, mW + nGam*wW )
    else:               a,b = solveq( mW, wW, mW - nGam*wW, mW + nGam*wW )
    mehndl.setIntegrationLimits( mehndl.getIWP(), a, b )
    mehndl.setIntegrationLimits( mehndl.getIWM(), a, b )
    if not useNarrowWidthApprox:
        # integrate +/- 5 GeV or +/- 10 Gamma_{H} (whichever is larger) about the Higgs pole 
        a,b = solveq( mH, wH, mH - max( 5, nGam * wH ), mH + max( 5, nGam * wH ) )
        mehndl.setIntegrationLimits( mehndl.getIH(), a, b )
else:
    if mH < (mW + mW): a,b = 0, mW + nGam*wW 
    else:              a,b = mW - nGam*wW, mW + nGam*wW
    mehndl.setIntegrationLimits( mehndl.getIWP(), a, b )
    mehndl.setIntegrationLimits( mehndl.getIWM(), a, b )
    if not useNarrowWidthApprox:
        # integrate +/- 5 GeV or +/- 10 Gamma_{H} (whichever is larger) about the Higgs pole 
        mehndl.setIntegrationLimits( mehndl.getIH(), mH - max( 5, nGam * wH ), mH + max( 5, nGam * wH ) )

mehndl.setIntegrationLimits( mehndl.getIY(), -3, 3 )

# Set the jet multiplicity bin (determines which ME calculations are setup)

mehndl.setName( 'XWWpBoosted' + str( int(mH) ) + '_0p' )
mehndl.setPDF( 'cteq6m' )

# apply boost correction using recoil estimate
applyBoostCorrection = False

if recoilTF.__class__.__name__ in [ 'SystemBoostHybridTF', 'SystemBoostTF' ]:
    boostCalculator = 'tmvaRecoilEstimator.py'

if selectHelicity:
    mehndl.setIHelicity( iHel )

if useJHUModel:
    mehndl.setUseJHU( True, spin, +1 )
