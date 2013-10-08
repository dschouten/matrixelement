import ROOT

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

recoilTF = None

allowNull            = True   # allow PH == 0, otherwise assume non-convergence and try with finer grid
useNarrowWidthApprox = True   # remove integration over mH
useSMKludge          = False  # flag to use the SM kludge
useJHUModel          = True   # flag to use JHU ME

spin = 0 # specify spin of resonance in pp -> X -> WW(lv,lv) !! JHU only !!
    
# Configure default/backup/failsafe integrators
myMEIntegrators = { 'default'  : getVegasIntegrator( 0.015, 1.0e5, 2.5e5 ),
                    'backup'   : getVegasIntegrator( 0.025, 2.5e5, 5.0e5 ),
                    'failsafe' : getVegasIntegrator( 0.050, 5.0e5, 1.0e6 ) }

# Configure the ME calculators
myMECalculation = toolFactory( 'HWW' )( integrator = myMEIntegrators['default'],
                                        args = ( mH, ROOT.HWW.kWMASS_4D, recoilTF, useNarrowWidthApprox ) )

# Set integration limits (use TF to set dynamic limits for jets)
myMECalculation.mehndl.setWidth( wH )

mehndl = myMECalculation.mehndl

# Gets integration limits for transformed Q^{2} variable
def solveq( m, w, qlow, qhigh ):
    tlow  = math.atan( (pow(qlow,2) - pow(m,2)) / (m * w) )
    thigh = math.atan( (pow(qhigh,2) - pow(m,2)) / (m * w) )
    return (tlow, thigh)

mW = ROOT.hepstd.wMass
wW = ROOT.hepstd.wWidth

# for low mass Higgs, one W is usually off-shell
a,b = 0,0
nGam = 10
if mH < (mW + mW): a,b = solveq( mW, wW, 0, mW + nGam*wW )
else:              a,b = solveq( mW, wW, mW - nGam*wW, mW + nGam*wW )

mehndl.setIntegrationLimits( mehndl.getIWP(), a, b )
mehndl.setIntegrationLimits( mehndl.getIWM(), a, b )

mehndl.setIntegrationLimits( mehndl.getIY(), -3, 3 )

if not useNarrowWidthApprox:
    a,b = solveq( mH, wH, mH - max( 10, nGam * wH ), mH + max( 10, nGam * wH ) )
    mehndl.setIntegrationLimits( mehndl.getIH(),  a, b )

# Use the SM kludge ?
mehndl.setUseSM( useSMKludge )

# Set the jet multiplicity bin (determines which ME calculations are setup)

mehndl.setName( 'HWWp' + str( int(mH) ) + '_0p' )

# apply boost correction using recoil estimate
applyBoostCorrection = False

allowNull = True

cuts = None

if useJHUModel:
    mehndl.setUseJHU( True, spin )
