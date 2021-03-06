import ROOT

import math

from matrixelement import toolFactory
from matrixelement import getDefaultIntegrator
from matrixelement import getAlternateIntegrator
from matrixelement import getVegasIntegrator
from matrixelement import getMiserIntegrator

# Configure default/backup/failsafe integrators
myMEIntegrators = { 'default'  : getVegasIntegrator( 0.025, 7.5e4, 1.5e5 ),
                    'backup'   : getVegasIntegrator( 0.050, 1.5e5, 5.0e5 ),
                    'failsafe' : getVegasIntegrator( 0.075, 2.5e5, 1.0e6 ) }

recoilTF = getattr( ROOT, 'SystemBoostTF' )( 'recoil_ww_nj.tf' )

# Configure the ME calculators
myMECalculation = toolFactory( 'WW' )( integrator = myMEIntegrators['default'],
                                       args = [ ROOT.WW.kWMASS_6D, recoilTF ] ) 

mehndl = myMECalculation.mehndl

# Gets integration limits for transformed Q^{2} variable
def solveq( m, w, qlow, qhigh ):
    tlow  = math.atan( (pow(qlow,2) - pow(m,2)) / (m * w) )
    thigh = math.atan( (pow(qhigh,2) - pow(m,2)) / (m * w) )
    return (tlow, thigh)

mW = ROOT.hepstd.wMass
wW = ROOT.hepstd.wWidth

a,b = solveq( mW, wW, 0, 160 )

# Set integration limits (use TF to set dynamic limits for jets)
mehndl.setIntegrationLimits( mehndl.getIWP(), a, b )
mehndl.setIntegrationLimits( mehndl.getIWM(), a, b )

s = 4 * pow(ROOT.hepstd.beamEnergy,2)

tmin = pow(2 * mW,2) / s
tmax = 1 ## 

mehndl.setIntegrationLimits( mehndl.getIT(), tmin, tmax )
mehndl.setIntegrationLimits( mehndl.getIY(), -3, 3 )

# Set the jet multiplicity bin (determines which ME calculations are setup)

mehndl.setName( 'WWpNj' )
mehndl.setPDF( 'ct10' )

# Set to use the strange PDF
myMECalculation.mehndl.setUseStrangePDF( True )

cuts = None

allowNull = True # whether to allow PWW == 0, otherwise assume non-convergence and try with finer grid

# apply boost correction using recoil estimate
boostCalculator = 'tmvaRecoilEstimator.py'
