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
wH = getHiggsWidth( mH ) + math.pow( getHiggsWidth( mH ) , 1/6.0  ) / 2.0
wH = max( 1, getHiggsWidth( mH ) )

recoilTF = getattr( ROOT, 'SystemBoostHybridTF' )( 'recoilhyb_h125.tf' ) 

# Configure default/backup/failsafe integrators
myMEIntegrators = { 'default'  : getDefaultIntegrator( 0.010, 1.0e4, 2.5e6 ),
                    'backup'   : getDefaultIntegrator( 0.025, 1.0e4, 5.0e6 ),
                    'failsafe' : getDefaultIntegrator( 0.100, 1.0e4, 1.0e7 ) }

# Configure the ME calculaskype:bernd.stelzertors
myMECalculation = toolFactory( 'HWW' )( integrator = myMEIntegrators['default'],
                                        args = [ mH, ROOT.HWW.kNEUTRINO_6D, recoilTF ] ) 

# Set integration limits (use TF to set dynamic limits for jets)
mehndl = myMECalculation.mehndl

mehndl.setWidth( wH )
mehndl.setIntegrationLimits( 0, 0, 250 )
mehndl.setIntegrationLimits( 1, -math.pi, math.pi )
mehndl.setIntegrationLimits( 2, -500, 500 )
mehndl.setIntegrationLimits( 3, -500, 500 )

# Use the SM kludge ?
mehndl.setUseSM( False )

# Set the jet multiplicity bin (determines which ME calculations are setup)

myMECalculation.mehndl.setName( 'HWW' + str( int(mH) ) )

# Choose to apply boost correction using recoil estimate
applyBoostCorrection = False

if recoilTF.__class__.__name__ in [ 'SystemBoostHybridTF', 'SystemBoostTF' ]:
    boostCalculator = 'tmvaRecoilEstimator.py'

cuts = None
