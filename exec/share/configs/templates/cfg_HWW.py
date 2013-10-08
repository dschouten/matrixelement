import ROOT

import math

from matrixelement import toolFactory
from matrixelement import getDefaultIntegrator
from matrixelement import getAlternateIntegrator
from matrixelement import getHiggsWidth

# Set the Higgs mass & width for H->WW matrix elements
mH = ??
wH = getHiggsWidth( mH ) + math.pow( getHiggsWidth( mH ) , 1/6.0  ) / 2.0

# Configure default/backup/failsafe integrators
myMEIntegrators = { 'default'  : getDefaultIntegrator( 0.010, 2.0e4, 2.5e6 ),
                    'backup'   : getDefaultIntegrator( 0.025, 2.0e4, 5.0e6 ),
                    'failsafe' : getDefaultIntegrator( 0.100, 2.0e4, 1.0e7 ) }

# Configure the ME calculators
myMECalculators = { 'HWW' : toolFactory( 'HWW' )( myMEIntegrators['default'],
                                                  args = [ mH, ROOT.HWW.kNEUTRINO_4D, None ] ) }

# Set integration limits (use TF to set dynamic limits for jets)
myMECalculators['HWW'].matrixElement.setWidth( wH )
myMECalculators['HWW'].matrixElement.setIntegrationLimits( 0, 0, 250 )
myMECalculators['HWW'].matrixElement.setIntegrationLimits( 1, -math.pi, math.pi )
myMECalculators['HWW'].matrixElement.setIntegrationLimits( 2, -500, 500 )
myMECalculators['HWW'].matrixElement.setIntegrationLimits( 3, -500, 500 )

# Set the jet multiplicity bin (determines which ME calculations are setup)
myMEJetN        = 0

# Set the ME calculation to perform
myMECalculation = myMECalculators['HWW']

myMECalculation.matrixElement.setName( 'HWW' + str( int(mH) ) )
