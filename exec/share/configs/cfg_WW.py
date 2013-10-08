import ROOT

import math

from matrixelement import toolFactory
from matrixelement import getDefaultIntegrator
from matrixelement import getVegasIntegrator
from matrixelement import getAlternateIntegrator
from matrixelement import getHiggsWidth

# Configure default/backup/failsafe integrators
myMEIntegrators = { 'default'  : getDefaultIntegrator( 0.010, 5.0e4, 1.0e06 ),
                    'backup'   : getDefaultIntegrator( 0.025, 5.0e4, 5.0e06 ),
                    'failsafe' : getDefaultIntegrator( 0.075, 5.0e4, 1.0e07 ) }

# Configure the ME calculators
myMECalculation = toolFactory( 'WW' )( integrator = myMEIntegrators['default'],
                                       args = [ ROOT.WW.kNEUTRINO_4D, None ] ) 

# Set integration limits (use TF to set dynamic limits for jets)

myMECalculation.mehndl.setIntegrationLimits( 0, 0, 500 )
myMECalculation.mehndl.setIntegrationLimits( 1, -math.pi, math.pi )
myMECalculation.mehndl.setIntegrationLimits( 2, -500, 500 )
myMECalculation.mehndl.setIntegrationLimits( 3, -500, 500 )

# Set the jet multiplicity bin (determines which ME calculations are setup)

myMECalculation.mehndl.setName( 'WW' )

# Set to use the strange PDF
myMECalculation.mehndl.setUseStrangePDF( True )

cuts = None
