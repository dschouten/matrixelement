import ROOT

import math

from matrixelement import toolFactory
from matrixelement import getLineIntegrator
from matrixelement import getVegasIntegrator

useTauDecay = True
maxRapidity = 5

#################################################################################

strategy = ROOT.DY.kJET_1D

# Configure default/backup/failsafe integrators

myMEIntegrators = { 'default'  : getLineIntegrator( 0.001, 1.0e3, 5.0e05 ),                
                    'backup'   : getLineIntegrator( 0.010, 1.0e3, 1.0e06 ),
                    'failsafe' : getLineIntegrator( 0.020, 1.0e3, 5.0e06 ) }

if useTauDecay:
    myMEIntegrators = { 'default'  : getVegasIntegrator( 0.01, 5.0e3, 2.5e04 ),
                        'backup'   : getVegasIntegrator( 0.05, 2.5e4, 1.0e05 ),
                        'failsafe' : getVegasIntegrator( 0.10, 1.0e5, 5.0e05 ) }
    
    strategy = ROOT.DY.kJET_TT_3D

# Configure the ME calculator
myMECalculation = toolFactory( 'DY' )( integrator = myMEIntegrators['default'], args = [ strategy ] )

# Set integration limits (use TF to set dynamic limits for jets)
myMECalculation.mehndl.setIntegrationLimits( 0, -maxRapidity, maxRapidity )

myMECalculation.mehndl.setPDF( 'ct10' )

