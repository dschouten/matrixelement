import ROOT

import math

from matrixelement import toolFactory
from matrixelement import getLineIntegrator
from matrixelement import getVegasIntegrator

#################################################################################

strategy = ROOT.DY2j.k4D

myMEIntegrators = { 'default'  : getVegasIntegrator( 0.025, 1.0e3, 5.0e3 ),
                    'backup'   : getVegasIntegrator( 0.050, 5.0e3, 1.0e4 ),
                    'failsafe' : getVegasIntegrator( 0.100, 1.0e4, 1.0e5 ) }

# Configure the ME calculator
myMECalculation = toolFactory( 'DY2j' )( integrator = myMEIntegrators['default'], args = [ strategy ] )

jetEnergyTF = None
if os.path.exists( 'gluonjetresolution.tf' ):
    jetEnergyTF = getattr( ROOT, 'JetEnergyResolutionTF' )( 'gluonjetresolution.tf' ) 
else:
    print 'ERROR [gluonjetresolution.tf] file not found'
    raise RuntimeError

myMECalculation.mehndl.setJetTF( jetEnergyTF )

myMECalculation.mehndl.setPDF( 'ct10' )

