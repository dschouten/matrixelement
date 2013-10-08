import ROOT

import math

from matrixelement import toolFactory
from matrixelement import getDefaultIntegrator
from matrixelement import getAlternateIntegrator
from matrixelement import getHiggsWidth

jetEnergyTF = None
if os.path.exists( 'gluonjetresolution.tf' ):
    jetEnergyTF = getattr( ROOT, 'JetEnergyResolutionTF' )( 'gluonjetresolution.tf' ) 
else:
    print 'ERROR [gluonjetresolution.tf] file not found'
    raise RuntimeError

# Configure default/backup/failsafe integrators
myMEIntegrators = { 'default'  : getDefaultIntegrator( 0.025, 5.0e4, 2.0e06 ),
                    'backup'   : getDefaultIntegrator( 0.050, 5.0e4, 5.0e06 ),
                    'failsafe' : getDefaultIntegrator( 0.100, 5.0e4, 1.0e07 ) }

# Configure the ME calculators

myMECalculation = toolFactory( 'WW1j' )( integrator = myMEIntegrators['default'],
                                         args = [ ROOT.WW1j.kNEUTRINO_5D, jetEnergyTF ] ) 

# Set integration limits (use TF to set dynamic limits for jets)
myMECalculation.mehndl.setIntegrationLimits( 0, 0, 250 )
myMECalculation.mehndl.setIntegrationLimits( 1, -math.pi, math.pi )
myMECalculation.mehndl.setIntegrationLimits( 2, -500, 500 )
myMECalculation.mehndl.setIntegrationLimits( 3, -500, 500 )

# Set the jet multiplicity bin (determines which ME calculations are setup)

# Set to use the strange PDF
myMECalculation.mehndl.setUseStrangePDF( True )

