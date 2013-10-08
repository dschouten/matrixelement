import ROOT

import math

from matrixelement import toolFactory
from matrixelement import getDefaultIntegrator
from matrixelement import getVegasIntegrator
from matrixelement import getAlternateIntegrator
from matrixelement import getHiggsWidth

jetEnergyTF = None
if os.path.exists( 'bjetresolution.tf' ):
    jetEnergyTF = getattr( ROOT, 'JetEnergyResolutionTF' )( 'bjetresolution.tf' ) 
else:
    print 'ERROR [bjetresolution.tf] file not found'
    raise RuntimeError

jetEfficiencyTF = None
if os.path.exists( 'bjetefficiency.tf' ):
    jetEfficiencyTF = getattr( ROOT, 'JetEfficiencyTF' )( 'bjetefficiency.tf' ) 
else:
    print 'ERROR [bjetefficiency.tf] file not found'
    raise RuntimeError

# Configure default/backup/failsafe integrators
myMEIntegrators = { 'default'  : getVegasIntegrator( 0.100, 1.0e5, 2.5e5 ),
                    'backup'   : getVegasIntegrator( 0.150, 5.0e5, 1.0e6 ),
                    'failsafe' : getVegasIntegrator( 0.200, 1.0e6, 2.0e6 ) }

# Configure the ME calculator
myMECalculation = toolFactory( 'TT1j' )( integrator = myMEIntegrators['default'], args = [ ROOT.TT1j.kJET_4D ] )

hndl = myMECalculation.mehndl

hndl.setObservedJetTF( jetEnergyTF )
hndl.setInvisibleJetTF( jetEfficiencyTF )

# Set integration limits (use TF to set dynamic limits for jets)
hndl.setIntegrationLimits( hndl.getIPTB(), 15, 80 )
hndl.setIntegrationLimits( hndl.getIETAB(), -2, 2 )
hndl.setIntegrationLimits( hndl.getIPHIB(), -math.pi, math.pi )

hndl.setPDF( 'ct10' )
