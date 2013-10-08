import ROOT

import os
import math

from matrixelement import toolFactory
from matrixelement import getVegasIntegrator
from matrixelement import getAlternateIntegrator
from matrixelement import getHiggsWidth

timeout = 600

jetEnergyTF = None
if os.path.exists( 'bjetresolution.tf' ): ## old jet resolution
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
myMEIntegrators = { 'default'  : getVegasIntegrator( 0.025, 1.0e4, 1.0e5 ),
                    'backup'   : getVegasIntegrator( 0.050, 1.0e5, 5.0e5 ),
                    'failsafe' : getVegasIntegrator( 0.100, 5.0e5, 1.0e6 ) }

# Configure the ME calculators
myMECalculation = toolFactory( 'TT1jSimple' )( integrator = myMEIntegrators['default'],
                                               args = [ ROOT.TT1jSimple.kJET_8D ] )

myMECalculation.mehndl.setObservedJetTF( jetEnergyTF )
myMECalculation.mehndl.setInvisibleJetTF( jetEfficiencyTF )

# Set integration limits (use TF to set dynamic limits for jets)

m = myMECalculation.mehndl
m.setIntegrationLimits( m.getIPTJB(), 0, 70 )
m.setIntegrationLimits( m.getIETAJB(), -2.5, 2.5 )
m.setIntegrationLimits( m.getIPHIJB(), -math.pi, math.pi )
m.setIntegrationLimits( m.getIPTA(), 0, 150 )
m.setIntegrationLimits( m.getIPHIA(), -math.pi, math.pi )
m.setIntegrationLimits( m.getIPZA(), -200, 200 )
m.setIntegrationLimits( m.getIPZB(), -200, 200 )

# Set the jet multiplicity bin (determines which ME calculations are setup)
