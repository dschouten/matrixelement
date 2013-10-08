import ROOT

import math

from matrixelement import toolFactory
from matrixelement import getDefaultIntegrator
from matrixelement import getAlternateIntegrator
from matrixelement import getHiggsWidth
from matrixelement import getVegasIntegrator
from matrixelement import getMiserIntegrator

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
myMEIntegrators = { 'default'  : getVegasIntegrator( 0.050, 1.0e5, 5.0e5 ),
                    'backup'   : getVegasIntegrator( 0.075, 2.5e5, 7.5e5 ),
                    'failsafe' : getVegasIntegrator( 0.100, 7.5e5, 1.5e6 ) }

# Configure the ME calculators
myMECalculators = { 'TT1j' : toolFactory( 'TT1j' )( myMEIntegrators['default'], args = [ ROOT.TT1j.kJET_WMASS_6D ] ) }

myMECalculators['TT1j'].mehndl.setObservedJetTF( jetEnergyTF )
myMECalculators['TT1j'].mehndl.setInvisibleJetTF( jetEfficiencyTF )

# Gets integration limits for transformed Q^{2} variable
def solveq( m, w, qlow, qhigh ):
    tlow  = math.atan( (pow(qlow,2) - pow(m,2)) / (m * w) )
    thigh = math.atan( (pow(qhigh,2) - pow(m,2)) / (m * w) )
    return (tlow, thigh)

mW = ROOT.hepstd.wMass
wW = ROOT.hepstd.wWidth

a,b = solveq( mW, wW, 70, 90 )

# Set integration limits (use TF to set dynamic limits for jets)
mehndl = myMECalculators['TT1j'].mehndl
mehndl.setIntegrationLimits( mehndl.getIPTB(), 10, 80 )
mehndl.setIntegrationLimits( mehndl.getIETAB(), -2.5, 2.5 )
mehndl.setIntegrationLimits( mehndl.getIPHIB(), -math.pi, math.pi )
mehndl.setIntegrationLimits( mehndl.getIWP(), a, b )
mehndl.setIntegrationLimits( mehndl.getIWM(), a, b )

# Set the jet multiplicity bin (determines which ME calculations are setup)

# Set the ME calculation to perform
myMECalculation = myMECalculators['TT1j']

mehndl.setPDF( 'cteq6m' )
