import ROOT

import math

from matrixelement import toolFactory
from matrixelement import getVegasIntegrator

jetEnergyTF = None
if os.path.exists( 'gluonjetresolution.tf' ):
    jetEnergyTF = getattr( ROOT, 'JetEnergyResolutionTF' )( 'gluonjetresolution.tf' ) 
else:
    print 'ERROR [gluonjetresolution.tf] file not found'
    raise RuntimeError

# Configure default/backup/failsafe integrators
myMEIntegrators = { 'default'  : getVegasIntegrator( 0.025, 5.0e4, 1.0e05 ),
                    'backup'   : getVegasIntegrator( 0.050, 1.0e5, 5.0e05 ),
                    'failsafe' : getVegasIntegrator( 0.100, 5.0e5, 1.0e06 ) }

# Configure the ME calculators

myMECalculation = toolFactory( 'WW1j' )( integrator = myMEIntegrators['default'],
                                         args = [ ROOT.WW1j.kWMASS_4D, jetEnergyTF ] ) 

mehndl = myMECalculation.mehndl

# Gets integration limits for transformed Q^{2} variable
def solveq( m, w, qlow, qhigh ):
    print 'INFO calculating integration limits [%0.2f,%0.2f] for particle with m,w = (%0.1f, %0.1f)'%( qlow, qhigh, m, w )
    tlow  = math.atan( (pow(qlow,2) - pow(m,2)) / (m * w) )
    thigh = math.atan( (pow(qhigh,2) - pow(m,2)) / (m * w) )
    return (tlow, thigh)
    
mW = ROOT.hepstd.wMass
wW = ROOT.hepstd.wWidth

nGam = 10

s = 4 * pow(ROOT.hepstd.beamEnergy,2)

tmin = 0 ## pow(2 * mW,2) / s
tmax = 1 ## 

a,b = solveq( mW, wW, mW - nGam*wW, mW + nGam*wW )

# Set integration limits (use TF to set dynamic limits for jets)
mehndl.setIntegrationLimits( mehndl.getIT(), tmin, tmax )
mehndl.setIntegrationLimits( mehndl.getIY(), -3, 3 )
mehndl.setIntegrationLimits( mehndl.getIWP(), a, b ) 
mehndl.setIntegrationLimits( mehndl.getIWM(), a, b )

# Set the jet multiplicity bin (determines which ME calculations are setup)

# Set to use the strange PDF
mehndl.setUseStrangePDF( False )

mehndl.setPDF( 'ct10' )

cuts = None

myMECalculation.mehndl.setName( 'WWp1j' )
