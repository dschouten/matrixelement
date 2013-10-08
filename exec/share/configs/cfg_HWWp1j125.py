import ROOT

import math

from matrixelement import toolFactory
from matrixelement import getHiggsWidth
from matrixelement import getVegasIntegrator

# Set the Higgs mass & width for H->WW matrix elements
mH = 125
wH = getHiggsWidth( mH )

jetEnergyTF = None
if os.path.exists( 'gluonjetresolution.tf' ):
    jetEnergyTF = getattr( ROOT, 'JetEnergyResolutionTF' )( 'gluonjetresolution.tf' ) 
else:
    print 'ERROR [gluonjetresolution.tf] file not found'
    raise RuntimeError

useNarrowWidth = True

# Configure default/backup/failsafe integrators
myMEIntegrators = { 'default'  : getVegasIntegrator( 0.015, 2.5e5, 5.0e5 ),
                    'backup'   : getVegasIntegrator( 0.025, 5.0e5, 1.0e6 ),
                    'failsafe' : getVegasIntegrator( 0.050, 1.0e6, 2.0e6 ) }

# Configure the ME calculators

myMECalculator = toolFactory( 'HWW1j' )( myMEIntegrators['default'],
                                         args = [ mH, ROOT.HWW1j.kWMASS_5D, jetEnergyTF, useNarrowWidth ] ) 

mehndl = myMECalculator.mehndl

# Gets integration limits for transformed Q^{2} variable
def solveq( m, w, qlow, qhigh ):
    print 'INFO calculating integration limits [%0.2f,%0.2f] for particle with m,w = (%0.1f, %0.1f)'%( qlow, qhigh, m, w )
    tlow  = math.atan( (pow(qlow,2) - pow(m,2)) / (m * w) )
    thigh = math.atan( (pow(qhigh,2) - pow(m,2)) / (m * w) )
    return (tlow, thigh)
    
mW = ROOT.hepstd.wMass
wW = ROOT.hepstd.wWidth

# for low mass Higgs, one W is usually off-shell
a,b = 0,0
nGam = 10
if mH < (mW + mW): a,b = solveq( mW, wW, 0, mW + nGam*wW )
else:              a,b = solveq( mW, wW, mW - nGam*wW, mW + nGam*wW )

mehndl.setIntegrationLimits( mehndl.getIWP(), a, b )
mehndl.setIntegrationLimits( mehndl.getIWM(), a, b )

if not useNarrowWidth:
    a,b = solveq( mH, wH, mH - max( 10, nGam*wH ), mH + max( 10, nGam*wH ) )
    mehndl.setIntegrationLimits( mehndl.getIH(),  a, b )
    
mehndl.setIntegrationLimits( mehndl.getIY(), -3, 3 )

# Set the jet multiplicity bin (determines which ME calculations are setup)

# Set the ME calculation to perform
myMECalculation = myMECalculator

myMECalculation.mehndl.setName( 'HWWp1j' + str( mH ) )

mehndl.setPDF( 'ct10' )

# Set to use the strange PDF
mehndl.setUseStrangePDF( False )

cuts = None

allowNull = True 
