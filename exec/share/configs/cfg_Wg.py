import ROOT

import math

from matrixelement import toolFactory
from matrixelement import getDefaultIntegrator
from matrixelement import getVegasIntegrator
from matrixelement import getAlternateIntegrator
from matrixelement import getLineIntegrator
from matrixelement import getHiggsWidth

strategy = ROOT.WFake.k2D

myMEIntegrators = None

# Configure default/backup/failsafe integrators
if strategy == ROOT.WFake.k1D:
    myMEIntegrators = { 'default'  : getLineIntegrator( 0.001, 1.0e3, 1.0e06 ),
                        'backup'   : getLineIntegrator( 0.010, 1.0e3, 5.0e06 ),
                        'failsafe' : getLineIntegrator( 0.020, 1.0e3, 1.0e07 ) }
else:
    myMEIntegrators = { 'default'  : getVegasIntegrator( 0.010, 5.0e3, 1.0e04 ),
                        'backup'   : getVegasIntegrator( 0.025, 1.0e4, 5.0e04 ),
                        'failsafe' : getVegasIntegrator( 0.050, 5.0e4, 1.0e06 ) }

convTF = ROOT.PhotonConversionTF( "conversion.tf" )

# Configure the ME calculators
myMECalculation = toolFactory( 'WFake' )( integrator = myMEIntegrators['default'],
                                          args = [ strategy, convTF ] )

# Set the jet multiplicity bin (determines which ME calculations are setup)

# Gets integration limits for transformed Q^{2} variable
def solveq( m, w, qlow, qhigh ):
    tlow  = math.atan( (pow(qlow,2) - pow(m,2)) / (m * w) )
    thigh = math.atan( (pow(qhigh,2) - pow(m,2)) / (m * w) )
    return (tlow, thigh)

mW = ROOT.hepstd.wMass
wW = ROOT.hepstd.wWidth

a,b = solveq( mW, wW, 0, 300 )

myMECalculation.mehndl.setIntegrationLimits( 0, a, b )

if strategy == ROOT.WFake.k2D:
    myMECalculation.mehndl.setIntegrationLimits( 1, 0, 1 )

myMECalculation.mehndl.setName( 'Wg' )
myMECalculation.mehndl.setWgammaMode( )

# Choose to apply boost correction using recoil estimate
applyBoostCorrection = False

cuts = None

allowNull = True
