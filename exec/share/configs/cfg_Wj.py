import ROOT

import math

from matrixelement import toolFactory
from matrixelement import getDefaultIntegrator
from matrixelement import getVegasIntegrator
from matrixelement import getAlternateIntegrator
from matrixelement import getLineIntegrator
from matrixelement import getHiggsWidth

strategy = ROOT.WFake.k1D

myMEIntegrators = None

# Configure default/backup/failsafe integrators
if strategy == ROOT.WFake.k1D:
    myMEIntegrators = { 'default'  : getLineIntegrator( 0.01, 1.0e3, 5.0e03 ),
                        'backup'   : getLineIntegrator( 0.01, 5.0e3, 1.0e04 ),
                        'failsafe' : getLineIntegrator( 0.01, 1.0e4, 5.0e04 ) }
else:
    myMEIntegrators = { 'default'  : getVegasIntegrator( 0.010, 1.0e3, 2.5e03 ),
                        'backup'   : getVegasIntegrator( 0.025, 2.0e4, 1.0e05 ),
                        'failsafe' : getVegasIntegrator( 0.050, 1.0e5, 5.0e05 ) }

fakeTF = None

if strategy == ROOT.WFake.k2D:
    fakeTF = ROOT.JetEnergyResolutionTF( "fakeresolution.tf" )

# Configure the ME calculators
myMECalculation = toolFactory( 'WFake' )( integrator = myMEIntegrators['default'],
                                          args = [ strategy, fakeTF ] ) 

# Set the jet multiplicity bin (determines which ME calculations are setup)

# Gets integration limits for transformed Q^{2} variable
def solveq( m, w, qlow, qhigh ):
    tlow  = math.atan( (pow(qlow,2) - pow(m,2)) / (m * w) )
    thigh = math.atan( (pow(qhigh,2) - pow(m,2)) / (m * w) )
    return (tlow, thigh)

mW = ROOT.hepstd.wMass
wW = ROOT.hepstd.wWidth

a,b = solveq( mW, wW, 0, 400 )

mehndl = myMECalculation.mehndl

mehndl.setIntegrationLimits( 0, a, b )
mehndl.setName( 'Wj' )
mehndl.setWjetMode()
mehndl.setPDF( 'ct10' )

# Choose to apply boost correction using recoil estimate
applyBoostCorrection = False

allowNull = True

# Correct the parton energy with simple response function
# mehndl.setFakeResponse( '(1 - 4.4 * pow(x, -0.8))' ) ##

myMECalculation.addExec( "if fake.type == 'electron': self.mehndl.setFakeResponse( '0.2' )" )
myMECalculation.addExec( "if fake.type != 'electron': self.mehndl.setFakeResponse( '0.2' )" )

