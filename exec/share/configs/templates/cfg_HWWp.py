import ROOT

import math

from matrixelement import toolFactory
from matrixelement import getDefaultIntegrator
from matrixelement import getAlternateIntegrator
from matrixelement import getVegasIntegrator
from matrixelement import getHiggsWidth

# Set the Higgs mass & width for H->WW matrix elements
mH = ??
wH = getHiggsWidth( mH ) ## + math.pow( getHiggsWidth( mH ) , 1/6.0  ) / 2.0

recoilTF = None ## = getattr( ROOT, 'SystemBoostAverageTF' )( 'recoilavg_ggh125.tf' ) 

# Configure default/backup/failsafe integrators
myMEIntegrators = { 'default'  : getVegasIntegrator( 0.010, 2.0e5, 5.0e6 ),
                    'backup'   : getVegasIntegrator( 0.025, 2.0e5, 1.0e7 ),
                    'failsafe' : getVegasIntegrator( 0.100, 2.0e5, 5.0e7 ) }

# Configure the ME calculators
myMECalculators = { 'HWW' : toolFactory( 'HWW' )( myMEIntegrators['default'],
                                                  args = [ mH, ROOT.HWW.kWMASS_4D, recoilTF ] ) }

myMECalculators['HWW'].matrixElement.setUseBoostGammaXY( False )

# Set integration limits (use TF to set dynamic limits for jets)
myMECalculators['HWW'].matrixElement.setWidth( wH )

# Gets integration limits for transformed Q^{2} variable
def solveq( m, w, qlow, qhigh ):
    tlow  = math.atan( (pow(qlow,2) - pow(m,2)) / (m * w) )
    thigh = math.atan( (pow(qhigh,2) - pow(m,2)) / (m * w) )
    return (tlow, thigh)

db = ROOT.TDatabasePDG()

mW = db.GetParticle( 'W+' ).Mass() 
wW = db.GetParticle( 'W+' ).Width()

a,b = solveq( mW, wW, mH - mW - wW/0.1, mW + wW/0.1 )

theME = myMECalculators['HWW'].matrixElement

theME.setIntegrationLimits( theME.getIWP(), a, b )
theME.setIntegrationLimits( theME.getIWP(), a, b )

a,b = solveq( mH, wH, mH - wH/0.1, mH + wH/0.1 )
theME.setIntegrationLimits( theME.getIH(),  a, b )
theME.setIntegrationLimits( theME.getIY(), -3, 3 )

# Set the jet multiplicity bin (determines which ME calculations are setup)
myMEJetN        = 0

# Set the ME calculation to perform
myMECalculation = myMECalculators['HWW']

myMECalculation.matrixElement.setName( 'HWWp' + str( int(mH) ) )
