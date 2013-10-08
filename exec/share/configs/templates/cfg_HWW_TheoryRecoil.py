import ROOT

import math

from matrixelement import toolFactory
from matrixelement import getDefaultIntegrator
from matrixelement import getAlternateIntegrator
from matrixelement import getHiggsWidth

# Set the Higgs mass & width for H->WW matrix elements
mH = ??
wH = getHiggsWidth( mH ) + math.pow( getHiggsWidth( mH ) , 1/6.0  ) / 2.0

# Configure default/backup/failsafe integrators
myMEIntegrators = { 'default'  : getDefaultIntegrator( 0.010, 4.0e4, 1.0e6 ),
                    'backup'   : getDefaultIntegrator( 0.025, 4.0e4, 5.0e6 ),
                    'failsafe' : getDefaultIntegrator( 0.100, 4.0e4, 1.0e7 ) }

boostTF = None
if os.path.exists( 'simplerecoil_%d.tf'%( mH ) ):
    boostTF = getattr( ROOT, 'SystemBoostAverageTF' )( 'simplerecoil_%d.tf'%( mH ) ) 
else:
    print 'ERROR [simplerecoil_%d.tf] file not found'%( mH )
    raise RuntimeError

# Configure the ME calculators
myMECalculators = { 'HWW' : toolFactory( 'HWW' )( myMEIntegrators['default'],
                                                  args = [ mH, ROOT.HWW.kNEUTRINO_6D, boostTF ] ) }

# Set integration limits (use TF to set dynamic limits for jets)
myME = myMECalculators['HWW'].matrixElement

myME.setWidth( wH )
myME.setIntegrationLimits( myME.getIPTA(), 0, 250 )
myME.setIntegrationLimits( myME.getIPHIA(), -math.pi, math.pi )
myME.setIntegrationLimits( myME.getIPZA(), -500, 500 )
myME.setIntegrationLimits( myME.getIPZB(), -500, 500 )

# Set the jet multiplicity bin (determines which ME calculations are setup)
myMEJetN        = 0

# Set the ME calculation to perform
myMECalculation = myMECalculators['HWW']

myMECalculation.matrixElement.setName( 'HWWb' + str( int(mH) ) )
