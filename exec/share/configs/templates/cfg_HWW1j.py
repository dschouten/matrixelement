import ROOT

import math

from matrixelement import toolFactory
from matrixelement import getDefaultIntegrator
from matrixelement import getAlternateIntegrator
from matrixelement import getHiggsWidth

# Set the Higgs mass & width for H->WW matrix elements
mH = ??
wH = getHiggsWidth( mH ) + math.pow( getHiggsWidth( mH ) , 1/6.0  ) / 2.0

jetEnergyTF = None
if os.path.exists( 'gluonjetresolution.tf' ):
    jetEnergyTF = getattr( ROOT, 'JetEnergyResolutionTF' )( 'gluonjetresolution.tf' ) 
else:
    print 'ERROR [gluonjetresolution.tf] file not found'
    raise RuntimeError

# Configure default/backup/failsafe integrators
myMEIntegrators = { 'default'  : getDefaultIntegrator( 0.025, 5.0e4, 5.0e6 ),
                    'backup'   : getDefaultIntegrator( 0.050, 5.0e4, 1.0e7 ),
                    'failsafe' : getDefaultIntegrator( 0.100, 5.0e4, 1.0e7 ) }

# Configure the ME calculators

myMECalculators = { 'HWW1j' : toolFactory( 'HWW1j' )( myMEIntegrators['default'],
                                                      args = [ mH, ROOT.HWW1j.kNEUTRINO_5D, None, jetEnergyTF ] ) }

# myMECalculators = { 'HWW1j' : toolFactory( 'HWW1j' )( myMEIntegrators['default'],
#                                                       args = [ mH, ROOT.HWW1j.kNEUTRINO_4D, None, None ] ) }

# Set integration limits (use TF to set dynamic limits for jets)
myMECalculators['HWW1j'].matrixElement.setWidth( wH )
myMECalculators['HWW1j'].matrixElement.setIntegrationLimits(0, 0, 250)
myMECalculators['HWW1j'].matrixElement.setIntegrationLimits( 1, -math.pi, math.pi )
myMECalculators['HWW1j'].matrixElement.setIntegrationLimits( 2, -500, 500 )
myMECalculators['HWW1j'].matrixElement.setIntegrationLimits( 3, -500, 500 )

# Set the jet multiplicity bin (determines which ME calculations are setup)
myMEJetN        = 1

# Set the ME calculation to perform
myMECalculation = myMECalculators['HWW1j']

myMECalculation.matrixElement.setName( 'HWW1j' + str( mH ) )
