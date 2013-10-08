import ROOT
import os
import sys
import math

from matrixelement import toolFactory
from matrixelement import getDefaultIntegrator
from matrixelement import getAlternateIntegrator
from matrixelement import getVegasIntegrator
from matrixelement import getMiserIntegrator
from matrixelement import getHiggsWidth

# Set the Higgs mass & width for H->WW matrix elements
mH = 125
wH = getHiggsWidth( mH )

recoilTF = None

allowNull            = True   # allow PH == 0, otherwise assume non-convergence and try with finer grid
useNarrowWidthApprox = True   # remove integration over mH
useBoxPhaseSpace     = True   # if true, use q -> atan(t) transformation to sample BW
useRSModel           = False  # flag to use RS graviton model (spin 2)
useJHUModel          = True   # flag to use JHU model
selectHelicity       = False  # flag to select particular helicity combination (has no effect if not implemented in ME)

iHel = -1
if os.environ.has_key( 'MYOPTS' ):
    opts = os.environ['MYOPTS']
    if 'iHel' in opts: # integrate only using this helicity combination (index starts @ 1)
        iHel = int( filter( lambda t: 'iHel' in t, opts.split( ':' ) )[0].replace( 'iHel=', '' ) )
        selectHelicity = True
        print 'INFO using helicity i =', iHel
else:
    selectHelicity = False
    iHel = -1

# Configure default/backup/failsafe integrators
myMEIntegrators = { 'default'  : getVegasIntegrator( 0.015, 1.00e4, 5.00e4 ),
                    'backup'   : getVegasIntegrator( 0.050, 5.00e4, 1.00e5 ),
                    'failsafe' : getVegasIntegrator( 0.075, 1.00e5, 5.00e5 ) }

# Configure the ME calculators
myMECalculation = toolFactory( 'HWW' )( integrator = myMEIntegrators['default'],
                                        args = ( mH, ROOT.HWW.kWMASS_4D, recoilTF, useNarrowWidthApprox ) )

# Set integration limits (use TF to set dynamic limits for jets)
myMECalculation.mehndl.setWidth( wH )

mehndl = myMECalculation.mehndl

# Gets integration limits for transformed Q^{2} variable
def solveq( m, w, qlow, qhigh ):
    tlow  = math.atan( (pow(qlow,2) - pow(m,2)) / (m * w) )
    thigh = math.atan( (pow(qhigh,2) - pow(m,2)) / (m * w) )
    return (tlow, thigh)

mW = ROOT.hepstd.wMass
wW = ROOT.hepstd.wWidth

a,b = 0,0
nGam = 10

# for low mass Higgs, one W is usually off-shell
if mH < (mW + mW): a,b = solveq( mW, wW, 0, mW + nGam*wW )
else:              a,b = solveq( mW, wW, mW - nGam*wW, mW + nGam*wW )
mehndl.setIntegrationLimits( mehndl.getIWP(), a, b )
mehndl.setIntegrationLimits( mehndl.getIWM(), a, b )
if not useNarrowWidthApprox:
    a,b = solveq( mH, wH, mH - max( 10, nGam * wH ), mH + max( 10, nGam * wH ) )
    mehndl.setIntegrationLimits( mehndl.getIH(),  a, b )

mehndl.setIntegrationLimits( mehndl.getIY(), -3, 3 )

mehndl.setUseSM( False )
mehndl.setUseRS( useRSModel )

mehndl.setUseBoxPS( useBoxPhaseSpace )

# Set the jet multiplicity bin (determines which ME calculations are setup)

mehndl.setName( 'GWWp' + str( int(mH) ) )

if selectHelicity and iHel >= 0:
    mehndl.setName( 'GWWp%d_i%d'%( int(mH), iHel ) )

# apply boost correction using recoil estimate
applyBoostCorrection = False

allowNull = True

cuts = None

if selectHelicity:
    mehndl.setIHelicity( iHel )

if useJHUModel:
    if os.environ.has_key( 'MYOPTS' ):
        opts = os.environ['MYOPTS']
        if 'qqFraction' in opts:
            qqFraction = float( filter( lambda t: 'qqFraction' in t, opts.split( ':' ) )[0].replace( 'qqFraction=', '' ) )
    else:
        qqFraction = 0
        
    mehndl.setUseJHU( True, 2 )
    mehndl.setJHUFractionQQ( qqFraction )

    n = mehndl.name()
    mehndl.setName( n + '_%02dqq'%( int( round( qqFraction * 100 ) ) ) )

mehndl.setPDF( 'ct10' )
