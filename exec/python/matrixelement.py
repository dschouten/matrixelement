
#_________________________________________________________________
import ROOT as root

import os
import sys
import msg
import math
import random
import array

try:
    import numpy
except:
    numpy = None

from array import array
from ROOT import gSystem
from ROOT import TLorentzVector

__log = msg.msglog( 'matrixelement', 'warning' )
__log.setPrintOnly( False )

__wmap = { }

from math import sqrt

tlv = TLorentzVector

#_________________________________________________________________
def loadModule( name ):
    stat = gSystem.Load( name )
    if stat < 0:
        __log.error( 'could not load [' + name + ']' )
        raise RuntimeError

# load all C++ tools and namespaces
try:
    # Standard Makefile build system
    loadModule( 'libdhelas.so' )
    loadModule( 'libttanalysis.so' )
    loadModule( 'libintegrator.so' )
    loadModule( 'libintegrator_dict.so' )
    loadModule( 'libmatrix.so' )
    loadModule( 'libmatrix_dict.so' )
except RuntimeError:
    __log.fatal( 'could not load necessary modules' )

#_________________________________________________________________
def fptr( ):
    return array( 'd', [0.] )

#_________________________________________________________________
def iptr( ):
    return array( 'i', [0] )

#_________________________________________________________________
class Divonne( root.CubaIntegrators.Divonne ):
    pass

#_________________________________________________________________
class Vegas( root.CubaIntegrators.Vegas ):
    pass

#_________________________________________________________________
class Cuhre( root.CubaIntegrators.Cuhre ):
    pass

#_________________________________________________________________
class GSLLine( root.GSLIntegrator1D ):
    pass

#_________________________________________________________________
class GSLMiser( root.GSLIntegratorMiser ):
    pass

#_________________________________________________________________
class GSLVegas( root.GSLIntegratorVegas ):
    pass

#_________________________________________________________________
class NullIntegrator( root.NullIntegrator ):
    pass
  
#_________________________________________________________________
myKnownMatrixElements = [ 'HWW' ,
                          'WW',
                          'TT1j',
                          'TT1jSimple',
                          'HWW1j',
                          'WW1j',
                          'DY',
                          'DY2j',
                          'WFake' ]

#_________________________________________________________________
class MatrixElementToolBase:
    def __init__( self,
                  name = '',
                  integrator = GSLVegas(),
                  meklass = '',
                  args = () ):
        self.log = msg.msglog( name, 'warning' )
        self.integrator = integrator
        self.meklass = meklass
        self.results = None
        self.__ctor_args = args
        self.__loadMatrixElement( )
        
    def __loadMatrixElement( self ):
        if self.meklass not in myKnownMatrixElements:
            self.log.warning( 'Matrix element [', self.meklass, '] is unknown' )
        try:
            self.mehndl = getattr( root, self.meklass )( self.integrator, *self.__ctor_args )
        except AttributeError:
            self.log.error( 'No matrix element [', self.meklass, '] available' )
            self.mehndl = None
        
    def setIntegrator( self, integrator ):
        self.integrator = integrator
        self.mehndl.setIntegrator( integrator )

    def getIntegrator( self ):
        return self.integrator
        
    def setMatrixElement( self, meklass ):
        self.meklass = meklass
        self.__loadMatrixElement( )
        
    def getMatrixElement( self ):
        return self.meklass

    def __call__( self, *args, **kwargs ):
        raise NotImplemented

    def evaluate( self, *args, **kwargs ):
        raise NotImplemented

#_________________________________________________________________
class WWMatrixElementToolBase( MatrixElementToolBase ):
    def __init__( self, meklass, integrator = GSLVegas(), args = () ):
        MatrixElementToolBase.__init__( self, 'WWMatrixElementTool', integrator, meklass, args )
        pass
    
    def __call__( self, lplus, lminus, jets, missingET, nu = TLorentzVector(), nubar = TLorentzVector(),
                  doIntegration = True, coords = None, **kwargs ):
        self.mehndl.measured.lm = lminus
        self.mehndl.measured.lp = lplus
        self.mehndl.measured.nul = nu
        self.mehndl.measured.nur = nubar
        self.mehndl.measured.met = missingET
        if doIntegration:
            self.mehndl.initialize( )
            retval, error, prob = fptr(), fptr(), fptr()
            fail, neval = iptr(), iptr()
            self.mehndl.integrator().doIntegral( retval, error, fail, neval, prob )
            self.mehndl.finalize()
            return {
                'me'     : retval[0],
                'error'  : error[0],
                'status' : fail[0],
                'neval'  : neval[0],
                'prob'   : prob[0]  }
        else:
            if coords == None:
                self.log.error( 'integration mode disabled, but no phase space coordinates given' )
                raise RuntimeError
            self.mehndl.initialize( )
            self.results = [ self.mehndl( array( 'd', coords ) ) ]
            self.mehndl.finalize()
            return { 'me' : self.results[ 0 ], 
                     'error' : 0, 'status' : 0, 'neval' : 0, 'prob' : 0 }
        return None

    def evaluate( self, lplus, lminus, jets, missingET, nu = TLorentzVector(), nubar = TLorentzVector(),
                  doIntegration = True, coords = None, **kwargs ):
        return self.__call__( lplus, lminus, jets, missingET, nu, nubar, doIntegration, coords )

#_________________________________________________________________
class HWWMatrixElementTool( WWMatrixElementToolBase ):
    def __init__( self, args = (), integrator = GSLVegas() ):
        WWMatrixElementToolBase.__init__( self, 'HWW', integrator, args )
        pass

#_________________________________________________________________
class WWMatrixElementTool( WWMatrixElementToolBase ):
    def __init__( self, integrator = GSLVegas(), args = () ):
        WWMatrixElementToolBase.__init__( self, 'WW', integrator, args )
        pass

#_________________________________________________________________
class tTMatrixElementTool( MatrixElementToolBase ):
    def __init__( self, integrator = GSLVegas(), args = () ):
        MatrixElementToolBase.__init__( self, 'tTMatrixElementTool', integrator, 'TT1j', args )
        pass
        
    def __call__( self, lplus, lminus, jets, missingET, nu = TLorentzVector(), nubar = TLorentzVector(),
                  doIntegration = True, coords = None, **kwargs ):     
        self.mehndl.measured.lm = lminus
        self.mehndl.measured.lp = lplus
        self.mehndl.measured.met = missingET
        if len( jets ) > 0: self.mehndl.measured.jeta = jets[0]
        if len( jets ) > 1: self.mehndl.measured.jetb = jets[1]
        self.mehndl.measured.nul = nu
        self.mehndl.measured.nur = nubar
        if doIntegration:
            self.mehndl.initialize( )
            retval, error, prob = fptr(), fptr(), fptr()
            fail, neval = iptr(), iptr()
            self.mehndl.integrator().doIntegral( retval, error, fail, neval, prob )
            self.mehndl.finalize()
            return {
                'me'     : retval[0],
                'error'  : error[0],
                'status' : fail[0],
                'neval'  : neval[0],
                'prob'   : prob[0]  }
        else:
            if coords == None:
                self.log.error( 'integration mode disabled, but no phase space coordinates given' )
                raise RuntimeError
            self.mehndl.initialize( )
            self.results = [ self.mehndl( array( 'd', coords ) ) ]
            self.mehndl.finalize()
            return { 'me' : self.results[ 0 ], 
                     'error' : 0, 'status' : 0, 'neval' : 0, 'prob' : 0 }
        return None

    def evaluate( self, lplus, lminus, jets, missingET, nu = TLorentzVector(), nubar = TLorentzVector(),
                  doIntegration = True, coords = None, **kwargs ):
        return self.__call__( lplus, lminus, jets, missingET, nu, nubar, doIntegration, coords )

#_________________________________________________________________
class tTSimpleMatrixElementTool( tTMatrixElementTool ):
    def __init__( self, integrator = GSLVegas(), args = () ):
        MatrixElementToolBase.__init__( self, 'tTMatrixElementTool', integrator, 'TT1jSimple', args )
        pass

#_________________________________________________________________
class WW1jMatrixElementToolBase( MatrixElementToolBase ):
    def __init__( self, meklass, integrator = GSLVegas(), args = () ):
        MatrixElementToolBase.__init__( self, 'WW1jMatrixElementTool', integrator, meklass, args )
        pass
    
    def __call__( self, lplus, lminus, jets, missingET, nu = TLorentzVector(), nubar = TLorentzVector(),
                  doIntegration = True, coords = None, **kwargs ):     
        self.mehndl.measured.lm = lminus
        self.mehndl.measured.lp = lplus
        self.mehndl.measured.nul = nu
        self.mehndl.measured.nur = nubar
        self.mehndl.measured.met = missingET
        self.mehndl.measured.jet = jets[0]
        if doIntegration:
            self.mehndl.initialize( )
            retval, error, prob = fptr(), fptr(), fptr()
            fail, neval = iptr(), iptr()
            self.mehndl.integrator().doIntegral( retval, error, fail, neval, prob )
            self.mehndl.finalize()
            return {
                'me'     : retval[0],
                'error'  : error[0],
                'status' : fail[0],
                'neval'  : neval[0],
                'prob'   : prob[0]  }
        else:
            if coords == None:
                self.log.error( 'integration mode disabled, but no phase space coordinates given' )
                raise RuntimeError
            self.mehndl.initialize( )
            self.results = [ self.mehndl( array( 'd', coords ) ) ]
            self.mehndl.finalize()
            return { 'me' : self.results[ 0 ], 
                     'error' : 0, 'status' : 0, 'neval' : 0, 'prob' : 0 }
        return None

    def evaluate( self, lplus, lminus, jets, missingET, nu = TLorentzVector(), nubar = TLorentzVector(),
                  doIntegration = True, coords = None, **kwargs ):
        return self.__call__( lplus, lminus, jets, missingET, nu, nubar, doIntegration, coords )

#_________________________________________________________________
class HWW1jMatrixElementTool( WW1jMatrixElementToolBase ):
    def __init__( self, integrator = GSLVegas(), args = () ):
        WW1jMatrixElementToolBase.__init__( self, 'HWW1j', integrator, args )
        pass

#_________________________________________________________________
class WW1jMatrixElementTool( WW1jMatrixElementToolBase ):
    def __init__( self, integrator = GSLVegas(), args = () ):
        WW1jMatrixElementToolBase.__init__( self, 'WW1j', integrator, args )
        pass

#_________________________________________________________________
class DYMatrixElementTool( MatrixElementToolBase ):
    def __init__( self, integrator = GSLVegas(), args = () ):
        MatrixElementToolBase.__init__( self, 'DYMatrixElementTool', integrator, 'DY', args )
        pass
    
    def __call__( self, lplus, lminus, jets, missingET, doIntegration = True, coords = None, **kwargs ):     
        self.mehndl.measured.lm  = lminus
        self.mehndl.measured.lp  = lplus
        self.mehndl.measured.jet = -(lminus + lplus) 
        if doIntegration:
            self.mehndl.initialize( )
            retval, error, prob = fptr(), fptr(), fptr()
            fail, neval = iptr(), iptr()
            self.mehndl.integrator().doIntegral( retval, error, fail, neval, prob )
            self.mehndl.finalize()
            return {
                'me'     : retval[0],
                'error'  : error[0],
                'status' : fail[0],
                'neval'  : neval[0],
                'prob'   : prob[0]  }
        else:
            if coords == None:
                self.log.error( 'integration mode disabled, but no phase space coordinates given' )
                raise RuntimeError
            self.mehndl.initialize( )
            self.results = [ self.mehndl( array( 'd', coords ) ) ]
            self.mehndl.finalize()
            return { 'me' : self.results[ 0 ], 
                     'error' : 0, 'status' : 0, 'neval' : 0, 'prob' : 0 }
        return None

    def evaluate( self, lplus, lminus, jets, missingET,
                  doIntegration = True, coords = None, **kwargs ):
        return self.__call__( lplus, lminus, missingET, doIntegration, coords )
    
#_________________________________________________________________
class DY2jMatrixElementTool( MatrixElementToolBase ):
    def __init__( self, integrator = GSLVegas(), args = () ):
        MatrixElementToolBase.__init__( self, 'DY2jMatrixElementTool', integrator, 'DY2j', args )
        pass
    
    def __call__( self, lplus, lminus, jets, missingET, doIntegration = True, coords = None, **kwargs ):     
        self.mehndl.measured.lm  = lminus
        self.mehndl.measured.lp  = lplus
        self.mehndl.measured.jeta = jets[0]
        self.mehndl.measured.jetb = jets[1]
        if doIntegration:
            self.mehndl.initialize( )
            retval, error, prob = fptr(), fptr(), fptr()
            fail, neval = iptr(), iptr()
            self.mehndl.integrator().doIntegral( retval, error, fail, neval, prob )
            self.mehndl.finalize()
            return { 'me'     : retval[0],
                     'error'  : error[0],
                     'status' : fail[0],
                     'neval'  : neval[0],
                     'prob'   : prob[0]  }
        else:
            if coords == None:
                self.log.error( 'integration mode disabled, but no phase space coordinates given' )
                raise RuntimeError
            self.mehndl.initialize( )
            self.results = [ self.mehndl( array( 'd', coords ) ) ]
            self.mehndl.finalize()
            return { 'me' : self.results[ 0 ], 
                     'error' : 0, 'status' : 0, 'neval' : 0, 'prob' : 0 }
        return None

    def evaluate( self, lplus, lminus, jets, missingET,
                  doIntegration = True, coords = None, **kwargs ):
        return self.__call__( lplus, lminus, jets, missingET, doIntegration, coords )
    
#_________________________________________________________________
class WFakeMatrixElementTool( MatrixElementToolBase ):
    def __init__( self, integrator = GSLVegas(), args = () ):
        MatrixElementToolBase.__init__( self, 'WFakeMatrixElementTool', integrator, 'WFake', args )
        self.__exec_lines = []
        self.__lepton_decay = None
        pass

    def addExec( self, l ):
        self.__exec_lines.append( l )

    def forceLeptonDecay( self, d ):
        self.__lepton_decay = d
    
    def __call__( self, lp, lm, jets, met, doIntegration = True, coords = None, **kwargs ):
        real, fake = None, None
        if self.mehndl.isWgammaMode():
            # logic: >= 1 lepton(s) has to be an electron to be considered a fake lepton candidate event
            #        if both are electrons, the lowest pT lepton is set as the fake
            if lp.type == 'electron':
                if lm.type == 'electron':
                    if lp.Pt() > lm.Pt():
                        real, fake = lp, lm
                    else:
                        real, fake = lm, lp
                else:
                    real, fake = lm, lp
            elif lm.type == 'electron':
                real, fake = lp, lm
        elif self.mehndl.isWjetMode():
            # logic: just set sub-leading leptons as the fake unless
            # a specific topology is specified with forceLeptonDecay( 'real:fake' )
            if self.__lepton_decay == None \
                   or lp.type == lm.type   \
                   or self.__lepton_decay.lower() in ['e:e', 'm:m']:
                # logic: set the sub-leading lepton as the fake 
                if lp.Pt() > lm.Pt():
                    real, fake = lp, lm
                if lm.Pt() > lp.Pt():
                    real, fake = lm, lp
            else:
                if self.__lepton_decay.lower() == 'e:m' and lp.type != lm.type:
                    # the fake is muon, real is electron
                    if lp.type == 'electron': ## << assume DF here
                        real, fake = lp, lm
                    else:
                        real, fake = lm, lp
                elif self.__lepton_decay.lower() == 'm:e' and lp.type != lm.type:
                    # the fake is electron, real is muon
                    if lp.type == 'muon': ## << assume DF here
                        real, fake = lp, lm
                    else:
                        real, fake = lm, lp
            
        if real == None or fake == None: # no fake leptons, so not a fake candidate
            return {
                'me'     : 0,    # << ME is 0 since object ID means can't be a W+fake event
                'error'  : 0,    # error is +/- 0
                'status' : 0,    # set status == 0 to indicate that this is not a failed convergence
                'neval'  : 0,    # no ME evaluations ...
                'prob'   : -1 }  # flag this event with probability == -1 to indicate ME is set 'by hand'
        self.mehndl.measured.real  = real
        self.mehndl.measured.fake  = fake
        self.mehndl.measured.nu    = tlv()
        self.mehndl.measured.fakecharge = fake.charge
        self.mehndl.measured.realcharge = real.charge

        for l in self.__exec_lines:
            exec l
        
        if doIntegration:
            self.mehndl.initialize( )
            retval, error, prob = fptr(), fptr(), fptr()
            fail, neval = iptr(), iptr()
            self.mehndl.integrator().doIntegral( retval, error, fail, neval, prob )
            self.mehndl.finalize()
            return {
                'me'     : retval[0],
                'error'  : error[0],
                'status' : fail[0],
                'neval'  : neval[0],
                'prob'   : prob[0]  }
        else:
            if coords == None:
                self.log.error( 'integration mode disabled, but no phase space coordinates given' )
                raise RuntimeError
            self.mehndl.initialize( )
            self.results = [ self.mehndl( array( 'd', coords ) ) ]
            self.mehndl.finalize()
            return { 'me' : self.results[0],
                     'error' : 0, 'status' : 0, 'neval' : 0, 'prob' : 0 }
        return None

    def evaluate( self, lp, lm, jets, met = TLorentzVector(),
                  doIntegration = True, coords = None, **kwargs ):
        return self.__call__( lp, lm, jets, met, doIntegration, coords )


#_________________________________________________________________
__matrixElements_map = { 'HWW'        : HWWMatrixElementTool,
                         'WW'         : WWMatrixElementTool,
                         'TT1j'       : tTMatrixElementTool,
                         'TT1jSimple' : tTSimpleMatrixElementTool,
                         'HWW1j'      : HWW1jMatrixElementTool,
                         'WW1j'       : WW1jMatrixElementTool,
                         'DY'         : DYMatrixElementTool,
                         'DY2j'       : DY2jMatrixElementTool,
                         'WFake'      : WFakeMatrixElementTool }

def toolFactory( name = None ):
    global __matrixElements_map
    try:
        return __matrixElements_map[ name ]
    except KeyError:
        pass
    return None

#_________________________________________________________________
def generateRandomVector( ):
    x = TLorentzVector( )
    x.SetPtEtaPhiM( random.uniform( 0, 50 ),
                    random.uniform( -5, 5 ),
                    random.uniform( -math.pi, math.pi ),
                    0. )
    return x

#_________________________________________________________________
def getDefaultIntegrator( rel, minNEval, maxNEval ):
    myint = Divonne()
    myint.setPartitioningRule( 47 )
    myint.setIntegrationRule( 1 )
    myint.setRefinementRule( 1 )
    myint.setMaxPass( 7 )
    myint.setBorder( 0 )
    myint.setChiSqr( 10.0 )
    myint.setMinDev( 0.25 )
    myint.setNComp( 1 )
    myint.setEpsilon( rel, 0 )
    myint.setVerbose( 0 )
    myint.setSampleSet( True )
    myint.setRandomSeed( 1 )
    myint.setPseudoRandom( True )
    myint.setNEval( int(minNEval), int(maxNEval) )
    return myint

#_________________________________________________________________
def getVegasIntegrator( rel, minNEval, maxNEval ):
    myint = root.GSLIntegratorVegas()
    myint.setNComp( 1 )
    myint.setEpsilon( rel, 0 )
    myint.setVerbose( 0 )
    myint.setNSteps( 2 ) 
    myint.setMaxChiSquare( 1.0e6 ) 
    myint.setNEval( int(minNEval), int(maxNEval) )
    return myint

#_________________________________________________________________
def getMiserIntegrator( rel, minNEval, maxNEval ):
    myint = root.GSLIntegratorMiser()
    myint.setNComp( 1 )
    myint.setEpsilon( rel, 0 )
    myint.setVerbose( 0 )
    myint.setNEval( int(minNEval), int(maxNEval) )
    return myint

#_________________________________________________________________
def getLineIntegrator( rel, minNEval = int(1e3), maxNEval = int(1e10) ):
    myint = GSLLine()
    myint.setNComp( 1 )
    myint.setEpsilon( rel, 0 )
    myint.setVerbose( 0 )
    myint.setNEval( int(minNEval), int(maxNEval) )
    return myint

#_________________________________________________________________
def getAlternateIntegrator( rel, minNEval, maxNEval ):
    return getMiserIntegrator( rel, minNEval, maxNEval )

#_________________________________________________________________
def getDummyIntegrator( rel=1, minNEval=0, maxNEval=0 ):
    myint = NullIntegrator( )
    myint.setEpsilon( rel, 0 )
    myint.setNEval( int(minNEval), int(maxNEval) )
    return myint

#_________________________________________________________________
def getHiggsWidth( mass, fname = 'hwidth.txt' ):
    global __wmap
    if len( __wmap ) == 0:
        ifile = open( fname, 'r' )
        for il, l in enumerate( ifile.readlines( ) ):
            m, w = l.split()
            __wmap[int( float( m ) )] = float( w )
        ifile.close()
    return __wmap[int(mass)]

