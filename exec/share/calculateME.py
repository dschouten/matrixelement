#!/usr/bin/python

#---------------------------------------------------------------------
#
# This script reads MVA ntuples and calculates ME
#
# Is somewhat agnostic about specific ME configuration but still
# hard-coded for 0 & 1 jet lvlv final states
#
# Input:
#  cfg=   : python module(s) that configure the calculation
#  input= : ntuple file with the events to loop over
#  ouput= : the file to store the calculations in
#  
#---------------------------------------------------------------------

import sys
import os
import random
import math
import array
import getopt
import glob
import signal
import time
import pickle
import string

from commands import getstatusoutput as shell
from glob import glob
from array import array

#_________________________________________________________________

append      = False
verbose     = False
ifilen      = None
ofilen      = None
begevent    =  0
endevent    = -1
randomsel   = -1
randomseed  =  1
cfgmodules  = list() 
cfgexec     = ''
seltree     = 'HWWTree'

cuts = None

try:
    opts, args = getopt.getopt( sys.argv[1:], 'vat:c:x:i:o:b:e:r:s:q:',
                                ['timeout=','cfg=','exec=','input=','output=','begin=','end=','random=','seed=','tree='] )
except getopt.GetoptError, error:
    print str( error )
    print 'USAGE: [-v] [-a] --cfg=cfg.py [--exec="mehndl.setX()"] --input=ww.root --output=me_ww.root'
    print '       [--begin=0] [--end=1] [--random=-1] [--seed=1]'
    sys.exit( -1 )

for o, a in opts:
    if o in ( '-v' ):
        verbose = True
    if o in ( '-a' ):
        append = True
    if o in ( '-o', '--output' ):
        ofilen = a
    if o in ( '-i', '--input' ):
        ifilen = sorted( a.split( ',' ) )
    if o in ( '-c', '--cfg' ):
        cfgmodules += [ a ]
    if o in ( '-x', '--exec' ):
        cfgexec = a
    if o in ( '-b', '--begin' ):
        begevent = int( a )
    if o in ( '-e', '--end' ):
        endevent = int( a )
    if o in ( '-r', '--random' ):
        randomsel = int( a )
    if o in ( '-s', '--seed' ):
        randomseed = int( a )
    if o in ( '-q', '--tree' ):
        seltree = a

if len(cfgmodules)==0 or ofilen==None or ifilen==None:
    print 'USAGE: [-v] [-a] --cfg=cfg.py --input=ww.root --output=me_ww.root'
    print '       [--begin=0] [--end=1] [--random=10] [--seed=1]'
    sys.exit( -1 )

#_________________________________________________________________

from ROOT import gROOT
gROOT.SetBatch( True )

import ROOT
import matrixelement
import msg

#_________________________________________________________________
class Electron( ROOT.TLorentzVector ):
    def __init__( self, tlv, charge ):
        ROOT.TLorentzVector.__init__( self, tlv.X(), tlv.Y(), tlv.Z(), tlv.E() )
        self.charge = charge
        self.type = 'electron'
        
#_________________________________________________________________
class Muon( ROOT.TLorentzVector ):
    def __init__( self, tlv, charge ):
        ROOT.TLorentzVector.__init__( self, tlv.X(), tlv.Y(), tlv.Z(), tlv.E() )
        self.charge = charge
        self.type = 'muon'

#_________________________________________________________________
def fbranch():
    return array('d', [0])

def ibranch():
    return array('i', [0])

def farrbranch():
    return ROOT.std.vector('double')()

def normalize( s ):
    n = s
    for char in ',./<>?;\':"[]\{}\!@#$%^&*()+-=':
        n = n.replace( char, '_' )
    return n

def buffer_result( r, t=None ):
    rfile = open( 'result_buffer.pickle', 'w' )
    pickle.dump( r, rfile )
    rfile.close()
    if t!=None:
        t.Fill()

def id_generator( size=6, chars=string.ascii_uppercase + string.digits ):
    return ''.join( random.choice(chars) for x in range(size) )

#_________________________________________________________________

pdg = ROOT.TDatabasePDG( )

lepMasses = { 13 : pdg.GetParticle( 'mu-' ).Mass( ), 11 : pdg.GetParticle( 'e-' ).Mass( ) }
lepTypes  = { 13 : Muon, 11 : Electron }

GeV = 1000.

allowNull = True

applyBoostCorrection = False
boostCalculator = None

for module in cfgmodules:
    if os.path.exists( module ):
        execfile( module )
    else:
        print 'ERROR: configuration [%s] does not exist'%( module )
        sys.exit( -1 )

if len( cfgexec ) > 0:
    print
    print
    print cfgexec
    print
    print
    exec( cfgexec )

title = normalize( myMECalculation.mehndl.name( ) )

log = None
if verbose:
    log = msg.msglog( 'main', 'debug', useColor = False )
else:
    log = msg.msglog( 'main', 'info',  useColor = False )

ROOT.GlobalFlags.debug = verbose
ROOT.GlobalFlags.max_debug = 10
ROOT.gErrorIgnoreLevel = ROOT.kError

stat, host = shell( 'hostname' )
stat, time = shell( 'date' )
log.info( 'running on [%s] @ %s'%( host, time ) )
log.info( 'will perform [%s] ME calculation on file(s) %s'%( title, str(ifilen)  ) )
    
# ------------------------------------------------------------------------
#
# merge and load the input file(s)
#
# ------------------------------------------------------------------------

hadded = False
if len( ifilen ) > 1:
    ifilen = '%s.root'%( id_generator() )
    stat, out = shell( ' '.join( ['hadd -f %s.root'%( )] + ifilen ) )
    log.debug( out )
    hadded = True
    log.info( 'merged files into single input [%s]'%( ifilen ) )
else:
    ifilen = ifilen.pop()

ifile = ROOT.TFile.Open( ifilen, 'read' )
if ifile == None or ifile.IsZombie():
    log.fatal( 'file [%s] not found'%( ifilen ) )
    
itree = ifile.Get( seltree )
itree.SetBranchStatus( '*', True )

for key in ifile.GetListOfKeys():
    if key.ReadObj().Class() == 'TTree':
        itree.AddFriend( key.GetName() )

# ------------------------------------------------------------------------
#
# load the output file
#
# ------------------------------------------------------------------------

ofile = ROOT.TFile.Open( ofilen, 'recreate' )
if ofile == None or ofile.IsZombie():
    log.fatal( 'file [%s] can not be accessed'%( ofilen ) )

otreename = title
otreedesc = title

spectatortrees = []

if append:
    keynames = []
    for key in ifile.GetListOfKeys():
        keynames += [ key.GetName() ]
        obj = key.ReadObj()
        if obj.Class().InheritsFrom( ROOT.TTree.Class() ):
            t = ifile.Get( key.GetName() )
            t.SetBranchStatus( '*', True )
            ofile.cd( )
            c = t.CloneTree( 0 )
            spectatortrees.append( (t,c) )
        if obj.Class().InheritsFrom( ROOT.TH1.Class() ):
            ofile.cd( )
            buf = obj.Clone( obj.GetName() )
            buf.Write( )
    buf = title
    inc = 0
    while buf in keynames: 
        buf = title + '_%02d'%( inc )
        inc += 1
    otreename = buf

otree = ROOT.TTree( otreename, otreedesc )

me_value  = fbranch()
me_error  = fbranch()
me_stat   = fbranch()
me_neval  = fbranch()
me_prob   = fbranch()
me_ntries = ibranch()
me_time   = fbranch()
me_entry  = ibranch()
me_evtnum = ibranch()
me_runnum = ibranch()
me_norm   = fbranch()

me_value_br  = otree.Branch( 'me'     , me_value  , 'me/D' )
me_error_br  = otree.Branch( 'error'  , me_error  , 'error/D' )
me_stat_br   = otree.Branch( 'stat'   , me_stat   , 'stat/D' )
me_neval_br  = otree.Branch( 'neval'  , me_neval  , 'neval/D' )
me_prob_br   = otree.Branch( 'prob'   , me_prob   , 'prob/D' )
me_ntries_br = otree.Branch( 'ntries' , me_ntries , 'ntries/I' )
me_time_br   = otree.Branch( 'time'   , me_time   , 'time/D' )
me_entry_br  = otree.Branch( 'entry'  , me_entry  , 'entry/I' )
me_evtnum_br = otree.Branch( 'evtnum' , me_evtnum , 'evtnum/I' )
me_runnum_br = otree.Branch( 'runnum' , me_runnum , 'runnum/I' )
me_norm_br   = otree.Branch( 'norm'   , me_norm   , 'norm/D' )

# ------------------------------------------------------------------------
#
# loop over the events
#
# ------------------------------------------------------------------------

timer = ROOT.TStopwatch( )

a = min( begevent, itree.GetEntries() )
b = min( endevent, itree.GetEntries() )
if b <= 0:
    b = itree.GetEntries()

entries = range( a, b )

random.seed( randomseed ) 

if randomsel > 0:
    allowed = entries[:]
    entries = []
    while len(entries) < randomsel and len(allowed) > 0:
        index = int( random.uniform(0, len(allowed)) )
        entries += [ allowed.pop( index ) ]

tlv = ROOT.TLorentzVector

nprocessed = 0

elist = None
if cuts != None and len(cuts) > 0:
    itree.Draw( ">>entrylist", cuts )
    elist = ROOT.gROOT.FindObject( "entrylist" )
    if elist != None:
        log.info( 'using [%d] selected events of [%d] total'%( elist.GetN(), itree.GetEntries() ) )
    else:
        log.warning( 'failed to apply cuts [%s]'%( cuts ) )

for ientry in entries:

    measured_lp       = tlv() 
    measured_lm       = tlv() 
    measured_met      = tlv() 
    measured_jets     = []

    # ------------------------------------------------------------------------
    # initialize the results table 
    # ------------------------------------------------------------------------
    
    result = { 'me' : -1, 'error' : -1, 'neval' : -1, 'status' : -1, 'prob' : -1,
               'ievent' : -1, 'irun' : -1, 'ientry' : ientry, 'time' : -1, 'itry' : -1 }
    buffer_result( result )
    
    nb = itree.GetEntry( ientry )
    if nb <= 0:
        log.error( 'reading event #%d'%( ientry ) )
        sys.exit( -1 )
    
    log.debug( 'reading event #%d'%( ientry ) )

    ievent = getattr( itree, 'EventNumber' ) 
    irun   = getattr( itree, 'RunNumber' )

    log.info( "event number: ", ievent )
    
    if hasattr( itree, 'mc_channel_number' ):
        irun = getattr( itree, 'mc_channel_number' )
        
    result = { 'me' : -1, 'error' : -1, 'neval' : -1, 'status' : -1, 'prob' : -1,
               'ievent' : ievent, 'irun' : irun, 'ientry' : ientry, 'time' : -1, 'itry' : -1 }

    me_value[0]  = result['me']
    me_error[0]  = result['error']
    me_stat[0]   = result['status']
    me_neval[0]  = result['neval']
    me_prob[0]   = result['prob']
    me_ntries[0] = result['itry']
    me_entry[0]  = result['ientry']
    me_time[0]   = result['time']
    me_evtnum[0] = result['ievent']
    me_runnum[0] = result['irun']
    me_norm[0]   = 1 ## possible event-by-event ME normalization
    
    # ------------------------------------------------------------------------
    # buffer the results table (if integration times out, at least sensible
    # entries will be put in the tree if using wrapperME script)
    # ------------------------------------------------------------------------
    
    buffer_result( result )

    # ------------------------------------------------------------------------
    # cuts applied here; use TEntryList with proper TCuts definition
    # ------------------------------------------------------------------------
       
    if ( elist != None and (not elist.Contains( ientry )) ):
        ## otree.Fill() ##
        continue
    
    # ------------------------------------------------------------------------
    # get lepton flavours and charges
    # ------------------------------------------------------------------------

    lepID0 = itree.lepID0 
    lepID1 = itree.lepID1

    random.seed( ) ## uses system time ... randomizes w.r.t. seed used for entry selection above

    if lepID0 * lepID1 > 0: ## deal with same-sign charged leptons by randomizing the charge assignment
        s = -2 * (random.uniform(0,1) >= 0.5) + 1
        lepID0 = ( s ) * abs(lepID0)
        lepID1 = (-s ) * abs(lepID1)
    
    if lepID0 > 0: ## lepID contains charge & PDG code information
        measured_lp.SetPtEtaPhiM( itree.lepPt0 / GeV, itree.lepEta0, itree.lepPhi0, lepMasses[abs(lepID0)] )
        measured_lm.SetPtEtaPhiM( itree.lepPt1 / GeV, itree.lepEta1, itree.lepPhi1, lepMasses[abs(lepID1)] )
        measured_lp = lepTypes[abs(lepID0)]( measured_lp,  1)
        measured_lm = lepTypes[abs(lepID1)]( measured_lm, -1)
    else: 
        measured_lm.SetPtEtaPhiM( itree.lepPt0 / GeV, itree.lepEta0, itree.lepPhi0, lepMasses[abs(lepID0)] )
        measured_lp.SetPtEtaPhiM( itree.lepPt1 / GeV, itree.lepEta1, itree.lepPhi1, lepMasses[abs(lepID1)] )
        measured_lm = lepTypes[abs(lepID0)]( measured_lm, -1)
        measured_lp = lepTypes[abs(lepID1)]( measured_lp,  1)

    log.info( '\t lp : %0.2f, %0.2f, %0.2f'%( measured_lp.Pt(), measured_lp.Eta(), measured_lp.Phi() ) )
    log.info( '\t lm : %0.2f, %0.2f, %0.2f'%( measured_lm.Pt(), measured_lm.Eta(), measured_lm.Phi() ) )
    log.info( '\t ll : %0.2f'%( itree.Ptll / 1000 ) )

    # ------------------------------------------------------------------------
    # get the jets for the ME
    # ------------------------------------------------------------------------
       
    for ijet in xrange( itree.m_jet_n ):
        measured_jets.append( tlv() )
        thejet = measured_jets[-1]
        try:
            thejet.SetPtEtaPhiM( itree.m_jet_pt[ijet] / GeV, itree.m_jet_eta[ijet], itree.m_jet_phi[ijet], 0. )
        except AttributeError:
            if ijet < 2:
                thejet.SetPtEtaPhiM( getattr( itree, 'jetPt%d'%( ijet ) ) / GeV,
                                     getattr( itree, 'jetEta%d'%( ijet ) ),
                                     getattr( itree, 'jetPhi%d'%( ijet ) ), 0 )
            else:
                continue
        log.info( '\t jet #%d: %0.2f, %0.2f, %0.2f'%( ijet+1, thejet.Pt(), thejet.Eta(), thejet.Phi() ) )
    
    # ------------------------------------------------------------------------
    # get the MET & recoil estimate
    # ------------------------------------------------------------------------
    
    measured_met.SetPxPyPzE( itree.MET_x / GeV, itree.MET_y / GeV, 0., itree.MET / GeV )
    log.info( '\t MET: %0.2f, %0.2f'%( measured_met.Pt(), measured_met.Phi() ) )

    if boostCalculator != None:
        execfile( 'configs/%s'%( boostCalculator ) ) # << put the boost snippet in another file        
       
    itry = 0

    result['itry'] = itry ## indicate that integration attempted
    buffer_result( result )
    
    timer.Start( )

    # ------------------------------------------------------------------------
    # start integrating the ME calculation
    # ------------------------------------------------------------------------
    
    integrators = [ 'default', 'backup', 'failsafe' ]

    while itry < len( integrators ) and (result['status'] != 0 or (result['me']==0 and not allowNull)):
        ikey = integrators[itry]
        myMECalculation.setIntegrator( myMEIntegrators[ikey] )
        myMECalculation.mehndl.setIteration( itry )
        log.debug( '\t integration using [%s] integrator (%s)'%( myMEIntegrators[ikey].getName(), ikey ) )
        pvalue = myMECalculation( measured_lp, measured_lm, measured_jets, measured_met )
        result.update( pvalue )
        itry += 1
        result['itry'] = itry
        buffer_result( result )
        log.info( 'attempt #%d:'%( itry ), result )
    
    timer.Stop( )
    result['time'] = timer.CpuTime()
    
    log.info( 'final result:', result )
    
    me_value[0]  = result['me']
    me_error[0]  = result['error']
    me_stat[0]   = result['status']
    me_neval[0]  = result['neval']
    me_prob[0]   = result['prob']
    me_ntries[0] = result['itry']
    me_entry[0]  = result['ientry']
    me_time[0]   = result['time']
    me_evtnum[0] = result['ievent']
    me_runnum[0] = result['irun']
    me_norm[0]   = 1 ## possible event-by-event ME normalization

    log.info( 'completed event #%d in %0.2f (s)'%( nprocessed + 1, me_time[0] ) )

    sys.stdout.flush()

    if append:
        for pair in spectatortrees:
            orig, clone = pair
            orig.GetEntry( ientry )
            if orig.Class().GetName() == 'TNtuple':
                clone.Fill( array( 'f', ( getattr( orig, br.GetName() ) for br in orig.GetListOfBranches() ) ) )
            else:
                clone.Fill( )
    
    otree.Fill( )
    
    buffer_result( result )

    nprocessed += 1
    
if append:
    for pair in spectatortrees:
        orig, clone = pair
        clone.Write( )

otree.Write( )
ofile.Close( )

if hadded:
    shell( 'rm -rf %s'%( ifilen ) )
