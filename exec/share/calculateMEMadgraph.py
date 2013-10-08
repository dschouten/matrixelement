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

cuts = None

try:
    opts, args = getopt.getopt( sys.argv[1:], 'vat:c:x:i:o:b:e:r:s:',
                                ['timeout=','cfg=','exec=','input=','output=','begin=','end=','random=','seed='] )
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
class Tauon( ROOT.TLorentzVector ):
    def __init__( self, tlv, charge ):
        ROOT.TLorentzVector.__init__( self, tlv.X(), tlv.Y(), tlv.Z(), tlv.E() )
        self.charge = charge
        self.type = 'tauon'

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

def buffer_result( r ):
    rfile = open( 'result_buffer.pickle', 'w' )
    pickle.dump( r, rfile )
    rfile.close()

def id_generator( size=6, chars=string.ascii_uppercase + string.digits ):
    return ''.join( random.choice(chars) for x in range(size) )

#_________________________________________________________________

pdg = ROOT.TDatabasePDG( )

ROOT.gSystem.Load("libExRootAnalysis.so")

electronPDG = 11
muonPDG     = 13
tauPDG      = 15

lepMasses = { 13 : pdg.GetParticle( 'mu-' ).Mass( ),
              11 : pdg.GetParticle( 'e-' ).Mass( ),
              15 : pdg.GetParticle( 'tau-' ).Mass( ) }

lepTypes  = { 13 : Muon, 11 : Electron, 15 : Tauon }

GeV = 1000.

allowNull = True

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
ROOT.gErrorIgnoreLevel = ROOT.kError

log.info( 'will perform [%s] ME calculation in [%d] jet bin'%( title, myMEJetN ) )
    
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
else:
    ifilen = ifilen.pop()

ifile = ROOT.TFile.Open( ifilen, 'read' )
if ifile == None or ifile.IsZombie():
    log.fatal( 'file [%s] not found'%( ifilen ) )
    
itree = ifile.Get( 'LHEF' )
itree.SetBranchStatus( '*', True )

# ------------------------------------------------------------------------
#
# load the output file
#
# ------------------------------------------------------------------------

ofile = ROOT.TFile.Open( ofilen, 'recreate' )
if ofile == None or ofile.IsZombie():
    log.fatal( 'file [%s] can not be accessed'%( ofilen ) )

otree = ROOT.TTree( title, title )

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

spectatortrees = []
if append:
    for key in ifile.GetListOfKeys():
        if key.GetName() == title:
            continue
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

nselected = len(entries)

itree.GetEntry(0)

elist = None
if cuts != None:
    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    if ROOT.TTreeFormula( "cuts", cuts, itree ).GetNdim() != 0: # << valid TCuts specified
        itree.Draw( ">>entrylist", cuts )
        elist = ROOT.gROOT.FindObject( "entrylist" )
        if elist != None:
            nselected = elist.GetN()
    ROOT.gErrorIgnoreLevel = ROOT.kWarning

log.info( 'will run over [%d] selected events of [%d] total'%( nselected, len(entries) ) )

for ientry in entries:

    if elist != None and not elist.Contains( ientry ):
        continue

    measured_lp       = tlv() 
    measured_lm       = tlv() 
    inferred_met      = tlv()

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

    ievent = itree.Event[0].Number 
    irun   = itree.Event[0].ProcessID
    
    result = { 'me' : -1, 'error' : -1, 'neval' : -1, 'status' : -1, 'prob' : -1,
               'ievent' : ievent, 'irun' : irun, 'ientry' : ientry, 'time' : -1, 'itry' : -1 }

    # ------------------------------------------------------------------------
    # buffer the results table (if integration times out, at least sensible
    # entries will be put in the tree if using wrapperME script)
    # ------------------------------------------------------------------------
    
    buffer_result( result )
    
    # ------------------------------------------------------------------------
    # get lepton flavours and charges
    # ------------------------------------------------------------------------

    for ipart, part in enumerate( itree.Particle ):
        if abs( part.PID ) == electronPDG or abs( part.PID ) == muonPDG:
            if part.PID < 0:
                measured_lp.SetPxPyPzE( part.Px, part.Py, part.Pz, part.E )
                measured_lp = lepTypes[abs(part.PID)]( measured_lp,  1) 
            else:
                measured_lm.SetPxPyPzE( part.Px, part.Py, part.Pz, part.E )
                measured_lm = lepTypes[abs(part.PID)]( measured_lm, -1)
        if abs( part.PID ) in [ 1, 2, 3, 4, 5, 21 ] and part.Status == 1:
            jets += [ tlv( part.Px, part.Py, part.Pz, part.E ) ]
    
    log.info( '\t lp : %0.2f, %0.2f, %0.2f'%( measured_lp.Pt(), measured_lp.Eta(), measured_lp.Phi() ) )
    log.info( '\t lm : %0.2f, %0.2f, %0.2f'%( measured_lm.Pt(), measured_lm.Eta(), measured_lm.Phi() ) )

    # ------------------------------------------------------------------------
    # get the total four-vector for the outgoing branches
    # ------------------------------------------------------------------------
        
    inferred_met = -(measured_lp + measured_lm)
    for j in jets:
        inferred_met -= j
    inferred_met.SetZ( 0 )
    inferred_met.SetE( inferred_met.Pt() )
    
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
        result.update( myMECalculation( measured_lp, measured_lm, jets, inferred_met ) )            
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
