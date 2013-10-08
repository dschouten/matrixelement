#!/usr/bin/python

#---------------------------------------------------------------------
#
# This script reads MVA ntuples and calculates ME using calculateME
# and imposes a timeout on the calculation by polling
#
#---------------------------------------------------------------------

from ROOT import gROOT
gROOT.SetBatch( True )

import ROOT

import time
import os
import signal
import sys
import getopt
import glob
import pickle
import array

from commands import getstatusoutput as shell
from subprocess import Popen as popen
from glob import glob
from array import array

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

# ------------------------------------------------------------------------

polling     = 0.25 # << timestep for polling
timeout     = 360 # << total timeout, in seconds
verbose     = False # << verbose output
append      = False # << append non-conflicting objects from input file
ifilen      = None # << input filename
ofilen      = None # << output filename
cfgmodules  = [] # << configuration modules
begevent    =  0 # << first event to process
endevent    = -1 # << stop upon reaching this event

try:
    opts, args = getopt.getopt(sys.argv[1:],
                               'vat:c:i:o:b:e:',
                               ['timeout=','cfg=','input=','output=','begin=','end='])
except getopt.GetoptError, error:
    print str( error )
    print 'USAGE: [-v] [-a] [-t,--timeout=300] --cfg=cfg.py --input=ww.root --output=me_ww.root'
    print '       [-b,--begin=] [-e,--end]'
    sys.exit( -1 )

for o, a in opts:
    if o in ( '-v' ):
        verbose = True
    if o in ( '-a' ):
        append = True
    if o in ( '-t', '--timeout' ):
        timeout = abs(int( a )) + 5 
    if o in ( '-c', '--cfg' ):
        cfgmodules += [ a ]
    if o in ( '-i', '--input' ):
        ifilen = sorted( a.split( ',' ) )
    if o in ( '-o', '--output' ):
        ofilen = a
    if o in ( '-b', '--begin' ):
        begevent = int( a )
    if o in ( '-e', '--end' ):
        endevent = int( a )

if len(cfgmodules)==0 or ofilen==None or ifilen==None:
    print 'USAGE: [-v] [-t,--timeout=300] --cfg=cfg.py --input=ww.root --output=me_ww.root'
    sys.exit( -1 )

# ------------------------------------------------------------------------
#
# determine configuration
#
# ------------------------------------------------------------------------

import matrixelement
import msg

ROOT.GlobalFlags.debug = verbose
ROOT.gErrorIgnoreLevel = ROOT.kFatal

pdg = ROOT.TDatabasePDG( )

GeV = 1000.

for module in cfgmodules:
    if os.path.exists( module ):
        execfile( module )
    else:
        print 'ERROR: configuration [%s] does not exist'%( module )
        sys.exit( -1 )

log = None
if verbose: log = msg.msglog( 'wrapper', 'debug', useColor = False )
else:       log = msg.msglog( 'wrapper', 'info',  useColor = False )

title = normalize( myMECalculation.mehndl.name( ) )

## log.info( 'will perform [%s] ME calculation in [%d] jet bin'%( title, myMEJetN ) )

# ------------------------------------------------------------------------
#
# inspect the input file(s)
#
# ------------------------------------------------------------------------

log.info( 'peak %s'%( ifilen ) )

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
    
totalNumEvents = ifile.Get( 'HWWTree' ).GetEntries()

log.info( 'found [%d] events total'%( totalNumEvents ) )

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
# build the calculation command
#
# ------------------------------------------------------------------------

defcmd  = [ 'python', 'calculateME.py',
            '-t %d'%( int(timeout + polling * 50) ),
            '--input=%s'%( ifilen ) ]

for cfg in cfgmodules:
    defcmd += ['--cfg=%s'%( cfg )]

# ------------------------------------------------------------------------
#
# run the calculation for each event 
#
# ------------------------------------------------------------------------

a = min( begevent, totalNumEvents )
b = min( endevent, totalNumEvents )
if b <= 0:
    b = totalNumEvents

err = open( 'my.err', 'w' )

for ievt in xrange( a, b ):
    cmd = defcmd[:]
    
    if verbose and ievt == 0:
        cmd += ['-v']
        
    cmd += ['--output=result_buffer.root']
    cmd += ['--begin=%d'%( ievt )]
    cmd += ['--end=%d'%( ievt + 1 )]

    result = { 'me' : -1, 'error' : -1, 'neval' : -1, 'status' : -1, 'prob' : -1,
               'ievent' : -1, 'irun' : -1, 'ientry' : ievt, 'time' : -1, 'itry' : -1 }
    
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
        
    log.info( ">>>>>>>>>>>>> start of event <<<<<<<<<<<<<" )
    log.info( 'forking "%s" ...'%( ' '.join( cmd ) ) )

    proc = popen( cmd, stderr = err )

    pstatus = proc.poll()

    ttime = polling
    while ttime < timeout and pstatus == None:
        pstatus = proc.poll()
        if pstatus == None:
            ttime += polling
            time.sleep( polling )

    if pstatus == None:
        log.warning( 'calculation timed out for event #%d [%d seconds]'%( ievt, timeout ) )
        os.kill( proc.pid, signal.SIGKILL )
        ## stat, out = commands.getstatusoutput( 'kill -KILL %d'%( proc.pid ) )
    else:
        log.info( 'status = %d'%( pstatus ) )
    
    try:
        rfile  = open( 'result_buffer.pickle', 'r' )
        result = pickle.load( rfile )
        rfile.close()
    except:
        log.warning( 'could not load output file, ME calculation set to -1' )

    try:
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
    except:
        log.warning( 'could not set results from ME calculation' )
        
    otree.Fill() 

    if append:
        for pair in spectatortrees:
            orig, clone = pair
            orig.GetEntry( ievt )
            clone.Fill( )

    log.info( ">>>>>>>>>>>>> end of event <<<<<<<<<<<<<" )
    log.info( "" )

ofile.cd( )

if append:
    for pair in spectatortrees:
        orig, clone = pair
        clone.Write( )
        
otree.Write( )
nout = otree.GetEntries()

log.info( 'total read in / processed: ', b - a, '/', nout )

ifile.Close()
ofile.Close()

if (b-a) != nout:
    log.error( 'number of output events does not match input' )
    shell( 'rm -rf %s'%( ofilen ) ) # << remove output to force job failure
   
if hadded:
    shell( 'rm -rf %s'%( ifilen ) )
    
log.info( 'finished [%s/%s]'%( ifilen, title ) )

err.close()

stat, out = shell( 'cat my.err | grep -v TEnvRec::ChangeValue' ) # << ignore ROOT warnings
print out
