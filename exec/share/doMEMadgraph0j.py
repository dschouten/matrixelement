#!/usr/bin/python

# Validation of 0j matrix elements

import commands
import os
import sys
import random
import time
import array
import getopt
import glob
import math
import time

def fbranch():
    return array.array('d', [0])

def farrbranch():
    return root.std.vector('double')()

def read( args ):
    myargs = []
    for item in args:
        try:
            myargs.append( float(item) )
            continue
        except:
            pass
        try:
            myargs.append( int(item) )
            continue
        except:
            pass
        myargs.append( str(item) )
    return myargs

usage = 'USAGE: doMEMadgraph.py [-v] [-i] --me=HWW --files=ntuple.root,...,ntuple.root [--treename=physics] ' + \
        '                       --begin=0 --end=100 --output=file.root [--args="[a,b,c]"] '

if __name__ == '__main__':

    DELIM = ','
    
    #
    # todo: specify transfer functions and integration mode
    #
    
    iname   = 'events'
    fmode   = False
    verbose = False
    files   = None
    begin   = None
    end     = None
    outfn   = None
    matrix  = None
    config  = list()
            
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'vim:f:n:b:e:o:a:',
                                   ['me=', 'files=', 'treename=', 'begin=', 'end=', 'output=', 'args='])
    except getopt.GetoptError, error:
        print str( error )
        print usage
        sys.exit( -1 )
    for o, a in opts:
        if o in ( '-v' ):
            verbose = True
        if o in ( '-i' ):
            fmode = True
        if o in ( '--me', '-m' ):
            matrix = a
        if o in ( '--files', '-f' ):
            files = []
            for item in a.split( DELIM ):
                for fn in glob.glob( item ):
                    files.append( fn )
        if o in ( '--treename', '-n' ):
            iname = a
        if o in ( '--begin', '-b' ):
            begin = int( a )
        if o in ( '--end', '-e' ):
            end = int( a )
        if o in ( '--output', '-o' ):
            outfn = a
        if o in ( '--args', '-a' ):
            config = eval( a )

    if None in ( files, matrix ):
        print usage
        sys.exit( -1 )

    from ROOT import gROOT
    
    gROOT.SetBatch( True )

    import ROOT as root
    import msg
    import matrixelement

    log = None
    if verbose: log = msg.msglog( 'doMEMadgraph', 'debug', useColor = False )
    else:       log = msg.msglog( 'doMEMadgraph', 'info',  useColor = False )
    
    # ---------------------------------------------------------------------
    
    ifile = None
    itree = None
    
    #
    # if no output file is specified, and only one input file,
    # then edit the input file "in-place"
    #
        
    otree = None
    
    if outfn != None:
        outdir = os.path.dirname( outfn )
        outfn = os.path.basename( outfn )
        if not os.environ.has_key( 'HEPPYGRID' ):
            try:
                outfn = os.environ['PBS_JOBID'] + '_' + outfn
            except KeyError:
                pass
            try:
                outfn = os.environ['PBS_JOBDESCR'] + '_' + outfn
            except KeyError:
                pass
        if outdir == '':
            outdir = '.'
        outfn = os.path.sep.join( (outdir,outfn) )
        inplace = False
        log.debug( 'Inputs:', files )
        log.debug( 'Output:', outfn )
        itree = root.TChain( iname )
        for fn in files:
            itree.Add( fn )
        ofile = root.TFile.Open( outfn, 'RECREATE' )
        ofile.cd( )
        otree = itree.CloneTree( 0, 'fast SortBasketsByEntry' )
        otree.SetName( matrix )
        otree.SetTitle( matrix )
    elif len( files ) == 1 and outfn == None:
        inplace = True
        ifile = root.TFile.Open( files[0], 'UPDATE' )
        itree = ifile.Get( iname )
        otree = ifile.Get( iname )
        begin, end = 0, itree.GetEntries()        
    else:
        log.error( "must specifiy output file" )
        sys.exit( -1 )

    me_factors = farrbranch()
    me         = fbranch()
    me_error   = fbranch()
    me_stat    = fbranch()
    me_neval   = fbranch()
    me_prob    = fbranch()

    kin_metrel = fbranch()
    kin_met    = fbranch()
    kin_dphi   = fbranch()

    br_suff = '_Weight_%s'%( matrix )
    
    if len( config ) != 0:
        br_suff += '_'.join( [ str(cfg) for cfg in config ] )

    me_factors_br = otree.Branch( 'ME_factors' + br_suff, me_factors )
    me_br         = otree.Branch( 'ME' + br_suff, me, 'ME' + br_suff + '/D' )
    me_error_br   = otree.Branch( 'ME_error' + br_suff, me_error, 'ME_error' + br_suff + '/D' )
    me_stat_br    = otree.Branch( 'ME_stat' + br_suff, me_stat, 'ME_stat' + br_suff + '/D' )
    me_neval_br   = otree.Branch( 'ME_neval' + br_suff, me_neval, 'ME_neval' + br_suff + '/D' )
    me_prob_br    = otree.Branch( 'ME_prob' + br_suff, me_prob, 'ME_prob' + br_suff + '/D' )

    kin_metrel_br = otree.Branch( 'KIN_metrel' + br_suff, kin_metrel, 'KIN_metrel' + br_suff + '/D' )
    kin_met_br    = otree.Branch( 'KIN_met' + br_suff, kin_met, 'KIN_met' + br_suff + '/D' )
    kin_dphi_br   = otree.Branch( 'KIN_dphi' + br_suff, kin_dphi, 'KIN_dphi' + br_suff + '/D' )

    # ---------------------------------------------------------------------

    #
    # configure the integrator
    #

    myint = matrixelement.Divonne()
    myint.setPartitioningRule( 47 )
    myint.setIntegrationRule( 1 )
    myint.setRefinementRule( 1 )
    myint.setMaxPass( 7 )
    myint.setBorder( 0 )
    myint.setChiSqr( 10.0 )
    myint.setMinDev( 0.25 )  
    myint.setNComp( 1 )
    myint.setEpsilon( 0.025, 0 )
    myint.setVerbose( 0 )
    myint.setSampleSet( True )
    myint.setRandomSeed( 1 )
    myint.setPseudoRandom( True )
    myint.setNEval( int(1.0e4), int(1.0e8) )
    
    # ---------------------------------------------------------------------

    #
    # configure the matrix element (integrand)
    #

    if not os.path.exists( 'cteq6l.tbl' ):
        log.debug( 'Fetching CTEQ6 PDF' )
        commands.getstatusoutput( 'ln -s $HEPPY/higgs/mymeanalysis/exec/run/cteq6l.tbl .' )

    # ! FIXME ! the integration strategy is hard-coded here
    config += [ root.WWLeptonsBaseIntegrand.kNEUTRINO_4D, None ]
    
    myme = matrixelement.toolFactory( matrix )( myint, args = config )

    integrand = myme.matrixElement

    if 'hww' in matrix.lower():
        integrand.setWidth( matrixelement.getHiggsWidth( integrand.mass() ) )
        integrand.setUseSM( False )
    
    integrand.setIntegrationLimits( integrand.getIPTA(), 0, 500 )
    integrand.setIntegrationLimits( integrand.getIPHIA(), -math.pi, math.pi )
    integrand.setIntegrationLimits( integrand.getIPZA(), -500, 500 )
    integrand.setIntegrationLimits( integrand.getIPZB(), -500, 500 )
    
    integrand.initialize()

    # ---------------------------------------------------------------------

    #
    # main event loop
    #
        
    log.info( 'Begin processing events', (begin,end) )
    
    if begin == None: begin = 0
    if end == None: end = itree.GetEntries()
    
    if begin > itree.GetEntries():
        sys.exit( 0 )
    
    for ievt in xrange( begin, min( end, itree.GetEntries() ) ):
        nb = itree.GetEntry( ievt )
        if nb <= 0:
            break

        log.info( 'Event #', ievt - begin )

        lp = root.TLorentzVector( itree.px_ep, itree.py_ep, itree.pz_ep, itree.E_ep )
        lm = root.TLorentzVector( itree.px_em, itree.py_em, itree.pz_em, itree.E_em )

        nul = root.TLorentzVector( itree.px_ve, itree.py_ve, itree.pz_ve, itree.E_ve )
        nur = root.TLorentzVector( itree.px_vem, itree.py_vem, itree.pz_vem, itree.E_vem )
        
        met = root.TLorentzVector( (nul+nur).Px(), (nul+nur).Py(), 0, (nul+nur).Pt() )

        minDeltaPhi = math.pi / 2

        for obj in [lp, lm]:
            if abs( obj.DeltaPhi( met ) ) < minDeltaPhi:
                minDeltaPhi = abs( obj.DeltaPhi( met ) )

        kin_metrel[0] = met.Pt() * math.sin( minDeltaPhi )
        kin_met[0]    = met.Pt()
        kin_dphi[0]   = abs( lp.DeltaPhi( lm ) )

        log.debug( '\tlepton pT:', lp.Pt(), lm.Pt() )
        log.debug( '\tmissing ET:', met.Pt() )
        log.debug( '\tmetREL:', met.Pt() * math.sin( minDeltaPhi ) )

        result = { 'me' : 0, 'error' : 1, 'status' : 1, 'neval' : 0, 'prob' : 0 }

        result = myme.evaluate( lp, lm, met, nul, nur,
                                doIntegration = (not fmode),
                                coords = [ nur.Pt(), nur.Phi(), nur.Pz(), nul.Pz() ] )
        
        me[0]        = result['me']
        me_error[0]  = result['error']
        me_stat[0]   = result['status']
        me_neval[0]  = result['neval']
        me_prob[0]   = result['prob']

        me_factors.clear()
        if myme.results != None and len( myme.results ) > 0:
            for xfact in myme.results: me_factors.push_back( xfact )
                
        log.debug( '\tme =', '%0.5g'%( result['me'] ) )
        log.debug( '\terror =', '%0.5g'%( result['error'] ) )
        log.debug( '\tstat =', '%0.5g'%( result['status'] ) )
        log.debug( '\tprob =', '%0.5g'%( result['prob'] ) )
        log.debug( '\tneval =', result['neval'] )

        if not inplace:
            log.debug( 'filling output tree' )
            otree.Fill( )
        else:
            me_br.Fill( )
            me_error_br.Fill()
            me_stat_br.Fill()
            me_neval_br.Fill()
            me_prob_br.Fill()
            me_factors_br.Fill()
    
    if not inplace:
        log.debug( 'writing output tree to file' )
        otree.Write()
        ofile.Close()
    else:
        otree.Write( '', root.TObject.kOverwrite)

    log.info( 'Done processing events' )
