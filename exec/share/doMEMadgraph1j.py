#!/usr/bin/python

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
    
    itree   = 'events'
    verbose = False
    files   = None
    begin   = None
    end     = None
    outfn   = None
    matrix  = None
    config  = list()
            
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'vtm:f:n:b:e:o:a:',
                                   [ 'me=', 'files=', 'treename=',
                                     'begin=', 'end=', 'output=', 'args=' ])
    except getopt.GetoptError, error:
        print str( error )
        print usage
        sys.exit( -1 )
    for o, a in opts:
        if o in ( '-v' ):
            verbose = True
        if o in ( '--me', '-m' ):
            matrix = a
        if o in ( '--files', '-f' ):
            files = []
            for item in a.split( DELIM ):
                for fn in glob.glob( item ):
                    files.append( fn )
        if o in ( '--treename', '-n' ):
            itree = a
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

    import ROOT as root
    root.gROOT.SetBatch( True )
    
    from tools import rmap

    import defaults, units, msg, style
    import matrixelement

    log = None
    if verbose: log = msg.msglog( 'doMEMadgraph1J', 'debug' )
    else:       log = msg.msglog( 'doMEMadgraph1J', 'info' )
    
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
        iplace = False
        log.debug( 'Inputs:', files )
        log.debug( 'Output:', outfn )
        itree = root.TChain( itree )
        for fn in files:
            itree.Add( fn )
        ofile = root.TFile.Open( outfn, 'RECREATE' )
        ofile.cd( )
        otree = itree.CloneTree( 0, 'fast SortBasketsByEntry' )
        otree.SetName( matrix )
        otree.SetTitle( matrix )
    elif len( files ) == 1 and outfn == None:
        iplace = True
        ifile = root.TFile.Open( files[0], 'UPDATE' )
        itree = ifile.Get( itree )
        otree = ifile.Get( itree )
        begin, end = 0, itree.GetEntries()        
    else:
        log.error( "must specifiy input/output file" )
        sys.exit( -1 )
        
    me_factors = farrbranch()
    me        = fbranch()
    me_error  = fbranch()
    me_stat   = fbranch()
    me_neval  = fbranch()
    me_prob   = fbranch()
    jet_prob  = fbranch()

    br_suff = '_Weight_%s'%( matrix )
    if len( config ) != 0: br_suff += '_'.join( [ str(cfg) for cfg in config ] )

    me_factors_br = otree.Branch( 'ME_factors' + br_suff, me_factors )
    me_br         = otree.Branch( 'ME' + br_suff, me, 'ME' + br_suff + '/D' )
    me_error_br   = otree.Branch( 'ME_error' + br_suff, me_error, 'ME_error' + br_suff + '/D' )
    me_stat_br    = otree.Branch( 'ME_stat' + br_suff, me_stat, 'ME_stat' + br_suff + '/D' )
    me_neval_br   = otree.Branch( 'ME_neval' + br_suff, me_neval, 'ME_neval' + br_suff + '/D' )
    me_prob_br    = otree.Branch( 'ME_prob' + br_suff, me_prob, 'ME_prob' + br_suff + '/D' )
    jet_prob_br   = otree.Branch( 'JET_prob' + br_suff, jet_prob, 'JET_prob' + br_suff + '/D' )
    
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
    myint.setNEval( int(1.0e4), int(1.0e9) )
    
    # ---------------------------------------------------------------------

    #
    # configure the matrix element (integrand)
    #

    if not os.path.exists( 'cteq6l.tbl' ):
        log.debug( 'Fetching CTEQ6 PDF' )
        commands.getstatusoutput( 'ln -s $HEPPY/higgs/mymeanalysis/exec/run/cteq6l.tbl .' )

    if not os.path.exists( 'jetefficiency.tf' ):
        commands.getstatusoutput( 'ln -s /home/dschouten/code/higgs/mymeanalysis/exec/run/jetefficiency.tf .' )
        
    if 'tt' in matrix.lower():
        jetEffTF = root.JetEfficiencyTF( 'bjetefficiency.tf' )
        config += [ root.TT1j.kJET_3D ] 
        myme = matrixelement.toolFactory( matrix )( myint, args = config ) 
        integrand = myme.matrixElement
        integrand.setInvisibleJetTF( jetEffTF )
        integrand.setIntegrationLimits( integrand.getIPTB(), 1, 100 )
        integrand.setIntegrationLimits( integrand.getIPHIB(), -math.pi, math.pi )
        integrand.setIntegrationLimits( integrand.getIETAB(), -4, 4 )
        mymelist = [ myme ]
    else:
        config += [ root.WW1jLeptonsBaseIntegrand.kNEUTRINO_4D, None, None ]
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

        jet = root.TLorentzVector( itree.px_j, itree.py_j, itree.pz_j, itree.E_j )
        jet.SetPtEtaPhiE( jet.Pt(), jet.Eta(), jet.Phi(), jet.Rho() )
        
        met = root.TLorentzVector( (nul+nur).Px(), (nul+nur).Py(), 0, (nul+nur).Pt() )

        log.debug( '\t lepton pT:', lp.Pt(), lm.Pt() )
        log.debug( '\t jet pT:', jet.Pt() )
        log.debug( '\t missing ET:', met.Pt() )
        log.debug( '\t higgs mass:', (lp + nul + lm + nur).M() )

        result = rmap( )

        myme.matrixElement.setDynamicLimits( ) 

        buff = myme.evaluate( lp, lm, met, jet, doIntegration = True, coords = [ ] )
        print buff
        result += buff
        
        me[0]        = result['me']
        me_error[0]  = result['error']
        me_stat[0]   = result['status']
        me_neval[0]  = result['neval']
        me_prob[0]   = result['prob']

        log.debug( '\tme =', '%0.5g'%( result['me'] ) )
        log.debug( '\terror =', '%0.5g'%( result['error'] ) )
        log.debug( '\tstat =', '%0.5g'%( result['status'] ) )
        log.debug( '\tprob =', '%0.5g'%( result['prob'] ) )
        log.debug( '\tneval =', result['neval'] )

        if not iplace:
            otree.Fill( )
        else:
            me_br.Fill( )
            me_error_br.Fill()
            me_stat_br.Fill()
            me_neval_br.Fill()
            me_prob_br.Fill()

    if not iplace:
        otree.Write()
        ofile.Close()
    else:
        otree.Write( '', root.TObject.kOverwrite )

    log.info( 'Done processing events' )
