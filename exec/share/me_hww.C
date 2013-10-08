#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TStopwatch.h"

gROOT->Reset();

void me_hww( unsigned BEGIN, unsigned END, 
	     const char* INPUT = "", const char* OUTPUT = "", double HMASS = 120.0 ) {
  // -------------------------------------------------------------------------------

  const bool      USE_DIVONNE	= true;
  const unsigned  MAX_EVENTS	= 50000;
    
  // from HDECAY ... 

  std::map<int,double> HWIDTH;
  HWIDTH[120]	=	0.3779E-02;
  HWIDTH[121]	=	0.3875E-02;
  HWIDTH[122]	=	0.3979E-02;
  HWIDTH[123]	=	0.4091E-02;
  HWIDTH[124]	=	0.4212E-02;
  HWIDTH[125]	=	0.4342E-02;
  HWIDTH[126]	=	0.4483E-02;
  HWIDTH[127]	=	0.4636E-02;
  HWIDTH[128]	=	0.4801E-02;
  HWIDTH[129]	=	0.4980E-02;
  HWIDTH[130]	=	0.5175E-02;
  HWIDTH[131]	=	0.5385E-02;
  HWIDTH[132]	=	0.5615E-02;
  HWIDTH[133]	=	0.5864E-02;
  HWIDTH[134]	=	0.6135E-02;
  HWIDTH[135]	=	0.6432E-02;
  HWIDTH[136]	=	0.6755E-02;
  HWIDTH[137]	=	0.7108E-02;
  HWIDTH[138]	=	0.7495E-02;
  HWIDTH[139]	=	0.7920E-02;
  HWIDTH[140]	=	0.8387E-02;

  // -------------------------------------------------------------------------------

  std::vector<std::string> fileList;

  if( std::string( INPUT ) == "" )
  {
    std::string args;
    std::ifstream ifs("inputs.txt");
    std::getline(ifs,args);
    
    for (size_t i=0,n; i <= args.length(); i=n+1)
    {
      n = args.find_first_of(',',i);
      if (n == string::npos)
	n = args.length();
      string tmp = args.substr(i,n-i);
      fileList.push_back( tmp );
    }
  }
  else
  {
    fileList.push_back( std::string( INPUT ) );
  }
  
  // -------------------------------------------------------------------------------

  gSystem->Load( "libintegrator.so" );
  gSystem->Load( "libintegrator_dict.so" );
  gSystem->Load( "libdhelas.so" );
  gSystem->Load( "libmatrix.so" );
  gSystem->Load( "libmatrix_dict.so" );
  gSystem->Load( "libnr.so" );
  gSystem->Load( "libttanalysis.so" );
  gSystem->Load( "libcuba.so" );

  // -------------------------------------------------------------------------------

  Integrator* myint = NULL;

  if( USE_DIVONNE )
  {
    myint = new CubaIntegrators::Divonne();

    CubaIntegrators::Divonne* int_ptr = dynamic_cast<CubaIntegrators::Divonne*>( myint );

    int_ptr->setPartitioningRule( 47 );
    int_ptr->setIntegrationRule( 1 );
    int_ptr->setRefinementRule( 1 );
    int_ptr->setMaxPass( 7 );
    int_ptr->setBorder( 0 );
    int_ptr->setChiSqr( 10.0 );
    int_ptr->setMinDev( 0.25 );  
  }
  else
  {
    myint = new CubaIntegrators::Vegas();
  }
  myint->setNComp( 1 );
  myint->setEpsilon( 0.025, 0 );
  myint->setVerbose( 0 );
  myint->setSampleSet( true );
  myint->setRandomSeed( 1 );
  myint->setPseudoRandom( true );
  myint->setNEval( int(1.0e4), int(1.0e8) );

  HWW myHWWME( myint, HMASS, HWW::kNEUTRINO_4D, 0x0 );
  myHWWME.setUseSM( false );
  
  myHWWME.setWidth( HWIDTH[(int)HMASS] / 0.1 );

  myHWWME.setIntegrationLimits( myHWWME.getIPTA(), 0, 400 );
  myHWWME.setIntegrationLimits( myHWWME.getIPHIA(), -acos(-1), acos(-1) );
  myHWWME.setIntegrationLimits( myHWWME.getIPZA(), -500, 500 );
  myHWWME.setIntegrationLimits( myHWWME.getIPZB(), -500, 500 );

  myHWWME.initialize();

  // -------------------------------------------------------------------------------
  
  std::cout << std::endl;
  TChain* itree = new TChain( "HWWTree", "HWWTree" );
  for (int iFile=0; iFile<fileList.size(); ++iFile)
  {
    std::cout << "INFO opening " << fileList[iFile].c_str() << std::endl;
    itree->Add( fileList[iFile].c_str() );
  }

  Float_t lepPt0;
  Float_t lepEta0;
  Float_t lepPhi0;
  Float_t lepCharge0;

  Float_t lepPt1;
  Float_t lepEta1;
  Float_t lepPhi1;
  Float_t lepCharge1;

  Int_t jetN;

  Double_t higgsPtEventWeight;
  Double_t MVAEventWeight;

  UInt_t EventNumber;

  itree->SetBranchStatus( "*", false ); 
  
  itree->SetBranchStatus( "lepPt0", true );
  itree->SetBranchAddress( "lepPt0", &lepPt0 );
  itree->SetBranchStatus( "lepEta0", true );
  itree->SetBranchAddress( "lepEta0", &lepEta0 );
  itree->SetBranchStatus( "lepPhi0", true );
  itree->SetBranchAddress( "lepPhi0", &lepPhi0 );
  itree->SetBranchStatus( "lepCharge0", true );
  itree->SetBranchAddress( "lepCharge0", &lepCharge0 );

  itree->SetBranchStatus( "lepPt1", true );
  itree->SetBranchAddress( "lepPt1", &lepPt1 );
  itree->SetBranchStatus( "lepEta1", true );
  itree->SetBranchAddress( "lepEta1", &lepEta1 );
  itree->SetBranchStatus( "lepPhi1", true );
  itree->SetBranchAddress( "lepPhi1", &lepPhi1 );
  itree->SetBranchStatus( "lepCharge1", true );
  itree->SetBranchAddress( "lepCharge1", &lepCharge1 );

  itree->SetBranchStatus( "m_jet_n", true );
  itree->SetBranchAddress( "m_jet_n", &jetN );
 
  itree->SetBranchStatus( "higgsPtEventWeight", true );
  itree->SetBranchAddress( "higgsPtEventWeight", &higgsPtEventWeight );

  itree->SetBranchStatus( "MVAEventWeight", true );
  itree->SetBranchAddress( "MVAEventWeight", &MVAEventWeight );

  itree->SetBranchStatus( "EventNumber", true );
  itree->SetBranchAddress( "EventNumber", &EventNumber );

  // -------------------------------------------------------------------------------

  TFile* ofile = NULL;
  if( std::string( OUTPUT ) != "" )
  {
    ofile = TFile::Open( OUTPUT, "recreate" );
  }
  else
  {
    ofile = TFile::Open( "me.root" , "recreate" );
  }

  TTree* otree = itree->CloneTree( 0, "fast SortBasketsByEntry" );

  stringstream buff;

  buff << "HWW_ME_" << (int)HMASS;
  
  otree->SetName( buff.str().c_str() );

  buff.flush();
  
  buff << "HWW ME (mH = " << (int)HMASS << " GeV)";

  otree->SetTitle( buff.str().c_str() );
    
  Double_t retval, error, prob, time;
  Int_t fail, neval;
  
  otree->Branch( "me", &retval, "me/D" );
  otree->Branch( "error", &error, "error/D" );
  otree->Branch( "prob", &prob, "prob/D" );
  otree->Branch( "status", &fail, "status/I" );
  otree->Branch( "neval", &neval, "neval/I" );  
  otree->Branch( "time", &time, "time/D" );

  // -------------------------------------------------------------------------------

  TStopwatch timer;
  
  if( END != 0 )
  {
    END = std::min( END, BEGIN + MAX_EVENTS );
  }
  else
  {
    END = MAX_EVENTS;
  }

  for( unsigned int ievt = BEGIN; ievt < std::min( END, (unsigned)(itree->GetEntries()) ); ++ievt )
  {
    Int_t nb = itree->GetEntry( ievt );
    if( nb <= 0 )
    {
      std::cout << "ERROR retrieving event #" << ievt << std::endl;
      break;
    }

    if( jetN != 0 )
      continue;

    float lepMass;

    if( lepCharge0 > 0 )
    {
      myHWWME.measured.lp.SetPtEtaPhiM( lepPt0 / 1000., lepEta0, lepPhi0, lepMass );
      myHWWME.measured.lm.SetPtEtaPhiM( lepPt1 / 1000., lepEta1, lepPhi1, lepMass );
    }
    else
    {
      myHWWME.measured.lm.SetPtEtaPhiM( lepPt0 / 1000., lepEta0, lepPhi0, lepMass );
      myHWWME.measured.lp.SetPtEtaPhiM( lepPt1 / 1000., lepEta1, lepPhi1, lepMass );
    }
      
    myHWWME.measured.met.SetPx( -(myHWWME.measured.lp + myHWWME.measured.lm).Px() );
    myHWWME.measured.met.SetPy( -(myHWWME.measured.lp + myHWWME.measured.lm).Py() );
    
    
    std::cout << std::endl << "INFO dumping event kinematics" << std::endl;
    myHWWME.measured.print();
        
    timer.Start( true );

    std::cout << std::endl << "INFO calculating ME for event #" << ievt << std::endl;
    myHWWME.integrator()->doIntegral( &retval, &error, &fail, &neval, &prob );

    timer.Stop( );

    std::cout << "INFO \t ME = " << retval << " +/- " << error << std::endl;
    std::cout << "INFO \t status = " << fail << ", # of integrand evaluations = " << neval << std::endl;

    time = timer.CpuTime();

    otree->Fill( );
  }

  // -------------------------------------------------------------------------------

  otree->Write( );
  ofile->Close( );
}
