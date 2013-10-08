#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TStopwatch.h"

gROOT->Reset();

void me_wgamma( unsigned BEGIN, unsigned END, const char* INPUT = "" ) {
  // -------------------------------------------------------------------------------

  const bool      USE_DIVONNE	= true;
  const unsigned  MAX_EVENTS	= 50000;
    
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

  if( USE_DIVONNE ){
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
  else{
    myint = new CubaIntegrators::Vegas();
  }
  myint->setNComp( 1 );
  myint->setEpsilon( 0.025, 0 );
  myint->setVerbose( 0 );
  myint->setSampleSet( true );
  myint->setRandomSeed( 1 );
  myint->setPseudoRandom( true );
  // Do not want timeouts since exception is not caught
  // So rely on max no of integrations instead
  myint->setTimeout( 6000 );
  myint->setNEval( int(1.0e4), int(1.0e7) );

  //  PhotonConversionTF *TF = new PhotonConversionTF("photonconversion.tf");
  Wgamma myWgammaME( myint, Wgamma::kBOOST, NULL );

  myWgammaME.setIntegrationLimits(myWgammaME.getIPZ(), -500, 500);
  //  myWgammaME.setIntegrationLimits(myWgammaME.getIBOOSTPT(), 0., 0.01);
  myWgammaME.setIntegrationLimits(myWgammaME.getIBOOSTPT(), 0., 30.);
  myWgammaME.setIntegrationLimits(myWgammaME.getIBOOSTPHI(), 0., 2*TMath::Pi());

  myWgammaME.initialize();

  // -------------------------------------------------------------------------------
  
  std::cout << std::endl;

  TFile *ifile = new TFile(INPUT, "UPDATE");
  TTree *itree = (TTree*)ifile->Get("HWWTree");

  Float_t lepPt0;
  Float_t lepEta0;
  Float_t lepPhi0;
  Float_t lepID0;

  Float_t lepPt1;
  Float_t lepEta1;
  Float_t lepPhi1;
  Float_t lepID1;

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
  itree->SetBranchStatus( "lepID0", true );
  itree->SetBranchAddress( "lepID0", &lepID0 );

  itree->SetBranchStatus( "lepPt1", true );
  itree->SetBranchAddress( "lepPt1", &lepPt1 );
  itree->SetBranchStatus( "lepEta1", true );
  itree->SetBranchAddress( "lepEta1", &lepEta1 );
  itree->SetBranchStatus( "lepPhi1", true );
  itree->SetBranchAddress( "lepPhi1", &lepPhi1 );
  itree->SetBranchStatus( "lepID1", true );
  itree->SetBranchAddress( "lepID1", &lepID1 );

  itree->SetBranchStatus( "m_jet_n", true );
  itree->SetBranchAddress( "m_jet_n", &jetN );
 
  itree->SetBranchStatus( "higgsPtEventWeight", true );
  itree->SetBranchAddress( "higgsPtEventWeight", &higgsPtEventWeight );
  itree->SetBranchStatus( "MVAEventWeight", true );
  itree->SetBranchAddress( "MVAEventWeight", &MVAEventWeight );

  itree->SetBranchStatus( "EventNumber", true );
  itree->SetBranchAddress( "EventNumber", &EventNumber );

  // -------------------------------------------------------------------------------

  ifile->Delete("Wgamma;*");
  TTree* otree = new TTree("Wgamma", "Wgamma ME");
    
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
    END = std::min( END, BEGIN + MAX_EVENTS );
  else
    END = MAX_EVENTS;

  for( unsigned int ievt = BEGIN; ievt < std::min( END, (unsigned)(itree->GetEntries()) ); ++ievt ){
    Int_t nb = itree->GetEntry( ievt );
    if( nb <= 0 ){
      std::cout << "ERROR retrieving event #" << ievt << std::endl;
      break;
    }

    if( jetN != 0 )
      continue;

    float lepMass = 0.;
    
    timer.Start( true );
    
    // What decay channel is the event?
    TString channel = "";
    if(fabs(lepID0)==11 && fabs(lepID1)==11)
      channel = "ee";
    else if(fabs(lepID0)==13 && fabs(lepID1)==13)
      channel = "mumu";
    else if((fabs(lepID0)==11 && fabs(lepID1)==13) ||
	    (fabs(lepID0)==13 && fabs(lepID1)==11))
      channel = "emu";
    else{
      std::cout << "ERROR Unrecognised decay channel in event #" << ievt << std::endl;
      continue;
    }

    if(channel == "mumu"){
      retval=1.0e-30; error=0.; fail=0.; neval=0.; prob=0.;
      std::cout << "INFO no ME calculated for mumu event #" << ievt << std::endl << std::endl;
      timer.Stop();
      time = timer.CpuTime();
      otree->Fill( );
      continue;
    }
    else if(channel == "ee"){
      if(lepPt0 > lepPt1){
	myWgammaME.measured.l.SetPtEtaPhiM( lepPt0/1000., lepEta0, lepPhi0, lepMass);
        myWgammaME.measured.ph.SetPtEtaPhiM(lepPt1/1000., lepEta1, lepPhi1, lepMass);
        myWgammaME.measured.lCharge = lepID0 > 0 ? -1. : +1.;
      }
      else{
	myWgammaME.measured.l.SetPtEtaPhiM( lepPt1/1000., lepEta1, lepPhi1, lepMass);
        myWgammaME.measured.ph.SetPtEtaPhiM(lepPt0/1000., lepEta0, lepPhi0, lepMass);
        myWgammaME.measured.lCharge = lepID1 > 0 ? -1. : +1.;
      }
    }
    else if(channel == "emu"){
      if(fabs(lepID0)==11){
	myWgammaME.measured.l.SetPtEtaPhiM( lepPt0/1000., lepEta0, lepPhi0, lepMass);
	myWgammaME.measured.ph.SetPtEtaPhiM(lepPt1/1000., lepEta1, lepPhi1, lepMass);
	myWgammaME.measured.lCharge = lepID0 > 0 ? -1. : +1.;
      }
      else{
	myWgammaME.measured.l.SetPtEtaPhiM( lepPt1/1000., lepEta1, lepPhi1, lepMass);
	myWgammaME.measured.ph.SetPtEtaPhiM(lepPt0/1000., lepEta0, lepPhi0, lepMass);
	myWgammaME.measured.lCharge = lepID1 > 0 ? -1. : +1.;
      }
    }
    myWgammaME.measured.nu.SetPx(-(myWgammaME.measured.l + myWgammaME.measured.ph).Px());
    myWgammaME.measured.nu.SetPy(-(myWgammaME.measured.l + myWgammaME.measured.ph).Py());

    std::cout << std::endl << "INFO dumping event kinematics for event #" << ievt << std::endl;
    myWgammaME.measured.print();
    myWgammaME.WmassConstraint();
    std::cout << std::endl << "INFO calculating ME for event #" << ievt << std::endl;
    myWgammaME.integrator()->doIntegral(&retval, &error, &fail, &neval, &prob);    
    std::cout << std::endl;

    timer.Stop( );
    std::cout << "INFO \t ME = " << retval << " +/- " << error << std::endl;
    std::cout << "INFO \t status = " << fail << ", # of integrand evaluations = " << neval << std::endl;

    time = timer.CpuTime();

    otree->Fill( );
  }

  // -------------------------------------------------------------------------------

  otree->Write( );
  ifile->Close( );
}
