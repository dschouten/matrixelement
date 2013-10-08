#include "EventNtuple.h"

//Matrix element namespaces
using namespace me::tools;

void EventNtuple::map_tree_branches(TTree *tree)
{
	//Disable all branches
	tree->SetBranchStatus("*", false); 

	//Map lepton 0 information
	tree->SetBranchStatus("lepPt0", true);
	tree->SetBranchAddress("lepPt0", &lepPt0);
	tree->SetBranchStatus("lepEta0", true);
	tree->SetBranchAddress("lepEta0", &lepEta0);
	tree->SetBranchStatus("lepPhi0", true);
	tree->SetBranchAddress("lepPhi0", &lepPhi0);
	tree->SetBranchStatus("lepID0", true);
	tree->SetBranchAddress("lepID0", &lepID0);

	//Map lepton 1 information
	tree->SetBranchStatus("lepPt1", true);
	tree->SetBranchAddress("lepPt1", &lepPt1);
	tree->SetBranchStatus("lepEta1", true);
	tree->SetBranchAddress("lepEta1", &lepEta1);
	tree->SetBranchStatus("lepPhi1", true);
	tree->SetBranchAddress("lepPhi1", &lepPhi1);
	tree->SetBranchStatus("lepID1", true);
	tree->SetBranchAddress("lepID1", &lepID1);

	//Map jet information
	tree->SetBranchStatus("m_jet_n", true);
	tree->SetBranchAddress("m_jet_n", &jetN);

	//Map event information
	tree->SetBranchStatus("higgsPtEventWeight", true);
	tree->SetBranchAddress("higgsPtEventWeight", &higgsPtEventWeight);
	tree->SetBranchStatus("MVAEventWeight", true);
	tree->SetBranchAddress("MVAEventWeight", &MVAEventWeight);
	tree->SetBranchStatus("EventNumber", true);
	tree->SetBranchAddress("EventNumber", &EventNumber);
}
