#include "OutputNtuple.h"

//Standard includes
#include <iostream>

//Standard namespaces
using namespace std;

//Matrix element namespaces
using namespace me::tools;

void OutputNtuple::map_tree_branches(TTree *tree)
{
	//Create branches in the tree
	tree->Branch( "me", &retval, "me/D" );
	tree->Branch( "error", &error, "error/D" );
	tree->Branch( "prob", &prob, "prob/D" );
	tree->Branch( "status", &fail, "status/I" );
	tree->Branch( "neval", &neval, "neval/I" );  
	tree->Branch( "time", &time, "time/D" );
}

void OutputNtuple::print()
{
	cout 
		<< "\tME: " << retval << " +/- " << error << endl
    	<< "\tStatus: " << (fail == 0 ? "Success" : "Failure") << endl
    	<< "\t# of Integrand Evaluations: " << neval << endl;
}
