#pragma once

//ROOT includes
#include <TTree.h>

namespace me
{
	namespace tools
	{
		struct OutputNtuple
		{
			//Maps the tree branches for the output
			//data to the appropriate variables
			//as described below.
			void map_tree_branches(TTree *tree);

			//Print information on the integration
			void print();

			
			/* Integration Information */

			Double_t retval;
			Double_t error;
			Double_t prob;
			Double_t time;
			Int_t fail;
			Int_t neval;
		};
	}
}
