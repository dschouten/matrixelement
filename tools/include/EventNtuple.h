#pragma once

//ROOT includes
#include <TTree.h>

namespace me
{
	namespace tools
	{
		struct EventNtuple
		{
			//Maps the tree branches for the event
			//kinematics to the appropriate variables
			//as described below.
			void map_tree_branches(TTree *tree);

			
			/* Lepton 0 Information */

			Float_t lepPt0;
			Float_t lepEta0;
			Float_t lepPhi0;
			Float_t lepID0;


			/* Lepton 1 Information */

			Float_t lepPt1;
			Float_t lepEta1;
			Float_t lepPhi1;
			Float_t lepID1;


			/* Jet Information */

			Int_t jetN;


			/* Event Information */

			Double_t higgsPtEventWeight;
			Double_t MVAEventWeight;
			UInt_t EventNumber;
		};
	}
}
