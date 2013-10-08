//Standard includes
#include <iostream>
#include <string>

//Boost includes
#include <boost/program_options.hpp>

//ROOT includes
#include <TChain.h>
#include <TFile.h>
#include <TStopwatch.h>

//Tools includes
#include "EventNtuple.h"
#include "OutputNtuple.h"
#include "Integration.h"

//Standard namespaces
using namespace std;

//Boost namespace aliases
namespace po = boost::program_options;

//Boost namespaces
using namespace boost;

//Tools namespaces
using namespace me::tools;

//Helper functions
po::variables_map parse_command_line_options(int argc, char * argv[])
{
	//Create the options specifier.  This code is correct, it's just
	//Boost using the Boost.Assign voodoo.
	po::options_description desc("Allowed Options (non-global options are combinatorially combined for form different pipelines)");
	desc.add_options()
		("input,i", po::value<string>()->default_value("input.root"), "The path (either a file or directory) to the input data.")
		("output,o", po::value<string>()->default_value("output.root"), "The path to the output ntuple.")
		("integrator,I", po::value<string>()->default_value("cdivonne"), "The integrator to use (can be cdivonne (cuba Divonne), cvegas (cuba Vegas), or gvegas (GSL Vegas)).")
		("matrix,m", po::value<string>()->default_value("WW"), "The name of the matrix element to evaluate.")
		("options,O", po::value<string>()->default_value(""), "The options to pass into matrix element construction.")
		("verbose,v", "Make the program print more detailed output to command line (global option).")
		("help,h", "Print a description of the program options.")
	;
	
	//Do the actual parsing
	po::variables_map vm;
	try 
	{
		po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
		po::notify(vm);
		
		//Print help if necessary
		if(vm.count("help"))
		{
			cout << desc << endl;
			exit(1);
		}
	}
	catch(std::exception& e)
	{
		cerr << "Couldn't parse command line options: " << e.what() << endl;
		cout << desc << endl;
		exit(1);
	}
	catch(...)
	{
		cerr << "Couldn't parse command line options, not sure why." << endl;
		cout << desc << endl;
		exit(1);
	}
	
	return vm;
}

TChain * create_input_tree(string input_path)
{
	//Create the tree
	TChain *tree = new TChain("HWWTree", "HWWTree");

	//Add files
	tree->Add(input_path.c_str());

	return tree;
}

int main(int argc, char * argv[])
{
	//Print program information
	cout << "Matrix Element Calculator" << endl;

	//Parse command line options
	po::variables_map options = parse_command_line_options(argc, argv);

	//Create an integrator
	IntegratorType integrator_type = integrator_type_by_name(options["integrator"].as<string>());
	Integrator *integrator = create_integrator(integrator_type);
	if(integrator == NULL)
	{
		exit(1);
	}

	//Create the matrix element integrand
	MatrixElementType matrix_element_type = matrix_element_type_by_name(options["matrix"].as<string>());
	MatrixIntegrand integrand = create_matrix_element(matrix_element_type, 
		options["options"].as<string>(), 
		integrator);

	//Open the input tree
	TChain* input_tree = create_input_tree(options["input"].as<string>());
	if(input_tree == NULL)
	{
		exit(1);
	}

	//Create the structure that we'll map the tree into
	EventNtuple *event = new EventNtuple;
	event->map_tree_branches(input_tree);

	//Create the output file/tree
	TFile *output_file = TFile::Open(options["output"].as<string>().c_str(), "recreate");
	if(output_file == NULL)
	{
		cerr << "Unable to open output file" << endl;
		exit(1);
	}
	TTree* output_tree = input_tree->CloneTree(0, "fast SortBasketsByEntry");
	if(output_tree == NULL)
	{
		cerr << "Unable to create output tree" << endl;
		exit(1);
	}
	output_tree->SetName(options["matrix"].as<string>().c_str());
	output_tree->SetTitle((options["matrix"].as<string>() + "Matrix Element").c_str());

	//Create the struct that the output data will go into
	OutputNtuple *output = new OutputNtuple;
	output->map_tree_branches(output_tree);

	//Create the integration processor
	IntegrationExecuter executer(event, output);

	//Now loop over the tree
	Int_t n_events = input_tree->GetEntries();
	cout << "There are " << n_events << " events" << endl;
	for(Int_t i_event = 0; i_event < n_events; i_event++)
	{
		//Load the event into the event structure
		Int_t status = input_tree->GetEntry(i_event);
		if(status <= 0)
		{
			cerr << "Error retrieving event " << i_event << endl;
			exit(1);
		}
		cout << endl << "Event " << i_event << endl;

		//Do the integration
		apply_visitor(executer, integrand);

		//Print the output
		output->print();

		//Store the output
		output_tree->Fill();
	}

	//Write data
	output_tree->Write();
	output_file->Close();

	//Cleanup
	delete output;
	delete output_file;
	delete event;
	delete input_tree;
	delete integrator;

	//Successfully finished
	cout << "Finished" << endl;

	return 0;
}
