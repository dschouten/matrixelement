#include "Higgs.h"

//Standard includes
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//Standard namespaces
using namespace std;

//Boost includes
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

//Boost namespaces
using namespace boost;

double me::tools::calculate_higgs_width(double higgs_mass)
{
	//Open the higgs width file
	ifstream input_file("hwidth.txt");
	if(!input_file.is_open())
	{
		cerr << "ERROR: Could not open higgs width file" << endl;
		return 1.0;
	}

	//Set the initial value
	double result = 1.0;

	//Read through the lines of the file, parsing each one
	string line;
	while(getline(input_file, line))
	{
		//Continue if the line is empty
		if(line.length() == 0)
		{
			continue;
		}
		
		//Split the data
		vector<string> columns;
		split(columns, line, is_any_of("\t"));

		//Convert data
		double mass = boost::lexical_cast<double>(columns[0]);
		double width = boost::lexical_cast<double>(columns[1]);

		//Check if we've found it
		if(mass == higgs_mass)
		{
			result = width;
			break;
		}
	}

	//Close the file
	input_file.close();

	return result;
}
