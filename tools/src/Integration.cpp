#include "Integration.h"

//Standard includes
#include <map>

//Boost includes
#include <boost/lexical_cast.hpp>

//ROOT includes
#include <TDatabasePDG.h>

//Integrator includes
#include "integrator/CubaIntegrator.hh"

//Tools includes
#include "Higgs.h"

//Matrix includes
#include "matrix/TauTF.hh"

//Standard namespaces
using namespace std;

//Boost namespaces
using namespace boost;

//Integrator namespaces
using namespace CubaIntegrators;

//Matrix element namespaces
using namespace me::tools;

IntegratorType me::tools::integrator_type_by_name(string name)
{
	if(name == "cdivonne")
	{
		return IntegratorTypeCubaDivonne;
	}
	else if(name == "cvegas")
	{
		return IntegratorTypeCubaVegas;
	}
	else if(name == "gvegas")
	{
		return IntegratorTypeGSLVegas;
	}
	
	return IntegratorTypeUnknown;
}

Integrator * me::tools::create_integrator(IntegratorType type)
{
	//Create the initial result
	Integrator *integrator = NULL;

	//Switch based on the type of integrator requested
	if(type == IntegratorTypeCubaDivonne)
	{
		//Create a cuba Divonne integrator and set parameters
		Divonne *divonne = new Divonne();
		divonne->setPartitioningRule(47);
		divonne->setIntegrationRule(1);
		divonne->setRefinementRule(1);
		divonne->setMaxPass(7);
		divonne->setBorder(0);
		divonne->setChiSqr(10.0);
		divonne->setMinDev(0.25);
		integrator = divonne;
	}
	else if(type == IntegratorTypeCubaVegas)
	{
		//Create a cuba Vegas integrator
		Vegas *vegas = new Vegas();
		integrator = vegas;
	}
	else if(type == IntegratorTypeGSLVegas)
	{
		//Create a GSL Vegas integrator
		//Not yet implemented
		cerr << "ERROR: GSL Vegas integrator is not yet implemented" << endl;
		return NULL;
	}
	else
	{
		//Unknown integrator type
		cerr << "Unknown integrator type" << endl;
		return NULL;
	}

	//Set common parameters
	integrator->setNComp(1);
	integrator->setEpsilon(0.025, 0);
	integrator->setVerbose(0);
	integrator->setSampleSet(true);
	integrator->setRandomSeed(1);
	integrator->setPseudoRandom(true);
	integrator->setNEval(int(1.0e4), int(1.0e8));

	//Return the result
	return integrator;
}

MatrixElementType me::tools::matrix_element_type_by_name(string name)
{
	if(name == "WW")
	{
		return MatrixElementTypeWW;
	}
	else if(name == "WW1j")
	{
		return MatrixElementTypeWW1j;
	}
	else if(name == "HWW")
	{
		return MatrixElementTypeHWW;
	}
	else if(name == "HWW1j")
	{
		return MatrixElementTypeHWW1j;
	}
	else if(name == "TT1j")
	{
		return MatrixElementTypeTT1j;
	}
	else if(name == "DY")
	{
		return MatrixElementTypeDY;
	}

	return MatrixElementTypeUnknown;
}

MatrixIntegrand me::tools::create_matrix_element(MatrixElementType type, string me_options, Integrator *integrator)
{
	//Switch based on the type of matrix element we want to create
	if(type == MatrixElementTypeWW)
	{
		//WW

		//Create a tau transfer function if requested
		TransferFunction *tau_tf = NULL;
		if(me_options.find("tau") != string::npos)
		{
			tau_tf = new TauTF;
		}

		//Create the matrix element integrand
		WWRef ww(new WW(integrator, WW::kNEUTRINO_4D, NULL, tau_tf));
		ww->setIntegrationLimits(ww->getIPTA(), 0, 250);
		ww->setIntegrationLimits(ww->getIPHIA(), -acos(-1), acos(-1));
		ww->setIntegrationLimits(ww->getIPZA(), -500, 500);
		ww->setIntegrationLimits(ww->getIPZB(), -500, 500);
		if(tau_tf != NULL)
		{
			ww->setIntegrationLimits(ww->getITAUE(), 0, 500);
		}
		ww->initialize();

		return ww;
	}
	else if(type == MatrixElementTypeWW1j)
	{

	}
	else if(type == MatrixElementTypeHWW)
	{
		//HWW

		//Parse the higgs mass
		double higgs_mass = 125.0;
		if(me_options != "")
		{
			higgs_mass = lexical_cast<double>(me_options);
		}
		double higgs_width = calculate_higgs_width(higgs_mass);
		cout << "Using a higgs mass of " << higgs_mass << " GeV with width " << higgs_width << endl;

		//Create the matrix element integrand
		HWWRef hww(new HWW(integrator, higgs_mass, HWW::kNEUTRINO_4D, NULL,  NULL));
		higgs_width = higgs_width + (pow(higgs_width, 1.0/6.0) / 2.0);
		hww->setWidth(higgs_width);
		hww->setIntegrationLimits(hww->getIPTA(), 0, 250);
		hww->setIntegrationLimits(hww->getIPHIA(), -acos(-1), acos(-1));
		hww->setIntegrationLimits(hww->getIPZA(), -500, 500);
		hww->setIntegrationLimits(hww->getIPZB(), -500, 500);
		hww->initialize();

		return hww;
	}
	else if(type == MatrixElementTypeHWW1j)
	{

	}
	else if(type == MatrixElementTypeTT1j)
	{

	}
	else if(type == MatrixElementTypeDY)
	{

	}

	//Unknown type
	cerr << "Unknown matrix element type" << endl;
	exit(1);
}

IntegrationExecuter::IntegrationExecuter(EventNtuple *e, OutputNtuple *o) :
event(e),
output(o)
{

}

void IntegrationExecuter::operator()(WWRef integrand)
{
	//Create a particle database
	TDatabasePDG pdg_database;

	//Fill in kinematics
	if(event->lepID0 < 0)
	{
		integrand->measured.lp.SetPtEtaPhiM(event->lepPt0/1000.0, event->lepEta0, event->lepPhi0, pdg_database.GetParticle(abs(event->lepID0))->Mass());
		integrand->measured.lm.SetPtEtaPhiM(event->lepPt1/1000.0, event->lepEta1, event->lepPhi1, pdg_database.GetParticle(abs(event->lepID1))->Mass());
	}
	else
	{
		integrand->measured.lm.SetPtEtaPhiM(event->lepPt0/1000.0, event->lepEta0, event->lepPhi0, pdg_database.GetParticle(abs(event->lepID0))->Mass());
		integrand->measured.lp.SetPtEtaPhiM(event->lepPt1/1000.0, event->lepEta1, event->lepPhi1, pdg_database.GetParticle(abs(event->lepID1))->Mass());
	}

	integrand->measured.met.SetPx(-(integrand->measured.lp + integrand->measured.lm).Px());
	integrand->measured.met.SetPy(-(integrand->measured.lp + integrand->measured.lm).Py());

	//Start the integration timer
	TStopwatch timer;
	timer.Start(true);
	integrand->integrator()->doIntegral(&output->retval, &output->error, &output->fail, &output->neval, &output->prob);
	timer.Stop();
	output->time = timer.CpuTime();
}

void IntegrationExecuter::operator()(WW1jRef integrand)
{
	
}

void IntegrationExecuter::operator()(HWWRef integrand)
{
	//Create a particle database
	TDatabasePDG pdg_database;

	//Fill in kinematics
	if(event->lepID0 < 0)
	{
		integrand->measured.lp.SetPtEtaPhiM(event->lepPt0/1000.0, event->lepEta0, event->lepPhi0, pdg_database.GetParticle(abs(event->lepID0))->Mass());
		integrand->measured.lm.SetPtEtaPhiM(event->lepPt1/1000.0, event->lepEta1, event->lepPhi1, pdg_database.GetParticle(abs(event->lepID1))->Mass());
	}
	else
	{
		integrand->measured.lm.SetPtEtaPhiM(event->lepPt0/1000.0, event->lepEta0, event->lepPhi0, pdg_database.GetParticle(abs(event->lepID0))->Mass());
		integrand->measured.lp.SetPtEtaPhiM(event->lepPt1/1000.0, event->lepEta1, event->lepPhi1, pdg_database.GetParticle(abs(event->lepID1))->Mass());
	}

	integrand->measured.met.SetPx(-(integrand->measured.lp + integrand->measured.lm).Px());
	integrand->measured.met.SetPy(-(integrand->measured.lp + integrand->measured.lm).Py());

	//Start the integration timer
	TStopwatch timer;
	timer.Start(true);
	integrand->integrator()->doIntegral(&output->retval, &output->error, &output->fail, &output->neval, &output->prob);
	timer.Stop();
	output->time = timer.CpuTime();
}

void IntegrationExecuter::operator()(HWW1jRef integrand)
{
	
}

void IntegrationExecuter::operator()(TT1jRef integrand)
{
	
}

void IntegrationExecuter::operator()(DYRef integrand)
{
	
}
