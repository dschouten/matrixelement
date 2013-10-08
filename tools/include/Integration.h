#pragma once

//Standard includes
#include <string>

//Boost includes
#include <boost/shared_ptr.hpp>
#include <boost/variant.hpp>

//Integrator includes
#include "integrator/Integrator.hh"

//Matrix includes
#include "matrix/FeynIntegrand_DYFeynDiagram.hh"
#include "matrix/FeynIntegrand_WW0jFeynDiagram.hh"
#include "matrix/FeynIntegrand_WW1jFeynDiagram.hh"
#include "matrix/FeynIntegrand_WW2jFeynDiagram.hh"
#include "matrix/DY.hh"
#include "matrix/WW.hh"
#include "matrix/WW1j.hh"
#include "matrix/HWW.hh"
#include "matrix/HWW1j.hh"
#include "matrix/TT1j.hh"

//Tools includes
#include "EventNtuple.h"
#include "OutputNtuple.h"

namespace me
{
	namespace tools
	{
		/* Integration Type Definitions */
		enum
		{
			IntegratorTypeCubaDivonne,
			IntegratorTypeCubaVegas,
			IntegratorTypeGSLVegas,
			IntegratorTypeUnknown
		};
		typedef int IntegratorType;

		/* Integrator Construction Routines */
		IntegratorType integrator_type_by_name(std::string name);
		Integrator * create_integrator(IntegratorType type);

		/* Matrix Element Type Definitions */
		enum
		{
			MatrixElementTypeWW,
			MatrixElementTypeWW1j,
			MatrixElementTypeHWW,
			MatrixElementTypeHWW1j,
			MatrixElementTypeTT1j,
			MatrixElementTypeDY,
			MatrixElementTypeUnknown
		};
		typedef int MatrixElementType;
		typedef boost::shared_ptr<WW> WWRef;
		typedef boost::shared_ptr<WW1j> WW1jRef;
		typedef boost::shared_ptr<HWW> HWWRef;
		typedef boost::shared_ptr<HWW1j> HWW1jRef;
		typedef boost::shared_ptr<TT1j> TT1jRef;
		typedef boost::shared_ptr<DY> DYRef;
		typedef boost::variant<WWRef, WW1jRef, HWWRef, HWW1jRef, TT1jRef, DYRef> MatrixIntegrand;

		/* Matrix Element Construction Routines */
		MatrixElementType matrix_element_type_by_name(std::string name);
		MatrixIntegrand create_matrix_element(MatrixElementType type, std::string me_options, Integrator *integrator);

		/* Matrix Element Evaluation Objects */
		class IntegrationExecuter : public boost::static_visitor<void>
		{
			public:
				/* Construction/Destruction */
				IntegrationExecuter(EventNtuple *e, OutputNtuple *o);

				/* Evaluation Callbacks */
				void operator()(WWRef integrand);
				void operator()(WW1jRef integrand);
				void operator()(HWWRef integrand);
				void operator()(HWW1jRef integrand);
				void operator()(TT1jRef integrand);
				void operator()(DYRef integrand);

			private:
				/* Private Data Members */
				EventNtuple *event;
				OutputNtuple *output;
		};
	}
}
