
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

//////////////////////////////////////////////////////
//
// be careful to add all template specializations
// that are used in any derived classes
//
//////////////////////////////////////////////////////

#pragma link C++ namespace hepstd;
#pragma link C++ class FeynIntegrand<DY2jFeynDiagram>;

#endif
