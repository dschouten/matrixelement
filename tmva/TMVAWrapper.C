
#include "TMVA/Reader.h"

#include "TString.h"
#include "TObjArray.h"

#include <vector>
#include <map>
#include <iostream>

#define MAX_NUM_READERS 10

namespace TMVATools
{
  TMVA::Reader* myReader[MAX_NUM_READERS] = { 0x0 };
  float p__data[100];

  std::string p__method[MAX_NUM_READERS] = { std::string("") };
  
  void loadReader( const char* method, const char* wts_file, const char* variables, unsigned ireader )
  {
    unsigned id = std::min( ireader, (unsigned)MAX_NUM_READERS );
    if( myReader[id] != 0x0 )
    {
      std::cerr << "WARNING overwriting existing TMVA::Reader object [" << id << "]" << std::endl;
      delete myReader[id];
    }
    myReader[id] = new TMVA::Reader( "!Color:!Silent" );
    TObjArray* list = TString( variables ).Tokenize( ":" );
    for( unsigned int il = 0; il < (unsigned) list->GetEntries(); ++il )
    {
      myReader[id]->AddVariable( ((TObjString*)list->At(il))->GetString().Data(), &( p__data[il] ) );
    }
    p__method[id] = method;
    myReader[id]->BookMVA( method, wts_file );
    delete list;
  }
  
  double evaluateRegressionMVA( unsigned id , double a = 0, double b = 0, double c = 0, 
				double d = 0, double e = 0, double f = 0, double g = 0, 
				double h = 0, double i = 0, double j = 0, double k = 0, 
				double l = 0, double m = 0, double n = 0, double o = 0 )
  {
    unsigned ic = 0;
    p__data[ic++] = float( a );    
    p__data[ic++] = float( b );
    p__data[ic++] = float( c );    
    p__data[ic++] = float( d );
    p__data[ic++] = float( e );    
    p__data[ic++] = float( f );
    p__data[ic++] = float( g );    
    p__data[ic++] = float( h );
    p__data[ic++] = float( i );    
    p__data[ic++] = float( j );
    p__data[ic++] = float( k );    
    p__data[ic++] = float( l );
    p__data[ic++] = float( m );    
    p__data[ic++] = float( n );
    p__data[ic++] = float( o );   

    // for( ic = 0; ic < 15; ++ic ) { std::cout << p__data[ic] << " "; }
    // std::cout << endl;

    return (myReader[id]->EvaluateRegression( p__method[id].c_str() ))[0];
  }
}
