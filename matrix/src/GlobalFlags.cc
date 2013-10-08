//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/GlobalFlags.hh"

#ifdef DEBUG
bool GlobalFlags::debug = true;
#else
bool GlobalFlags::debug = false;
#endif

long GlobalFlags::max_debug = 5;
