#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {w_U_0, w_U_1, w_U_2, phi1, phi2, sigma1, sigma2, test, NumScalars};
};

namespace InitIdx { 
         enum ETInitIndexes {NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {NumScalars}; 
};

#define NUM_GHOST_CELLS 2

#endif