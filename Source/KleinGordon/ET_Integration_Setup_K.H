#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {E_U_0, E_U_1, E_U_2, A, A_L_0, A_L_1, A_L_2, Z, Phi1, Phi2, Pi1, Pi2, NumScalars};
};

namespace InitIdx { 
         enum ETInitIndexes {NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {AbsPhi, E_U_0, E_U_1, E_U_2, A, A_L_0, A_L_1, A_L_2, Z, Phi1, Phi2, Pi1, Pi2, NumScalars}; 
};

#define NUM_GHOST_CELLS 4

#endif