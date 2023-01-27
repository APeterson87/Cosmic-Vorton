#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {u, w_U_0, w_U_1, w_U_2, chi, phi1, phi2, sigma1, sigma2, dDtsigma1, dDtsigma2, A_L_0, A_L_1, A_L_2, V, A, B, C, S_U_phi, NumScalars}; 
};

namespace InitIdx { 
         enum ETInitIndexes {NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {w_U_phi, NumScalars}; 
};

#define NUM_GHOST_CELLS 2

#endif