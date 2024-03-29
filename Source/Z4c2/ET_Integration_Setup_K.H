#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {chi, gamtilde_LL_00, gamtilde_LL_01, gamtilde_LL_02, gamtilde_LL_11, gamtilde_LL_12, gamtilde_LL_22, Khat, Atilde_LL_00, Atilde_LL_01, Atilde_LL_02, Atilde_LL_11, Atilde_LL_12, Atilde_LL_22, theta, Gamtilde_U_0, Gamtilde_U_1, Gamtilde_U_2, alpha, beta_U_0, beta_U_1, beta_U_2, B_U_0, B_U_1, B_U_2, E_U_0, E_U_1, E_U_2, A, A_L_0, A_L_1, A_L_2, Z, phi1, phi2, Pi1, Pi2, NumScalars}; 
};

namespace InitIdx { 
         enum ETInitIndexes {u, alphapsi, chi, V, A, B, NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {rhoADM, Psi4re, Psi4im, ModPsi4, NumScalars};  
         //enum DiagnosticIndexes {Psi4re, Psi4im, ModPsi4, AbsAlpha, Ksclr, NumScalars};
};

#define NUM_GHOST_CELLS 4

#endif