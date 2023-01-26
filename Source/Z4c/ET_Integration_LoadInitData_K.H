#ifndef ET_INTEGRATION_LOADINITDATA_K_H
#define ET_INTEGRATION_LOADINITDATA_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
initialize_from_data(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab, 
        amrex::Array4<amrex::Real const> const& initial_data,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x_U_0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x_U_1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x_U_2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        
        amrex::Real dDw_UL_00 = ((2.0/3.0)*initial_data(i + 1, j, k, InitIdx::w_U_0) - 1.0/12.0*initial_data(i + 2, j, k, InitIdx::w_U_0) - 2.0/3.0*initial_data(i - 1, j, k, InitIdx::w_U_0) + (1.0/12.0)*initial_data(i - 2, j, k, InitIdx::w_U_0))/dx[0];
        amrex::Real dDw_UL_01 = ((2.0/3.0)*initial_data(i, j + 1, k, InitIdx::w_U_0) - 1.0/12.0*initial_data(i, j + 2, k, InitIdx::w_U_0) - 2.0/3.0*initial_data(i, j - 1, k, InitIdx::w_U_0) + (1.0/12.0)*initial_data(i, j - 2, k, InitIdx::w_U_0))/dx[1];
        amrex::Real dDw_UL_02 = ((2.0/3.0)*initial_data(i, j, k + 1, InitIdx::w_U_0) - 1.0/12.0*initial_data(i, j, k + 2, InitIdx::w_U_0) - 2.0/3.0*initial_data(i, j, k - 1, InitIdx::w_U_0) + (1.0/12.0)*initial_data(i, j, k - 2, InitIdx::w_U_0))/dx[2];
        amrex::Real dDw_UL_10 = ((2.0/3.0)*initial_data(i + 1, j, k, InitIdx::w_U_1) - 1.0/12.0*initial_data(i + 2, j, k, InitIdx::w_U_1) - 2.0/3.0*initial_data(i - 1, j, k, InitIdx::w_U_1) + (1.0/12.0)*initial_data(i - 2, j, k, InitIdx::w_U_1))/dx[0];
        amrex::Real dDw_UL_11 = ((2.0/3.0)*initial_data(i, j + 1, k, InitIdx::w_U_1) - 1.0/12.0*initial_data(i, j + 2, k, InitIdx::w_U_1) - 2.0/3.0*initial_data(i, j - 1, k, InitIdx::w_U_1) + (1.0/12.0)*initial_data(i, j - 2, k, InitIdx::w_U_1))/dx[1];
        amrex::Real dDw_UL_12 = ((2.0/3.0)*initial_data(i, j, k + 1, InitIdx::w_U_1) - 1.0/12.0*initial_data(i, j, k + 2, InitIdx::w_U_1) - 2.0/3.0*initial_data(i, j, k - 1, InitIdx::w_U_1) + (1.0/12.0)*initial_data(i, j, k - 2, InitIdx::w_U_1))/dx[2];
        amrex::Real dDw_UL_20 = ((2.0/3.0)*initial_data(i + 1, j, k, InitIdx::w_U_2) - 1.0/12.0*initial_data(i + 2, j, k, InitIdx::w_U_2) - 2.0/3.0*initial_data(i - 1, j, k, InitIdx::w_U_2) + (1.0/12.0)*initial_data(i - 2, j, k, InitIdx::w_U_2))/dx[0];
        amrex::Real dDw_UL_21 = ((2.0/3.0)*initial_data(i, j + 1, k, InitIdx::w_U_2) - 1.0/12.0*initial_data(i, j + 2, k, InitIdx::w_U_2) - 2.0/3.0*initial_data(i, j - 1, k, InitIdx::w_U_2) + (1.0/12.0)*initial_data(i, j - 2, k, InitIdx::w_U_2))/dx[1];
        amrex::Real dDw_UL_22 = ((2.0/3.0)*initial_data(i, j, k + 1, InitIdx::w_U_2) - 1.0/12.0*initial_data(i, j, k + 2, InitIdx::w_U_2) - 2.0/3.0*initial_data(i, j, k - 1, InitIdx::w_U_2) + (1.0/12.0)*initial_data(i, j, k - 2, InitIdx::w_U_2))/dx[2];
        
        amrex::Real Abar_UU_00 = 2*dDw_UL_00;
        amrex::Real Abar_UU_01 = dDw_UL_01 + dDw_UL_10;
        amrex::Real Abar_UU_02 = dDw_UL_02 + dDw_UL_20;
        amrex::Real Abar_UU_11 = 2*dDw_UL_11;
        amrex::Real Abar_UU_12 = dDw_UL_12 + dDw_UL_21;
        amrex::Real Abar_UU_22 = 2*dDw_UL_22;
        
        amrex::Real psi = 1 + initial_data(i, j, k, InitIdx::u);
        
        state_fab(i, j, k, Idx::chi) = std::pow(psi,-4);

        state_fab(i, j, k, Idx::gamtilde_LL_00) = 1;
        state_fab(i, j, k, Idx::gamtilde_LL_01) = 0;
        state_fab(i, j, k, Idx::gamtilde_LL_02) = 0;
        state_fab(i, j, k, Idx::gamtilde_LL_11) = 1;
        state_fab(i, j, k, Idx::gamtilde_LL_12) = 0;
        state_fab(i, j, k, Idx::gamtilde_LL_22) = 1;

        state_fab(i, j, k, Idx::Khat) = 0;

        state_fab(i, j, k, Idx::Atilde_LL_00) = std::pow(psi,-6)*Abar_UU_00;
        state_fab(i, j, k, Idx::Atilde_LL_01) = std::pow(psi,-6)*Abar_UU_01;
        state_fab(i, j, k, Idx::Atilde_LL_02) = std::pow(psi,-6)*Abar_UU_02;
        state_fab(i, j, k, Idx::Atilde_LL_11) = std::pow(psi,-6)*Abar_UU_11;
        state_fab(i, j, k, Idx::Atilde_LL_12) = std::pow(psi,-6)*Abar_UU_12;
        state_fab(i, j, k, Idx::Atilde_LL_22) = std::pow(psi,-6)*Abar_UU_22;

        state_fab(i, j, k, Idx::theta) = 0;

        state_fab(i, j, k, Idx::Gamtilde_U_0) = 0;
        state_fab(i, j, k, Idx::Gamtilde_U_1) = 0;
        state_fab(i, j, k, Idx::Gamtilde_U_2) = 0;

        state_fab(i, j, k, Idx::alpha) = std::pow(psi,-2);

        state_fab(i, j, k, Idx::beta_U_0) = 0;
        state_fab(i, j, k, Idx::beta_U_1) = 0;
        state_fab(i, j, k, Idx::beta_U_2) = 0;
        
        
        state_fab(i, j, k, Idx::E_U_0) = 0;
        state_fab(i, j, k, Idx::E_U_1) = 0;
        state_fab(i, j, k, Idx::E_U_2) = 0;

        state_fab(i, j, k, Idx::A) = 0;
        


        state_fab(i, j, k, Idx::A_L_0) = initial_data(i, j, k, InitIdx::A_L_0);
        state_fab(i, j, k, Idx::A_L_1) = initial_data(i, j, k, InitIdx::A_L_1);
        state_fab(i, j, k, Idx::A_L_2) = initial_data(i, j, k, InitIdx::A_L_2);

        state_fab(i, j, k, Idx::Z) = 0;

        state_fab(i, j, k, Idx::phi1) = initial_data(i, j, k, InitIdx::phi1);

        state_fab(i, j, k, Idx::phi2) = initial_data(i, j, k, InitIdx::phi2);

        state_fab(i, j, k, Idx::Pi1) = 0;

        state_fab(i, j, k, Idx::Pi2) = 0;
        
        state_fab(i, j, k, Idx::sigma1) = initial_data(i, j, k, InitIdx::sigma1);

        state_fab(i, j, k, Idx::sigma2) = initial_data(i, j, k, InitIdx::sigma2);

        state_fab(i, j, k, Idx::delta1) = psi*psi*initial_data(i, j, k, InitIdx::dDtsigma1);

        state_fab(i, j, k, Idx::delta2) = psi*psi*initial_data(i, j, k, InitIdx::dDtsigma2);
        
        
}
#endif