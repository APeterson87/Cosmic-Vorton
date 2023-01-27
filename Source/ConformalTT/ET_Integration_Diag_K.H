#ifndef ET_INTEGRATION_DIAG_K_H
#define ET_INTEGRATION_DIAG_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_diagnostics(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& diag, 
        amrex::Array4<amrex::Real const> const& state_fab,
        const amrex::Real time_lev,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx, 
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        
        amrex::Real w_U_0 = state_fab(i, j, k, Idx::w_U_0);
        amrex::Real w_U_1 = state_fab(i, j, k, Idx::w_U_1);
        amrex::Real u = state_fab(i, j, k, Idx::u);
        //diag(i, j, k, Diag::w_U_phi) = 1.0/(x0*x0 + x1*x1)*(-x1*w_U_0 + x0*w_U_1);
        
        amrex::Real dDw_UL_00 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::w_U_0) - 1.0/12.0*state_fab(i + 2, j, k, Idx::w_U_0) - 2.0/3.0*state_fab(i - 1, j, k, Idx::w_U_0) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::w_U_0))/dx[0];
        amrex::Real dDw_UL_01 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::w_U_0) - 1.0/12.0*state_fab(i, j + 2, k, Idx::w_U_0) - 2.0/3.0*state_fab(i, j - 1, k, Idx::w_U_0) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::w_U_0))/dx[1];
        amrex::Real dDw_UL_02 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::w_U_0) - 1.0/12.0*state_fab(i, j, k + 2, Idx::w_U_0) - 2.0/3.0*state_fab(i, j, k - 1, Idx::w_U_0) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::w_U_0))/dx[2];
        amrex::Real dDw_UL_10 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::w_U_1) - 1.0/12.0*state_fab(i + 2, j, k, Idx::w_U_1) - 2.0/3.0*state_fab(i - 1, j, k, Idx::w_U_1) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::w_U_1))/dx[0];
        amrex::Real dDw_UL_11 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::w_U_1) - 1.0/12.0*state_fab(i, j + 2, k, Idx::w_U_1) - 2.0/3.0*state_fab(i, j - 1, k, Idx::w_U_1) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::w_U_1))/dx[1];
        amrex::Real dDw_UL_12 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::w_U_1) - 1.0/12.0*state_fab(i, j, k + 2, Idx::w_U_1) - 2.0/3.0*state_fab(i, j, k - 1, Idx::w_U_1) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::w_U_1))/dx[2];
        amrex::Real dDw_UL_20 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::w_U_2) - 1.0/12.0*state_fab(i + 2, j, k, Idx::w_U_2) - 2.0/3.0*state_fab(i - 1, j, k, Idx::w_U_2) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::w_U_2))/dx[0];
        amrex::Real dDw_UL_21 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::w_U_2) - 1.0/12.0*state_fab(i, j + 2, k, Idx::w_U_2) - 2.0/3.0*state_fab(i, j - 1, k, Idx::w_U_2) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::w_U_2))/dx[1];
        amrex::Real dDw_UL_22 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::w_U_2) - 1.0/12.0*state_fab(i, j, k + 2, Idx::w_U_2) - 2.0/3.0*state_fab(i, j, k - 1, Idx::w_U_2) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::w_U_2))/dx[2];
        
        amrex::Real psi = 1+u;
        
        amrex::Real Atilde_UU_00 = 2*dDw_UL_00;
        amrex::Real Atilde_UU_01 = dDw_UL_01 + dDw_UL_10;
        amrex::Real Atilde_UU_02 = dDw_UL_02 + dDw_UL_20;
        amrex::Real Atilde_UU_11 = 2*dDw_UL_11;
        amrex::Real Atilde_UU_12 = dDw_UL_12 + dDw_UL_21;
        amrex::Real Atilde_UU_22 = 2*dDw_UL_22;
        
        amrex::Real AA = std::pow(Atilde_UU_00, 2) + 2*std::pow(Atilde_UU_01, 2) + 2*std::pow(Atilde_UU_02, 2) + std::pow(Atilde_UU_11, 2) + 2*std::pow(Atilde_UU_12, 2) + std::pow(Atilde_UU_22, 2);
        
        diag(i, j, k, Diag::w_U_phi) = (1.0/8.0)*AA/std::pow(psi,7);
        
}

AMREX_GPU_DEVICE 
inline 
amrex::Real
sum_C_constraint(const int i, const int j, const int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx, 
        const amrex::GeometryData& geom) noexcept 
{
    
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x_U_0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x_U_1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x_U_2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        
        return 0;
        
}
#endif