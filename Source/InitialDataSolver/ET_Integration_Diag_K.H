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
        
        
        amrex::Real u = state_fab(i, j, k, Idx::u);
        
        amrex::Real psi = 1+u;
        amrex::Real V = state_fab(i, j, k, Idx::V);
        amrex::Real A = state_fab(i, j, k, Idx::A);
        amrex::Real B = state_fab(i, j, k, Idx::B);
        
        amrex::Real rho = V + A/std::pow(psi,4) + B*std::pow(psi,-8);
        
        diag(i, j, k, Diag::rhoADM) = rho;
        
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