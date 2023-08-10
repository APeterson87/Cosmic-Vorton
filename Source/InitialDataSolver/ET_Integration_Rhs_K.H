#ifndef ET_INTEGRATION_RHS_K_H
#define ET_INTEGRATION_RHS_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_rhs(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& rhs_fab, 
        amrex::Array4<amrex::Real const> const& state_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx, 
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        
        amrex::Real u = state_fab(i, j, k, Idx::u);
        
        amrex::Real alphapsi = state_fab(i, j, k, Idx::alphapsi);
        
        amrex::Real V = state_fab(i, j, k, Idx::V);
        amrex::Real A = state_fab(i, j, k, Idx::A);
        amrex::Real B = state_fab(i, j, k, Idx::B);

        amrex::Real dDu_L_0 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::u) - 1.0/12.0*state_fab(i + 2, j, k, Idx::u) - 2.0/3.0*state_fab(i - 1, j, k, Idx::u) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::u))/dx[0];
        amrex::Real dDu_L_1 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::u) - 1.0/12.0*state_fab(i, j + 2, k, Idx::u) - 2.0/3.0*state_fab(i, j - 1, k, Idx::u) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::u))/dx[1];
        amrex::Real dDu_L_2 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::u) - 1.0/12.0*state_fab(i, j, k + 2, Idx::u) - 2.0/3.0*state_fab(i, j, k - 1, Idx::u) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::u))/dx[2];

        amrex::Real dDDu_LL_00 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::u) - 1.0/12.0*state_fab(i + 2, j, k, Idx::u) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::u) - 1.0/12.0*state_fab(i - 2, j, k, Idx::u) - 5.0/2.0*state_fab(i, j, k, Idx::u))/std::pow(dx[0], 2);
        amrex::Real dDDu_LL_01 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::u) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::u) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::u) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::u) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::u) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::u) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::u) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::u) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::u) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::u) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::u) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::u) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::u) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::u) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::u) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::u))/(dx[0]*dx[1]);
        amrex::Real dDDu_LL_02 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::u) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::u) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::u) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::u) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::u) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::u) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::u) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::u) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::u) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::u) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::u) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::u) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::u) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::u) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::u) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::u))/(dx[0]*dx[2]);
        amrex::Real dDDu_LL_11 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::u) - 1.0/12.0*state_fab(i, j + 2, k, Idx::u) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::u) - 1.0/12.0*state_fab(i, j - 2, k, Idx::u) - 5.0/2.0*state_fab(i, j, k, Idx::u))/std::pow(dx[1], 2);
        amrex::Real dDDu_LL_12 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::u) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::u) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::u) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::u) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::u) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::u) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::u) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::u) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::u) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::u) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::u) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::u) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::u) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::u) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::u) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::u))/(dx[1]*dx[2]);
        amrex::Real dDDu_LL_22 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::u) - 1.0/12.0*state_fab(i, j, k + 2, Idx::u) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::u) - 1.0/12.0*state_fab(i, j, k - 2, Idx::u) - 5.0/2.0*state_fab(i, j, k, Idx::u))/std::pow(dx[2], 2);
        
        
        amrex::Real dDDalphapsi_LL_00 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::alphapsi) - 1.0/12.0*state_fab(i + 2, j, k, Idx::alphapsi) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::alphapsi) - 1.0/12.0*state_fab(i - 2, j, k, Idx::alphapsi) - 5.0/2.0*state_fab(i, j, k, Idx::alphapsi))/std::pow(dx[0], 2);
        amrex::Real dDDalphapsi_LL_11 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::alphapsi) - 1.0/12.0*state_fab(i, j + 2, k, Idx::alphapsi) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::alphapsi) - 1.0/12.0*state_fab(i, j - 2, k, Idx::alphapsi) - 5.0/2.0*state_fab(i, j, k, Idx::alphapsi))/std::pow(dx[1], 2);
        amrex::Real dDDalphapsi_LL_22 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::alphapsi) - 1.0/12.0*state_fab(i, j, k + 2, Idx::alphapsi) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::alphapsi) - 1.0/12.0*state_fab(i, j, k - 2, Idx::alphapsi) - 5.0/2.0*state_fab(i, j, k, Idx::alphapsi))/std::pow(dx[2], 2);
        
        
        amrex::Real psi = 1+u;
        
        amrex::Real rho = V*std::pow(psi,5) + A*psi + B*std::pow(psi,-3);
        
        

        rhs_fab(i, j, k, Idx::u) = dDDu_LL_00 + dDDu_LL_11 + dDDu_LL_22+2*M_PI*rho;
        
        rhs_fab(i, j, k, Idx::alphapsi) = dDDalphapsi_LL_00 + dDDalphapsi_LL_11 + dDDalphapsi_LL_22 - 2*M_PI*rho*alphapsi/psi;
        

}
#endif