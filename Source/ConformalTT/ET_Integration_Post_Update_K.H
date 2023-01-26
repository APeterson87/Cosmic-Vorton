#ifndef ET_INTEGRATION_POST_UPDATE_K_H
#define ET_INTEGRATION_POST_UPDATE_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_post_update(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x_U_0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x_U_1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x_U_2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2];
        
        amrex::Real R = 200.0;
        
        amrex::Real u = state_fab(i, j, k, Idx::u);
        state_fab(i, j, k, Idx::chi) = std::pow(1+u,-4);
        
        amrex::Real psi = 1 + u;
        
        amrex::Real dDtsigma1 = state_fab(i, j, k, Idx::dDtsigma1);
        amrex::Real dDtsigma2 = state_fab(i, j, k, Idx::dDtsigma2);
        
        
        state_fab(i, j, k, Idx::S_U_phi) = R/2*(dDtsigma1*dDtsigma1 + dDtsigma2*dDtsigma2)/(psi*psi*(x_U_0*x_U_0+x_U_1*x_U_1));
        
        
}
#endif