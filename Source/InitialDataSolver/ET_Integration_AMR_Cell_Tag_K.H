#ifndef ET_INTEGRATION_AMR_CELL_TAG_K_H
#define ET_INTEGRATION_AMR_CELL_TAG_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
bool
state_is_tagged(int i, int j, int k, 
        int lev,
        amrex::Array4<amrex::Real const> const& state_fab,
        amrex::Real error_threshold,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{

        const auto domain_xlo = geom.ProbLo(); 

        amrex::Real x_U_0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x_U_1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x_U_2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2];
        
        amrex::Real lambda = 2.0;
        amrex::Real eta = 0.05;
        amrex::Real e = std::sqrt(lambda/2.0);
        
        amrex::Real u = state_fab(i, j, k, Idx::u);
        amrex::Real V = state_fab(i, j, k, Idx::V);
        amrex::Real A = state_fab(i, j, k, Idx::A);
        amrex::Real B = state_fab(i, j, k, Idx::B);
        
        amrex::Real chi = state_fab(i, j, k, Idx::chi);
        
        /*amrex::Real v = 4*V/(lambda*eta*eta*eta*eta); */
        
        amrex::Real r = std::sqrt(x_U_0*x_U_0 + x_U_1*x_U_1 + x_U_2*x_U_2);
        
        //return (chi < error_threshold);
        return (lev <= 3)*(r < error_threshold) || (lev > 3)*(chi < error_threshold);
        //return (V*1000000 > error_threshold);
}
#endif