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

        amrex::Real v = 0.050000000000000003;

        amrex::Real Phi1 = state_fab(i, j, k, Idx::phi1);

        amrex::Real Phi2 = state_fab(i, j, k, Idx::phi2);
        
        amrex::Real sigma1 = state_fab(i, j, k, Idx::sigma1);

        amrex::Real sigma2 = state_fab(i, j, k, Idx::sigma2);
        
        amrex::Real V = (1.0/4.0)*lambda*std::pow((1.0/2.0)*std::pow(Phi1, 2) + (1.0/2.0)*std::pow(Phi2, 2) + (1.0/2.0)*std::pow(sigma1, 2) + (1.0/2.0)*std::pow(sigma2, 2) - std::pow(v, 2), 2);
        

        amrex::Real r = std::sqrt(x_U_0*x_U_0 + x_U_1*x_U_1 + x_U_2*x_U_2);
        
        return V*1000000.0 > error_threshold;
        
        //return (lev <= 1)*(r < error_threshold) || (lev > 1)*(std::fabs(V*1000000) > error_threshold);
}
#endif