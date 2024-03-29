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
        
        amrex::Real r = std::sqrt(x_U_0*x_U_0 + x_U_1*x_U_1 + x_U_2*x_U_2);
        
        amrex::Real lambda = 2.0;
        amrex::Real eta = 0.05;
        amrex::Real e = 1.0;
        amrex::Real OMEGA = 0.0141421;
        
        amrex::Real chi = state_fab(i, j, k, Idx::chi);
        
        amrex::Real psi = std::pow(chi, 0.25);
        
        amrex::Real phi1 = state_fab(i, j, k, Idx::phi1);
        amrex::Real phi2 = state_fab(i, j, k, Idx::phi2);
        
        amrex::Real sigma1 = state_fab(i, j, k, Idx::sigma1);
        amrex::Real sigma2 = state_fab(i, j, k, Idx::sigma2);


        amrex::Real V = 0.25*lambda*std::pow(((phi1*phi1+phi2*phi2+sigma1*sigma1 + sigma2*sigma2)/2.0-eta*eta),2) + OMEGA*OMEGA/2*(sigma1*sigma1 + sigma2*sigma2);
        
        
        
        /*
        amrex::Real r = std::sqrt(x0*x0 + x1*x1 + x2*x2);
        
        amrex::Real alpha = std::pow(state_fab(i, j, k, Idx::psi),-2);
        
        //return (r < error_threshold);
        
        return (lev <= 3)*(r < error_threshold) || (lev > 3)*(std::fabs(alpha) < error_threshold);
        */
        
        return (V*1000000 > error_threshold);
        
        //return (lev <= 1)*(r < error_threshold) || (lev > 1)*(V*1000000 > error_threshold);
        
        //return (r < error_threshold);
        
        //return (lev <= 3)*(r < error_threshold) || (lev > 3)*(std::fabs(alpha) < error_threshold);
       
}
#endif