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
        
        amrex::Real x_U_0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x_U_1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x_U_2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        
        amrex::Real x = x_U_0;
        amrex::Real y = x_U_1;
        amrex::Real z = x_U_2;
        
        amrex::Real u = state_fab(i, j, k, Idx::u);
        
        amrex::Real V = state_fab(i, j, k, Idx::V);
        amrex::Real A = state_fab(i, j, k, Idx::A);
        amrex::Real B = state_fab(i, j, k, Idx::B);
        amrex::Real C = state_fab(i, j, k, Idx::C);

        amrex::Real dDDu_LL_00 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::u) - 1.0/12.0*state_fab(i + 2, j, k, Idx::u) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::u) - 1.0/12.0*state_fab(i - 2, j, k, Idx::u) - 5.0/2.0*state_fab(i, j, k, Idx::u))/std::pow(dx[0], 2);
        amrex::Real dDDu_LL_01 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::u) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::u) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::u) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::u) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::u) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::u) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::u) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::u) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::u) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::u) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::u) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::u) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::u) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::u) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::u) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::u))/(dx[0]*dx[1]);
        amrex::Real dDDu_LL_02 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::u) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::u) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::u) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::u) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::u) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::u) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::u) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::u) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::u) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::u) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::u) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::u) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::u) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::u) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::u) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::u))/(dx[0]*dx[2]);
        amrex::Real dDDu_LL_11 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::u) - 1.0/12.0*state_fab(i, j + 2, k, Idx::u) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::u) - 1.0/12.0*state_fab(i, j - 2, k, Idx::u) - 5.0/2.0*state_fab(i, j, k, Idx::u))/std::pow(dx[1], 2);
        amrex::Real dDDu_LL_12 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::u) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::u) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::u) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::u) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::u) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::u) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::u) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::u) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::u) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::u) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::u) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::u) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::u) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::u) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::u) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::u))/(dx[1]*dx[2]);
        amrex::Real dDDu_LL_22 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::u) - 1.0/12.0*state_fab(i, j, k + 2, Idx::u) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::u) - 1.0/12.0*state_fab(i, j, k - 2, Idx::u) - 5.0/2.0*state_fab(i, j, k, Idx::u))/std::pow(dx[2], 2);
        
        
        amrex::Real dDw_UL_00 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::w_U_0) - 1.0/12.0*state_fab(i + 2, j, k, Idx::w_U_0) - 2.0/3.0*state_fab(i - 1, j, k, Idx::w_U_0) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::w_U_0))/dx[0];
        amrex::Real dDw_UL_01 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::w_U_0) - 1.0/12.0*state_fab(i, j + 2, k, Idx::w_U_0) - 2.0/3.0*state_fab(i, j - 1, k, Idx::w_U_0) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::w_U_0))/dx[1];
        amrex::Real dDw_UL_02 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::w_U_0) - 1.0/12.0*state_fab(i, j, k + 2, Idx::w_U_0) - 2.0/3.0*state_fab(i, j, k - 1, Idx::w_U_0) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::w_U_0))/dx[2];
        amrex::Real dDw_UL_10 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::w_U_1) - 1.0/12.0*state_fab(i + 2, j, k, Idx::w_U_1) - 2.0/3.0*state_fab(i - 1, j, k, Idx::w_U_1) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::w_U_1))/dx[0];
        amrex::Real dDw_UL_11 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::w_U_1) - 1.0/12.0*state_fab(i, j + 2, k, Idx::w_U_1) - 2.0/3.0*state_fab(i, j - 1, k, Idx::w_U_1) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::w_U_1))/dx[1];
        amrex::Real dDw_UL_12 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::w_U_1) - 1.0/12.0*state_fab(i, j, k + 2, Idx::w_U_1) - 2.0/3.0*state_fab(i, j, k - 1, Idx::w_U_1) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::w_U_1))/dx[2];
        amrex::Real dDw_UL_20 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::w_U_2) - 1.0/12.0*state_fab(i + 2, j, k, Idx::w_U_2) - 2.0/3.0*state_fab(i - 1, j, k, Idx::w_U_2) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::w_U_2))/dx[0];
        amrex::Real dDw_UL_21 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::w_U_2) - 1.0/12.0*state_fab(i, j + 2, k, Idx::w_U_2) - 2.0/3.0*state_fab(i, j - 1, k, Idx::w_U_2) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::w_U_2))/dx[1];
        amrex::Real dDw_UL_22 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::w_U_2) - 1.0/12.0*state_fab(i, j, k + 2, Idx::w_U_2) - 2.0/3.0*state_fab(i, j, k - 1, Idx::w_U_2) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::w_U_2))/dx[2];

        amrex::Real dDDw_ULL_000 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::w_U_0) - 1.0/12.0*state_fab(i + 2, j, k, Idx::w_U_0) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::w_U_0) - 1.0/12.0*state_fab(i - 2, j, k, Idx::w_U_0) - 5.0/2.0*state_fab(i, j, k, Idx::w_U_0))/std::pow(dx[0], 2);
        amrex::Real dDDw_ULL_001 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::w_U_0) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::w_U_0) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::w_U_0) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::w_U_0) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::w_U_0) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::w_U_0) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::w_U_0) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::w_U_0) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::w_U_0) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::w_U_0) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::w_U_0) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::w_U_0) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::w_U_0) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::w_U_0) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::w_U_0) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::w_U_0))/(dx[0]*dx[1]);
        amrex::Real dDDw_ULL_002 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::w_U_0) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::w_U_0) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::w_U_0) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::w_U_0) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::w_U_0) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::w_U_0) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::w_U_0) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::w_U_0) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::w_U_0) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::w_U_0) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::w_U_0) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::w_U_0) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::w_U_0) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::w_U_0) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::w_U_0) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::w_U_0))/(dx[0]*dx[2]);
        amrex::Real dDDw_ULL_011 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::w_U_0) - 1.0/12.0*state_fab(i, j + 2, k, Idx::w_U_0) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::w_U_0) - 1.0/12.0*state_fab(i, j - 2, k, Idx::w_U_0) - 5.0/2.0*state_fab(i, j, k, Idx::w_U_0))/std::pow(dx[1], 2);
        amrex::Real dDDw_ULL_012 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::w_U_0) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::w_U_0) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::w_U_0) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::w_U_0) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::w_U_0) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::w_U_0) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::w_U_0) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::w_U_0) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::w_U_0) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::w_U_0) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::w_U_0) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::w_U_0) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::w_U_0) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::w_U_0) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::w_U_0) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::w_U_0))/(dx[1]*dx[2]);
        amrex::Real dDDw_ULL_022 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::w_U_0) - 1.0/12.0*state_fab(i, j, k + 2, Idx::w_U_0) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::w_U_0) - 1.0/12.0*state_fab(i, j, k - 2, Idx::w_U_0) - 5.0/2.0*state_fab(i, j, k, Idx::w_U_0))/std::pow(dx[2], 2);
        amrex::Real dDDw_ULL_100 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::w_U_1) - 1.0/12.0*state_fab(i + 2, j, k, Idx::w_U_1) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::w_U_1) - 1.0/12.0*state_fab(i - 2, j, k, Idx::w_U_1) - 5.0/2.0*state_fab(i, j, k, Idx::w_U_1))/std::pow(dx[0], 2);
        amrex::Real dDDw_ULL_101 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::w_U_1) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::w_U_1) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::w_U_1) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::w_U_1) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::w_U_1) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::w_U_1) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::w_U_1) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::w_U_1) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::w_U_1) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::w_U_1) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::w_U_1) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::w_U_1) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::w_U_1) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::w_U_1) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::w_U_1) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::w_U_1))/(dx[0]*dx[1]);
        amrex::Real dDDw_ULL_102 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::w_U_1) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::w_U_1) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::w_U_1) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::w_U_1) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::w_U_1) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::w_U_1) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::w_U_1) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::w_U_1) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::w_U_1) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::w_U_1) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::w_U_1) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::w_U_1) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::w_U_1) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::w_U_1) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::w_U_1) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::w_U_1))/(dx[0]*dx[2]);
        amrex::Real dDDw_ULL_111 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::w_U_1) - 1.0/12.0*state_fab(i, j + 2, k, Idx::w_U_1) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::w_U_1) - 1.0/12.0*state_fab(i, j - 2, k, Idx::w_U_1) - 5.0/2.0*state_fab(i, j, k, Idx::w_U_1))/std::pow(dx[1], 2);
        amrex::Real dDDw_ULL_112 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::w_U_1) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::w_U_1) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::w_U_1) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::w_U_1) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::w_U_1) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::w_U_1) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::w_U_1) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::w_U_1) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::w_U_1) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::w_U_1) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::w_U_1) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::w_U_1) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::w_U_1) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::w_U_1) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::w_U_1) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::w_U_1))/(dx[1]*dx[2]);
        amrex::Real dDDw_ULL_122 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::w_U_1) - 1.0/12.0*state_fab(i, j, k + 2, Idx::w_U_1) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::w_U_1) - 1.0/12.0*state_fab(i, j, k - 2, Idx::w_U_1) - 5.0/2.0*state_fab(i, j, k, Idx::w_U_1))/std::pow(dx[2], 2);
        amrex::Real dDDw_ULL_200 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::w_U_2) - 1.0/12.0*state_fab(i + 2, j, k, Idx::w_U_2) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::w_U_2) - 1.0/12.0*state_fab(i - 2, j, k, Idx::w_U_2) - 5.0/2.0*state_fab(i, j, k, Idx::w_U_2))/std::pow(dx[0], 2);
        amrex::Real dDDw_ULL_201 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::w_U_2) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::w_U_2) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::w_U_2) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::w_U_2) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::w_U_2) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::w_U_2) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::w_U_2) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::w_U_2) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::w_U_2) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::w_U_2) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::w_U_2) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::w_U_2) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::w_U_2) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::w_U_2) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::w_U_2) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::w_U_2))/(dx[0]*dx[1]);
        amrex::Real dDDw_ULL_202 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::w_U_2) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::w_U_2) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::w_U_2) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::w_U_2) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::w_U_2) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::w_U_2) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::w_U_2) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::w_U_2) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::w_U_2) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::w_U_2) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::w_U_2) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::w_U_2) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::w_U_2) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::w_U_2) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::w_U_2) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::w_U_2))/(dx[0]*dx[2]);
        amrex::Real dDDw_ULL_211 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::w_U_2) - 1.0/12.0*state_fab(i, j + 2, k, Idx::w_U_2) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::w_U_2) - 1.0/12.0*state_fab(i, j - 2, k, Idx::w_U_2) - 5.0/2.0*state_fab(i, j, k, Idx::w_U_2))/std::pow(dx[1], 2);
        amrex::Real dDDw_ULL_212 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::w_U_2) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::w_U_2) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::w_U_2) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::w_U_2) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::w_U_2) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::w_U_2) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::w_U_2) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::w_U_2) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::w_U_2) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::w_U_2) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::w_U_2) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::w_U_2) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::w_U_2) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::w_U_2) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::w_U_2) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::w_U_2))/(dx[1]*dx[2]);
        amrex::Real dDDw_ULL_222 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::w_U_2) - 1.0/12.0*state_fab(i, j, k + 2, Idx::w_U_2) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::w_U_2) - 1.0/12.0*state_fab(i, j, k - 2, Idx::w_U_2) - 5.0/2.0*state_fab(i, j, k, Idx::w_U_2))/std::pow(dx[2], 2);

        
        amrex::Real psi = 1+u;
        
        amrex::Real rho_adj = V*std::pow(psi,5) + A*psi + B*psi + C*std::pow(psi,-3);
        
        amrex::Real Atilde_UU_00 = 2*dDw_UL_00;
        amrex::Real Atilde_UU_01 = dDw_UL_01 + dDw_UL_10;
        amrex::Real Atilde_UU_02 = dDw_UL_02 + dDw_UL_20;
        amrex::Real Atilde_UU_11 = 2*dDw_UL_11;
        amrex::Real Atilde_UU_12 = dDw_UL_12 + dDw_UL_21;
        amrex::Real Atilde_UU_22 = 2*dDw_UL_22;

        amrex::Real AA = std::pow(Atilde_UU_00, 2) + 2*std::pow(Atilde_UU_01, 2) + 2*std::pow(Atilde_UU_02, 2) + std::pow(Atilde_UU_11, 2) + 2*std::pow(Atilde_UU_12, 2) + std::pow(Atilde_UU_22, 2);
        
        rho_adj += (1.0/8.0)*AA/std::pow(psi,7);
        
        amrex::Real S_U_phi = state_fab(i, j, k, Idx::S_U_phi);
        
        amrex::Real S_U_0 = -y*S_U_phi;
        amrex::Real S_U_1 = x*S_U_phi;
        
        

        rhs_fab(i, j, k, Idx::u) = dDDu_LL_00 + dDDu_LL_11 + dDDu_LL_22+2*M_PI*rho_adj;
        
        rhs_fab(i, j, k, Idx::w_U_0) = dDDw_ULL_000 + dDDw_ULL_011 + dDDw_ULL_022-8*M_PI*std::pow(psi,10)*S_U_0;
        
        rhs_fab(i, j, k, Idx::w_U_1) = dDDw_ULL_100 + dDDw_ULL_111 + dDDw_ULL_122-8*M_PI*std::pow(psi,10)*S_U_1;
        

}
#endif