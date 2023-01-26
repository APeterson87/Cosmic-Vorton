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
        
        amrex::Real x = x_U_0;
        amrex::Real y = x_U_1;
        amrex::Real z = x_U_2;
        
        amrex::Real R = 100;
        
        amrex::Real lambda = 2.0;
        amrex::Real eta = 0.05;
        amrex::Real e = 1.0;
        
        amrex::Real rdomain = std::sqrt(x*x + y*y + z*z);
        
        amrex::Real rsq = std::pow((std::sqrt(x*x + y*y)-R),2) + z*z;
        amrex::Real r = std::sqrt(rsq);
        
        amrex::Real x2 = r/2000;
        
        
        amrex::Real QF, QA;
        
        if(r < 600.0)
        {
            QF = 0.29365056608285506*x2 - 0.14019194815266911*std::cos(2*std::acos(2*x2 - 1)) + 0.12974328582629518*std::cos(3*std::acos(2*x2 - 1)) - 0.11630572639855538*std::cos(4*std::acos(2*x2 - 1)) + 0.10088393888199518*std::cos(5*std::acos(2*x2 - 1)) - 0.08454159948694942*std::cos(6*std::acos(2*x2 - 1)) + 0.068286336621182595*std::cos(7*std::acos(2*x2 - 1)) - 0.052975862368292458*std::cos(8*std::acos(2*x2 - 1)) + 0.039255812505915776*std::cos(9*std::acos(2*x2 - 1)) - 0.027533511470604073*std::cos(10*std::acos(2*x2 - 1)) + 0.017984725821405416*std::cos(11*std::acos(2*x2 - 1)) - 0.010586085218123278*std::cos(12*std::acos(2*x2 - 1)) + 0.0051630118112582606*std::cos(13*std::acos(2*x2 - 1)) - 0.0014431769640796503*std::cos(14*std::acos(2*x2 - 1)) - 0.00089301642676042556*std::cos(15*std::acos(2*x2 - 1)) + 0.0021705730828666461*std::cos(16*std::acos(2*x2 - 1)) - 0.0026895298256690292*std::cos(17*std::acos(2*x2 - 1)) + 0.0027063556106901171*std::cos(18*std::acos(2*x2 - 1)) - 0.0024257108651442143*std::cos(19*std::acos(2*x2 - 1)) + 0.0020001609135479846*std::cos(20*std::acos(2*x2 - 1)) - 0.0015352034249258642*std::cos(21*std::acos(2*x2 - 1)) + 0.0010971870591888061*std::cos(22*std::acos(2*x2 - 1)) - 0.00072220826543304556*std::cos(23*std::acos(2*x2 - 1)) + 0.00042465340112844349*std::cos(24*std::acos(2*x2 - 1)) - 0.00020460115507671059*std::cos(25*std::acos(2*x2 - 1)) + 5.3735362414718962e-5*std::cos(26*std::acos(2*x2 - 1)) + 4.027221374220194e-5*std::cos(27*std::acos(2*x2 - 1)) - 9.0760781521259659e-5*std::cos(28*std::acos(2*x2 - 1)) + 0.00011021855916877599*std::cos(29*std::acos(2*x2 - 1)) - 0.00010929083808659118*std::cos(30*std::acos(2*x2 - 1)) + 9.6404833979421567e-5*std::cos(31*std::acos(2*x2 - 1)) - 7.7795137123415348e-5*std::cos(32*std::acos(2*x2 - 1)) + 5.7761171693156144e-5*std::cos(33*std::acos(2*x2 - 1)) - 3.9029788196205673e-5*std::cos(34*std::acos(2*x2 - 1)) + 2.3143996980449361e-5*std::cos(35*std::acos(2*x2 - 1)) - 1.0812943186021273e-5*std::cos(36*std::acos(2*x2 - 1)) + 2.2258459937800346e-6*std::cos(37*std::acos(2*x2 - 1)) + 2.7496727188214276e-6*std::cos(38*std::acos(2*x2 - 1)) - 2.1511680285413171e-6*std::cos(39*std::acos(2*x2 - 1)) + 0.77862490140340335;
        
        QA = 0.40615100484695388*x2 - 0.188895758406779*std::cos(2*std::acos(2*x2 - 1)) + 0.16703515572607788*std::cos(3*std::acos(2*x2 - 1)) - 0.13983746098933553*std::cos(4*std::acos(2*x2 - 1)) + 0.11002574513047381*std::cos(5*std::acos(2*x2 - 1)) - 0.080294758448274661*std::cos(6*std::acos(2*x2 - 1)) + 0.052951586068632372*std::cos(7*std::acos(2*x2 - 1)) - 0.029668520044072895*std::cos(8*std::acos(2*x2 - 1)) + 0.011374276197784451*std::cos(9*std::acos(2*x2 - 1)) + 0.0017187136171678066*std::cos(10*std::acos(2*x2 - 1)) - 0.009984272440596971*std::cos(11*std::acos(2*x2 - 1)) + 0.014193627996733443*std::cos(12*std::acos(2*x2 - 1)) - 0.015315436707385853*std::cos(13*std::acos(2*x2 - 1)) + 0.014340649423501599*std::cos(14*std::acos(2*x2 - 1)) - 0.012154981577788097*std::cos(15*std::acos(2*x2 - 1)) + 0.0094652128636321699*std::cos(16*std::acos(2*x2 - 1)) - 0.0067744218452604546*std::cos(17*std::acos(2*x2 - 1)) + 0.0043935746533618029*std::cos(18*std::acos(2*x2 - 1)) - 2.0293988959267997e-5*std::cos(19*std::acos(2*x2 - 1)) + 0.0010554086052825574*std::cos(20*std::acos(2*x2 - 1)) - 9.6260900268953925e-5*std::cos(21*std::acos(2*x2 - 1)) - 0.0004800472883749185*std::cos(22*std::acos(2*x2 - 1)) + 0.00076484745660935702*std::cos(23*std::acos(2*x2 - 1)) - 0.00084611674425869032*std::cos(24*std::acos(2*x2 - 1)) + 0.00079909020838346556*std::cos(25*std::acos(2*x2 - 1)) - 0.00068249656647022608*std::cos(26*std::acos(2*x2 - 1)) + 0.00053848584244806097*std::cos(27*std::acos(2*x2 - 1)) - 0.00039464636960387031*std::cos(28*std::acos(2*x2 - 1)) + 0.00026693256047260333*std::cos(29*std::acos(2*x2 - 1)) - 0.00016273310341603439*std::cos(30*std::acos(2*x2 - 1)) + 8.3638533688583372e-5*std::cos(31*std::acos(2*x2 - 1)) - 2.7705041975540495e-5*std::cos(32*std::acos(2*x2 - 1)) - 8.834065831849619e-6*std::cos(33*std::acos(2*x2 - 1)) + 3.03703442972314e-5*std::cos(34*std::acos(2*x2 - 1)) - 4.1153237984089402e-5*std::cos(35*std::acos(2*x2 - 1)) + 4.4868728346587894e-5*std::cos(36*std::acos(2*x2 - 1)) - 4.4482454308093083e-5*std::cos(37*std::acos(2*x2 - 1)) + 4.2343989252447561e-5*std::cos(38*std::acos(2*x2 - 1)) - 4.0587977918535995e-5*std::cos(39*std::acos(2*x2 - 1)) + 0.69292997035750858;
        }
        
        else
        {
            QF = 1.0;
            QA = 1.0;
        }
        
        //QA = 0;
        /*Black Hole binary puncture data from initial data solver */
        
        
        

        amrex::Real psi = 1 + initial_data(i, j, k, InitIdx::u);
        
        state_fab(i, j, k, Idx::chi) = std::pow(psi,-4);

        state_fab(i, j, k, Idx::gamtilde_LL_00) = 1;
        state_fab(i, j, k, Idx::gamtilde_LL_01) = 0;
        state_fab(i, j, k, Idx::gamtilde_LL_02) = 0;
        state_fab(i, j, k, Idx::gamtilde_LL_11) = 1;
        state_fab(i, j, k, Idx::gamtilde_LL_12) = 0;
        state_fab(i, j, k, Idx::gamtilde_LL_22) = 1;

        state_fab(i, j, k, Idx::Khat) = 0;

        state_fab(i, j, k, Idx::Atilde_LL_00) = 0;
        state_fab(i, j, k, Idx::Atilde_LL_01) = 0;
        state_fab(i, j, k, Idx::Atilde_LL_02) = 0;
        state_fab(i, j, k, Idx::Atilde_LL_11) = 0;
        state_fab(i, j, k, Idx::Atilde_LL_12) = 0;
        state_fab(i, j, k, Idx::Atilde_LL_22) = 0;

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
        
        amrex::Real A_L_theta = QA/e;
        
        amrex::Real dthetadx = -x*z/(std::sqrt(rdomain*rdomain - z*z)*(R*R + rdomain*rdomain - 2*R*std::sqrt(rdomain*rdomain - z*z)));
        amrex::Real dthetady = -y*z/(std::sqrt(rdomain*rdomain - z*z)*(R*R + rdomain*rdomain - 2*R*std::sqrt(rdomain*rdomain - z*z)));
        
        amrex::Real f1 = -R + std::sqrt(rdomain*rdomain - z*z);
        amrex::Real f2 = 1 + z*z/((R - std::sqrt(x*x+y*y))*(R - std::sqrt(x*x+y*y)));
        
        amrex::Real dthetadz = 1.0/(f1*f2);


        state_fab(i, j, k, Idx::A_L_0) = dthetadx*A_L_theta;
        state_fab(i, j, k, Idx::A_L_1) = dthetady*A_L_theta;
        state_fab(i, j, k, Idx::A_L_2) = dthetadz*A_L_theta;

        state_fab(i, j, k, Idx::Z) = 0;
        
        amrex::Real costheta = (rdomain*rdomain-R*R-r*r)/(2*R*r);//(std::sqrt(x*x+y*y) - R)/r;
        amrex::Real sintheta = z/r;

        state_fab(i, j, k, Idx::phi1) = std::sqrt(2)*eta*QF*costheta;

        state_fab(i, j, k, Idx::phi2) = std::sqrt(2)*eta*QF*sintheta;

        state_fab(i, j, k, Idx::Pi1) = 0;

        state_fab(i, j, k, Idx::Pi2) = 0;
        
        
}
#endif