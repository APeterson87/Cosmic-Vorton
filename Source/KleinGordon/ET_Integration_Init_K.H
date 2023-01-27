#ifndef ET_INTEGRATION_INIT_K_H
#define ET_INTEGRATION_INIT_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_init(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab, 
        amrex::Real time, const amrex::GeometryData& geom) noexcept 
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
            QF = 0.20896353130802303*x2 - 0.10211223689432683*std::cos(2*std::acos(2*x2 - 1)) + 0.098272275860871644*std::cos(3*std::acos(2*x2 - 1)) - 0.09311669593970287*std::cos(4*std::acos(2*x2 - 1)) + 0.086849616077820074*std::cos(5*std::acos(2*x2 - 1)) - 0.079708769475895697*std::cos(6*std::acos(2*x2 - 1)) + 0.071953768006534358*std::cos(7*std::acos(2*x2 - 1)) - 0.063848697197569859*std::cos(8*std::acos(2*x2 - 1)) + 0.055650299194329469*std::cos(9*std::acos(2*x2 - 1)) - 0.047593464787780641*std::cos(10*std::acos(2*x2 - 1)) + 0.039883225577238141*std::cos(11*std::acos(2*x2 - 1)) - 0.032686183643460286*std::cos(12*std::acos(2*x2 - 1)) + 0.02612819254295308*std::cos(13*std::acos(2*x2 - 1)) - 0.020292382052463247*std::cos(14*std::acos(2*x2 - 1)) + 0.01522213038425645*std::cos(15*std::acos(2*x2 - 1)) - 0.010924252356053772*std::cos(16*std::acos(2*x2 - 1)) + 0.0073752856296346869*std::cos(17*std::acos(2*x2 - 1)) - 0.0045273456337332485*std::cos(18*std::acos(2*x2 - 1)) + 0.002315289371541962*std::cos(19*std::acos(2*x2 - 1)) - 0.00066280550845689359*std::cos(20*std::acos(2*x2 - 1)) - 0.00051148664980081929*std::cos(21*std::acos(2*x2 - 1)) + 0.00128933346799983*std::cos(22*std::acos(2*x2 - 1)) - 0.0017488069689017283*std::cos(23*std::acos(2*x2 - 1)) + 0.0019616005999074405*std::cos(24*std::acos(2*x2 - 1)) - 0.0019911392109407018*std::cos(25*std::acos(2*x2 - 1)) + 0.0018916966509791022*std::cos(26*std::acos(2*x2 - 1)) - 0.0017081781765125453*std::cos(27*std::acos(2*x2 - 1)) + 0.0014764741605359892*std::cos(28*std::acos(2*x2 - 1)) - 0.0012242522803171327*std::cos(29*std::acos(2*x2 - 1)) + 0.0009718877236291959*std::cos(30*std::acos(2*x2 - 1)) - 0.00073370056042746171*std::cos(31*std::acos(2*x2 - 1)) + 0.00051890923985509901*std::cos(32*std::acos(2*x2 - 1)) - 0.00033300174679158981*std::cos(33*std::acos(2*x2 - 1)) + 0.00017827281487187499*std::cos(34*std::acos(2*x2 - 1)) - 5.5370134037530105e-5*std::cos(35*std::acos(2*x2 - 1)) - 3.7086889872764428e-5*std::cos(36*std::acos(2*x2 - 1)) + 0.00010017708061388353*std::cos(37*std::acos(2*x2 - 1)) - 0.00013703265962246997*std::cos(38*std::acos(2*x2 - 1)) + 7.3910347924175529e-5*std::cos(39*std::acos(2*x2 - 1)) + 0.84287701272714854;
        
            QA = 0.28988630310372554*x2 - 0.1398642903953197*std::cos(2*std::acos(2*x2 - 1)) + 0.13172005389585492*std::cos(3*std::acos(2*x2 - 1)) - 0.12095679058319467*std::cos(4*std::acos(2*x2 - 1)) + 0.10815285097346047*std::cos(5*std::acos(2*x2 - 1)) - 0.093960853691770124*std::cos(6*std::acos(2*x2 - 1)) + 0.079068052936334232*std::cos(7*std::acos(2*x2 - 1)) - 0.064136937213247092*std::cos(8*std::acos(2*x2 - 1)) + 0.049770808044987*std::cos(9*std::acos(2*x2 - 1)) - 0.036471066334833831*std::cos(10*std::acos(2*x2 - 1)) + 0.02462200101099675*std::cos(11*std::acos(2*x2 - 1)) - 0.014473879206189026*std::cos(12*std::acos(2*x2 - 1)) + 0.0061502141171301816*std::cos(13*std::acos(2*x2 - 1)) + 0.00034563655546149263*std::cos(14*std::acos(2*x2 - 1)) - 0.0051072823183213947*std::cos(15*std::acos(2*x2 - 1)) + 0.0083052499577140719*std::cos(16*std::acos(2*x2 - 1)) - 0.01015893777361454*std::cos(17*std::acos(2*x2 - 1)) + 0.010914690423277048*std::cos(18*std::acos(2*x2 - 1)) + 0.00031361786200028406*std::cos(19*std::acos(2*x2 - 1)) + 0.010120549698957444*std::cos(20*std::acos(2*x2 - 1)) - 0.0090210520737649749*std::cos(21*std::acos(2*x2 - 1)) + 0.0077039949719643272*std::cos(22*std::acos(2*x2 - 1)) - 0.0063128043827225047*std::cos(23*std::acos(2*x2 - 1)) + 0.0049550021659360774*std::cos(24*std::acos(2*x2 - 1)) - 0.003704621451918591*std::cos(25*std::acos(2*x2 - 1)) + 0.0026068686814357984*std::cos(26*std::acos(2*x2 - 1)) - 0.0016833342755716886*std::cos(27*std::acos(2*x2 - 1)) + 0.00093753344869824257*std::cos(28*std::acos(2*x2 - 1)) - 0.00036019513798869751*std::cos(29*std::acos(2*x2 - 1)) - 6.6336624994230153e-5*std::cos(30*std::acos(2*x2 - 1)) + 0.00036400158431543206*std::cos(31*std::acos(2*x2 - 1)) - 0.00055662793152901416*std::cos(32*std::acos(2*x2 - 1)) + 0.00066702157724501596*std::cos(33*std::acos(2*x2 - 1)) - 0.00071690663377086544*std::cos(34*std::acos(2*x2 - 1)) + 0.00072422878233856583*std::cos(35*std::acos(2*x2 - 1)) - 0.00070622835123692815*std::cos(36*std::acos(2*x2 - 1)) + 0.00067478185551544945*std::cos(37*std::acos(2*x2 - 1)) - 0.00064572911739346414*std::cos(38*std::acos(2*x2 - 1)) + 0.00062723572400056812*std::cos(39*std::acos(2*x2 - 1)) + 0.78172296862817137;
        }
        
        else
        {
            QF = 1.0;
            QA = 1.0;
        }
        
        
        
        
        
        amrex::Real costheta = (rdomain*rdomain-R*R-r*r)/(2*R*r);//(std::sqrt(x*x+y*y) - R)/r;
        amrex::Real sintheta = z/r;

        state_fab(i, j, k, Idx::phi1) = std::sqrt(2)*eta*QF*costheta;

        state_fab(i, j, k, Idx::phi2) = std::sqrt(2)*eta*QF*sintheta;
}
#endif