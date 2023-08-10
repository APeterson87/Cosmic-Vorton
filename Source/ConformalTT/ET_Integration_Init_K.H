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
        
        amrex::Real lambda = 2.0;
        amrex::Real eta = 0.05;
        amrex::Real e = 1.0;
        
        amrex::Real R = 200.0;
        amrex::Real omega = 2.0/R;
        amrex::Real OMEGA = 0.0141421;
        
        amrex::Real x = x_U_0;
        amrex::Real y = x_U_1;
        amrex::Real z = x_U_2;
        
        
        amrex::Real rdomain = std::sqrt(x*x + y*y + z*z);
        
        amrex::Real rsq = std::pow((std::sqrt(x*x + y*y)-R),2) + z*z;
        amrex::Real r = std::sqrt(rsq);
        
        amrex::Real x2 = r/2000;
        
        amrex::Real QF, QG, QA, QU, QW, QDF, QDG, QDA;
        
        //amrex::Real x2 = std::sqrt(x_U_0*x_U_0 + x_U_1*x_U_1 + x_U_2*x_U_2)/1000;
        
        if(rdomain < 2000.0)
        {
        
            QF = 0.22533462278456101*x2 - 0.10945612817914691*std::cos(2*std::acos(2*x2 - 1)) + 0.10431766367601869*std::cos(3*std::acos(2*x2 - 1)) - 0.097546655748109334*std::cos(4*std::acos(2*x2 - 1)) + 0.089510133550337254*std::cos(5*std::acos(2*x2 - 1)) - 0.080610682030641811*std::cos(6*std::acos(2*x2 - 1)) + 0.071251421427643916*std::cos(7*std::acos(2*x2 - 1)) - 0.061805500349823461*std::cos(8*std::acos(2*x2 - 1)) + 0.052594714193575406*std::cos(9*std::acos(2*x2 - 1)) - 0.043876656742942333*std::cos(10*std::acos(2*x2 - 1)) + 0.035841096098399958*std::cos(11*std::acos(2*x2 - 1)) - 0.028612465774098126*std::cos(12*std::acos(2*x2 - 1)) + 0.022257315041256571*std::cos(13*std::acos(2*x2 - 1)) - 0.016793596396653006*std::cos(14*std::acos(2*x2 - 1)) + 0.012200964243527259*std::cos(15*std::acos(2*x2 - 1)) - 0.0084302951006838336*std::cos(16*std::acos(2*x2 - 1)) + 0.0054124657032948616*std::cos(17*std::acos(2*x2 - 1)) - 0.0030656980640317542*std::cos(18*std::acos(2*x2 - 1)) + 0.0013018671382866893*std::cos(19*std::acos(2*x2 - 1)) - 3.1518884351589699e-5*std::cos(20*std::acos(2*x2 - 1)) - 0.00083206753436851008*std::cos(21*std::acos(2*x2 - 1)) + 0.0013699000503674119*std::cos(22*std::acos(2*x2 - 1)) - 0.0016551491559494945*std::cos(23*std::acos(2*x2 - 1)) + 0.0017519066839266559*std::cos(24*std::acos(2*x2 - 1)) - 0.001714629716017331*std::cos(25*std::acos(2*x2 - 1)) + 0.0015882049890543724*std::cos(26*std::acos(2*x2 - 1)) - 0.0014084469831101426*std::cos(27*std::acos(2*x2 - 1)) + 0.0012029268316908759*std::cos(28*std::acos(2*x2 - 1)) - 0.00099198307254462839*std::cos(29*std::acos(2*x2 - 1)) + 0.00078982592827502376*std::cos(30*std::acos(2*x2 - 1)) - 0.00060564030087112766*std::cos(31*std::acos(2*x2 - 1)) + 0.00044463438413823442*std::cos(32*std::acos(2*x2 - 1)) - 0.00030898808537635145*std::cos(33*std::acos(2*x2 - 1)) + 0.00019868051689338195*std::cos(34*std::acos(2*x2 - 1)) - 0.00011218541978350385*std::cos(35*std::acos(2*x2 - 1)) + 4.7034849422766117e-5*std::cos(36*std::acos(2*x2 - 1)) - 2.6103981979122559e-7*std::cos(37*std::acos(2*x2 - 1)) - 3.1274691030779012e-5*std::cos(38*std::acos(2*x2 - 1)) + 5.0643489226492172e-5*std::cos(39*std::acos(2*x2 - 1)) - 6.0696609002510588e-5*std::cos(40*std::acos(2*x2 - 1)) + 6.396318595381506e-5*std::cos(41*std::acos(2*x2 - 1)) - 6.2617324015270398e-5*std::cos(42*std::acos(2*x2 - 1)) + 5.8444225745590134e-5*std::cos(43*std::acos(2*x2 - 1)) - 5.28935475072774e-5*std::cos(44*std::acos(2*x2 - 1)) + 4.7039560228223642e-5*std::cos(45*std::acos(2*x2 - 1)) - 4.1745016398340301e-5*std::cos(46*std::acos(2*x2 - 1)) + 3.7485315785318329e-5*std::cos(47*std::acos(2*x2 - 1)) - 3.4806876540786208e-5*std::cos(48*std::acos(2*x2 - 1)) + 1.6823066280502799e-5*std::cos(49*std::acos(2*x2 - 1)) + 0.8304528057089281;


         QG = -0.21648378939715729*x2 + 0.098469202339645692*std::cos(2*std::acos(2*x2 - 1)) - 0.084070606851148103*std::cos(3*std::acos(2*x2 - 1)) + 0.067307725563056031*std::cos(4*std::acos(2*x2 - 1)) - 0.050436160128964957*std::cos(5*std::acos(2*x2 - 1)) + 0.03522220626410582*std::cos(6*std::acos(2*x2 - 1)) - 0.022696686037049257*std::cos(7*std::acos(2*x2 - 1)) + 0.013165864446329927*std::cos(8*std::acos(2*x2 - 1)) - 0.0064015780909980917*std::cos(9*std::acos(2*x2 - 1)) + 0.0018942768046816647*std::cos(10*std::acos(2*x2 - 1)) + 0.00092919467604489871*std::cos(11*std::acos(2*x2 - 1)) - 0.0025714706347877664*std::cos(12*std::acos(2*x2 - 1)) + 0.003414447505415846*std::cos(13*std::acos(2*x2 - 1)) - 0.0037252741482617226*std::cos(14*std::acos(2*x2 - 1)) + 0.0036859786240546464*std::cos(15*std::acos(2*x2 - 1)) - 0.0034236189562211537*std::cos(16*std::acos(2*x2 - 1)) + 0.0030309477883453968*std::cos(17*std::acos(2*x2 - 1)) - 0.002577227800809709*std::cos(18*std::acos(2*x2 - 1)) + 0.0021130017372718883*std::cos(19*std::acos(2*x2 - 1)) - 0.0016726685475709395*std::cos(20*std::acos(2*x2 - 1)) + 0.0012770205696199697*std::cos(21*std::acos(2*x2 - 1)) - 0.00093624471263307996*std::cos(22*std::acos(2*x2 - 1)) + 0.00065302428844877137*std::cos(23*std::acos(2*x2 - 1)) - 0.00042524308685525547*std::cos(24*std::acos(2*x2 - 1)) + 0.00024802502709461987*std::cos(25*std::acos(2*x2 - 1)) - 0.00011511318441405131*std::cos(26*std::acos(2*x2 - 1)) + 1.9743554614597318e-5*std::cos(27*std::acos(2*x2 - 1)) + 4.480796820542954e-5*std::cos(28*std::acos(2*x2 - 1)) - 8.4874637232357228e-5*std::cos(29*std::acos(2*x2 - 1)) + 0.0001061725884347806*std::cos(30*std::acos(2*x2 - 1)) - 0.00011366780042494894*std::cos(31*std::acos(2*x2 - 1)) + 0.00011152084601558788*std::cos(32*std::acos(2*x2 - 1)) - 0.00010309537808171745*std::cos(33*std::acos(2*x2 - 1)) + 9.1013534597605066e-5*std::cos(34*std::acos(2*x2 - 1)) - 7.7240553262412974e-5*std::cos(35*std::acos(2*x2 - 1)) + 6.3183077135100195e-5*std::cos(36*std::acos(2*x2 - 1)) - 4.9789618430494539e-5*std::cos(37*std::acos(2*x2 - 1)) + 3.7645756046392175e-5*std::cos(38*std::acos(2*x2 - 1)) - 2.7059901711171845e-5*std::cos(39*std::acos(2*x2 - 1)) + 1.813768302220993e-5*std::cos(40*std::acos(2*x2 - 1)) - 1.0844377163634259e-5*std::cos(41*std::acos(2*x2 - 1)) + 5.0556898746169085e-6*std::cos(42*std::acos(2*x2 - 1)) - 5.9771321187931419e-7*std::cos(43*std::acos(2*x2 - 1)) - 2.7227937402465638e-6*std::cos(44*std::acos(2*x2 - 1)) + 5.0964743125390325e-6*std::cos(45*std::acos(2*x2 - 1)) - 6.6966204249660554e-6*std::cos(46*std::acos(2*x2 - 1)) + 7.669655371039678e-6*std::cos(47*std::acos(2*x2 - 1)) - 8.1372843914242896e-6*std::cos(48*std::acos(2*x2 - 1)) + 4.0871401378874758e-6*std::cos(49*std::acos(2*x2 - 1)) + 0.16409535865306363;


         QA = 0.2917941657354306*x2 - 0.14072561186398394*std::cos(2*std::acos(2*x2 - 1)) + 0.13243370108074981*std::cos(3*std::acos(2*x2 - 1)) - 0.12148305771027593*std::cos(4*std::acos(2*x2 - 1)) + 0.10846433853509929*std::cos(5*std::acos(2*x2 - 1)) - 0.094049831230096095*std::cos(6*std::acos(2*x2 - 1)) + 0.07894088393987235*std::cos(7*std::acos(2*x2 - 1)) - 0.063817521184742376*std::cos(8*std::acos(2*x2 - 1)) + 0.049293058378657439*std::cos(9*std::acos(2*x2 - 1)) - 0.035879310973440215*std::cos(10*std::acos(2*x2 - 1)) + 0.02396308482514635*std::cos(11*std::acos(2*x2 - 1)) - 0.013796238973623354*std::cos(12*std::acos(2*x2 - 1)) + 0.0054973407027563738*std::cos(13*std::acos(2*x2 - 1)) + 0.00093580946002157564*std::cos(14*std::acos(2*x2 - 1)) - 0.0056064445791016784*std::cos(15*std::acos(2*x2 - 1)) + 0.0086941140529114037*std::cos(16*std::acos(2*x2 - 1)) - 0.010428693491543333*std::cos(17*std::acos(2*x2 - 1)) + 0.011065186895628621*std::cos(18*std::acos(2*x2 - 1)) - 0.010861815223881691*std::cos(19*std::acos(2*x2 - 1)) + 0.010062329958571247*std::cos(20*std::acos(2*x2 - 1)) - 0.0088834435784720042*std::cos(21*std::acos(2*x2 - 1)) + 0.0075071193956154827*std::cos(22*std::acos(2*x2 - 1)) - 0.0060773922311331725*std::cos(23*std::acos(2*x2 - 1)) + 0.0047007706915367133*std::cos(24*std::acos(2*x2 - 1)) - 0.0034493627577565598*std::cos(25*std::acos(2*x2 - 1)) + 0.0023656859807625087*std::cos(26*std::acos(2*x2 - 1)) - 0.0014683504293984927*std::cos(27*std::acos(2*x2 - 1)) + 0.00075787009367743463*std::cos(28*std::acos(2*x2 - 1)) - 0.00022211324503403421*std::cos(29*std::acos(2*x2 - 1)) - 0.00015895954194434722*std::cos(30*std::acos(2*x2 - 1)) + 0.00040942515509841473*std::cos(31*std::acos(2*x2 - 1)) - 0.00055444395230353979*std::cos(32*std::acos(2*x2 - 1)) + 0.00061816302915972121*std::cos(33*std::acos(2*x2 - 1)) - 0.0006223699751641866*std::cos(34*std::acos(2*x2 - 1)) + 0.00058575652376393825*std::cos(35*std::acos(2*x2 - 1)) - 0.00052364867584037827*std::cos(36*std::acos(2*x2 - 1)) + 0.00044806923918337335*std::cos(37*std::acos(2*x2 - 1)) - 0.00036801479075180416*std::cos(38*std::acos(2*x2 - 1)) + 0.00028985719739075058*std::cos(39*std::acos(2*x2 - 1)) - 0.00021779176589643498*std::cos(40*std::acos(2*x2 - 1)) + 0.00015429322229134138*std::cos(41*std::acos(2*x2 - 1)) - 0.00010052634987590958*std::cos(42*std::acos(2*x2 - 1)) + 5.6726279363716187e-5*std::cos(43*std::acos(2*x2 - 1)) - 2.2479144016110397e-5*std::cos(44*std::acos(2*x2 - 1)) - 2.9941138260901334e-6*std::cos(45*std::acos(2*x2 - 1)) + 2.0737540727223253e-5*std::cos(46*std::acos(2*x2 - 1)) - 3.1803645640635092e-5*std::cos(47*std::acos(2*x2 - 1)) + 3.7460730533993986e-5*std::cos(48*std::acos(2*x2 - 1)) - 1.9367680460562198e-5*std::cos(49*std::acos(2*x2 - 1)) + 0.78027563846425307;


         QU = -0.010576202918037571*x2 + 0.0049158484648276091*std::cos(2*std::acos(2*x2 - 1)) - 0.0043476546941215305*std::cos(3*std::acos(2*x2 - 1)) + 0.0036510468114268483*std::cos(4*std::acos(2*x2 - 1)) - 0.002901820194804759*std::cos(5*std::acos(2*x2 - 1)) + 0.0021711411472367381*std::cos(6*std::acos(2*x2 - 1)) - 0.0015152169044829356*std::cos(7*std::acos(2*x2 - 1)) + 0.00096945894743319112*std::cos(8*std::acos(2*x2 - 1)) - 0.00054764067133671469*std::cos(9*std::acos(2*x2 - 1)) + 0.00024524244637639153*std::cos(10*std::acos(2*x2 - 1)) - 4.5318379615411516e-5*std::cos(11*std::acos(2*x2 - 1)) - 7.5016437988135656e-5*std::cos(12*std::acos(2*x2 - 1)) + 0.00013903684892883959*std::cos(13*std::acos(2*x2 - 1)) - 0.00016669716221202129*std::cos(14*std::acos(2*x2 - 1)) + 0.00017293517850649356*std::cos(15*std::acos(2*x2 - 1)) - 0.00016762290100372283*std::cos(16*std::acos(2*x2 - 1)) + 0.0001565250979263642*std::cos(17*std::acos(2*x2 - 1)) - 0.00014262313149105034*std::cos(18*std::acos(2*x2 - 1)) + 0.00012733277994831236*std::cos(19*std::acos(2*x2 - 1)) - 0.00011137574712088311*std::cos(20*std::acos(2*x2 - 1)) + 9.5265922762152879e-5*std::cos(21*std::acos(2*x2 - 1)) - 7.949446679568027e-5*std::cos(22*std::acos(2*x2 - 1)) + 6.4539804546632901e-5*std::cos(23*std::acos(2*x2 - 1)) - 5.0812958869431718e-5*std::cos(24*std::acos(2*x2 - 1)) + 3.8607040414739756e-5*std::cos(25*std::acos(2*x2 - 1)) - 2.8078139633047628e-5*std::cos(26*std::acos(2*x2 - 1)) + 1.9256744034329973e-5*std::cos(27*std::acos(2*x2 - 1)) - 1.2076655967520002e-5*std::cos(28*std::acos(2*x2 - 1)) + 6.4078807450156464e-6*std::cos(29*std::acos(2*x2 - 1)) - 2.0848892851602414e-6*std::cos(30*std::acos(2*x2 - 1)) - 1.0728123079725384e-6*std::cos(31*std::acos(2*x2 - 1)) + 3.2472117308880848e-6*std::cos(32*std::acos(2*x2 - 1)) - 4.613394818354553e-6*std::cos(33*std::acos(2*x2 - 1)) + 5.3342909707103125e-6*std::cos(34*std::acos(2*x2 - 1)) - 5.5573164595291582e-6*std::cos(35*std::acos(2*x2 - 1)) + 5.4123264816095364e-6*std::cos(36*std::acos(2*x2 - 1)) - 5.0106566473464476e-6*std::cos(37*std::acos(2*x2 - 1)) + 4.4451399553703849e-6*std::cos(38*std::acos(2*x2 - 1)) - 3.7909502016114831e-6*std::cos(39*std::acos(2*x2 - 1)) + 3.1070609620702177e-6*std::cos(40*std::acos(2*x2 - 1)) - 2.4380795285656434e-6*std::cos(41*std::acos(2*x2 - 1)) + 1.8162317855278786e-6*std::cos(42*std::acos(2*x2 - 1)) - 1.2633287852501667e-6*std::cos(43*std::acos(2*x2 - 1)) + 7.9262328279254237e-7*std::cos(44*std::acos(2*x2 - 1)) - 4.1054722132961245e-7*std::cos(45*std::acos(2*x2 - 1)) + 1.1836464465185349e-7*std::cos(46*std::acos(2*x2 - 1)) + 8.6334464780438489e-8*std::cos(47*std::acos(2*x2 - 1)) - 2.0767142437350281e-7*std::cos(48*std::acos(2*x2 - 1)) + 1.2265058076128795e-7*std::cos(49*std::acos(2*x2 - 1)) + 0.0079969736601870848;


         QW = 0.0010502845245320962*x2 - 0.00047424806537354341*std::cos(2*std::acos(2*x2 - 1)) + 0.00039829187634152756*std::cos(3*std::acos(2*x2 - 1)) - 0.00030842359687639146*std::cos(4*std::acos(2*x2 - 1)) + 0.00021660529983876786*std::cos(5*std::acos(2*x2 - 1)) - 0.00013325542195377177*std::cos(6*std::acos(2*x2 - 1)) + 6.5554471688649202e-5*std::cos(7*std::acos(2*x2 - 1)) - 1.6744374021215435e-5*std::cos(8*std::acos(2*x2 - 1)) - 1.3577004399091652e-5*std::cos(9*std::acos(2*x2 - 1)) + 2.8441455188237235e-5*std::cos(10*std::acos(2*x2 - 1)) - 3.2192752270919155e-5*std::cos(11*std::acos(2*x2 - 1)) + 2.9267045537900673e-5*std::cos(12*std::acos(2*x2 - 1)) - 2.3340713022424632e-5*std::cos(13*std::acos(2*x2 - 1)) + 1.6937670016013351e-5*std::cos(14*std::acos(2*x2 - 1)) - 1.1434762864622835e-5*std::cos(15*std::acos(2*x2 - 1)) + 7.3168487622440586e-6*std::cos(16*std::acos(2*x2 - 1)) - 4.5194027501339433e-6*std::cos(17*std::acos(2*x2 - 1)) + 2.7357601113401332e-6*std::cos(18*std::acos(2*x2 - 1)) - 1.6267949384179162e-6*std::cos(19*std::acos(2*x2 - 1)) + 9.2539984341097565e-7*std::cos(20*std::acos(2*x2 - 1)) - 4.6152002544355246e-7*std::cos(21*std::acos(2*x2 - 1)) + 1.4388719977577337e-7*std::cos(22*std::acos(2*x2 - 1)) + 7.1170118434136188e-8*std::cos(23*std::acos(2*x2 - 1)) - 2.0598145767085363e-7*std::cos(24*std::acos(2*x2 - 1)) + 2.7679788704964439e-7*std::cos(25*std::acos(2*x2 - 1)) - 2.9943045696906068e-7*std::cos(26*std::acos(2*x2 - 1)) + 2.8938046693379158e-7*std::cos(27*std::acos(2*x2 - 1)) - 2.6020792292241089e-7*std::cos(28*std::acos(2*x2 - 1)) + 2.2228049320518512e-7*std::cos(29*std::acos(2*x2 - 1)) - 1.8253548245023164e-7*std::cos(30*std::acos(2*x2 - 1)) + 1.45021342981339e-7*std::cos(31*std::acos(2*x2 - 1)) - 1.1174116072256633e-7*std::cos(32*std::acos(2*x2 - 1)) + 8.3432870573577858e-8*std::cos(33*std::acos(2*x2 - 1)) - 6.0129328270964338e-8*std::cos(34*std::acos(2*x2 - 1)) + 4.149786794238191e-8*std::cos(35*std::acos(2*x2 - 1)) - 2.7026212008590692e-8*std::cos(36*std::acos(2*x2 - 1)) + 1.6123830743229749e-8*std::cos(37*std::acos(2*x2 - 1)) - 8.1834122387766739e-9*std::cos(38*std::acos(2*x2 - 1)) + 2.6227444681779583e-9*std::cos(39*std::acos(2*x2 - 1)) + 1.0868574376059125e-9*std::cos(40*std::acos(2*x2 - 1)) - 3.4042842347841957e-9*std::cos(41*std::acos(2*x2 - 1)) + 4.7132459814081889e-9*std::cos(42*std::acos(2*x2 - 1)) - 5.3253371243123665e-9*std::cos(43*std::acos(2*x2 - 1)) + 5.4884125576013118e-9*std::cos(44*std::acos(2*x2 - 1)) - 5.3969979275563877e-9*std::cos(45*std::acos(2*x2 - 1)) + 5.2016546705167162e-9*std::cos(46*std::acos(2*x2 - 1)) - 5.0151474394376526e-9*std::cos(47*std::acos(2*x2 - 1)) + 4.9156149738774846e-9*std::cos(48*std::acos(2*x2 - 1)) - 2.4225801557894737e-9*std::cos(49*std::acos(2*x2 - 1)) - 0.00079667276419180457;


         QDF = -0.010455066745063817*x2 + 0.0050593582905081105*std::cos(2*std::acos(2*x2 - 1)) - 0.004789708640902955*std::cos(3*std::acos(2*x2 - 1)) + 0.0044334519954988627*std::cos(4*std::acos(2*x2 - 1)) - 0.0040093350047312694*std::cos(5*std::acos(2*x2 - 1)) + 0.0035383502124445845*std::cos(6*std::acos(2*x2 - 1)) - 0.0030420063366992335*std::cos(7*std::acos(2*x2 - 1)) + 0.0025408298136973774*std::cos(8*std::acos(2*x2 - 1)) - 0.0020531178366578165*std::cos(9*std::acos(2*x2 - 1)) + 0.0015941244848603612*std::cos(10*std::acos(2*x2 - 1)) - 0.0011755842630321808*std::cos(11*std::acos(2*x2 - 1)) + 0.00080561997644330635*std::cos(12*std::acos(2*x2 - 1)) - 0.00048888474110406622*std::cos(13*std::acos(2*x2 - 1)) + 0.00022692949602540118*std::cos(14*std::acos(2*x2 - 1)) - 1.8663806887314557e-5*std::cos(15*std::acos(2*x2 - 1)) - 0.00013909961429497045*std::cos(16*std::acos(2*x2 - 1)) + 0.00025110577121937011*std::cos(17*std::acos(2*x2 - 1)) - 0.00032312354021896562*std::cos(18*std::acos(2*x2 - 1)) + 0.00036147095670710147*std::cos(19*std::acos(2*x2 - 1)) - 0.00037259451620934466*std::cos(20*std::acos(2*x2 - 1)) + 0.00036273171271154043*std::cos(21*std::acos(2*x2 - 1)) - 0.00033764766229243872*std::cos(22*std::acos(2*x2 - 1)) + 0.00030245608035755562*std::cos(23*std::acos(2*x2 - 1)) - 0.00026151076305032281*std::cos(24*std::acos(2*x2 - 1)) + 0.00021836451748329428*std::cos(25*std::acos(2*x2 - 1)) - 0.00017577923438368484*std::cos(26*std::acos(2*x2 - 1)) + 0.00013577781675913318*std::cos(27*std::acos(2*x2 - 1)) - 9.9723059267642776e-5*std::cos(28*std::acos(2*x2 - 1)) + 6.8413880502466061e-5*std::cos(29*std::acos(2*x2 - 1)) - 4.2188012292539325e-5*std::cos(30*std::acos(2*x2 - 1)) + 2.1024301111170907e-5*std::cos(31*std::acos(2*x2 - 1)) - 4.6382948636624143e-6*std::cos(32*std::acos(2*x2 - 1)) - 7.4323137019969709e-6*std::cos(33*std::acos(2*x2 - 1)) + 1.5754928967788635e-5*std::cos(34*std::acos(2*x2 - 1)) - 2.0942595605880331e-5*std::cos(35*std::acos(2*x2 - 1)) + 2.3607912279129951e-5*std::cos(36*std::acos(2*x2 - 1)) - 2.4329106457577498e-5*std::cos(37*std::acos(2*x2 - 1)) + 2.3627229235448358e-5*std::cos(38*std::acos(2*x2 - 1)) - 2.1952228750797814e-5*std::cos(39*std::acos(2*x2 - 1)) + 1.9677035100688809e-5*std::cos(40*std::acos(2*x2 - 1)) - 1.7096497602731732e-5*std::cos(41*std::acos(2*x2 - 1)) + 1.4432051229985118e-5*std::cos(42*std::acos(2*x2 - 1)) - 1.1836639755519902e-5*std::cos(43*std::acos(2*x2 - 1)) + 9.4058453027625118e-6*std::cos(44*std::acos(2*x2 - 1)) - 7.1820052475605402e-6*std::cos(45*std::acos(2*x2 - 1)) + 5.1722827654386992e-6*std::cos(46*std::acos(2*x2 - 1)) - 3.3414618186484205e-6*std::cos(47*std::acos(2*x2 - 1)) + 1.6486613198138318e-6*std::cos(48*std::acos(2*x2 - 1)) - 2.3148318570051458e-19*std::cos(49*std::acos(2*x2 - 1)) + 0.0078698798855119292;


         QDG = 0.0022610838262086923*x2 - 0.0009708196285129235*std::cos(2*std::acos(2*x2 - 1)) + 0.00073666490680725943*std::cos(3*std::acos(2*x2 - 1)) - 0.00046639573519408783*std::cos(4*std::acos(2*x2 - 1)) + 0.00019820283307177088*std::cos(5*std::acos(2*x2 - 1)) + 3.7966118276490607e-5*std::cos(6*std::acos(2*x2 - 1)) - 0.00022446385343084728*std::cos(7*std::acos(2*x2 - 1)) + 0.00035571988167206175*std::cos(8*std::acos(2*x2 - 1)) - 0.0004351177898990914*std::cos(9*std::acos(2*x2 - 1)) + 0.00047094834492425702*std::cos(10*std::acos(2*x2 - 1)) - 0.00047300334493550252*std::cos(11*std::acos(2*x2 - 1)) + 0.0004505060518301234*std::cos(12*std::acos(2*x2 - 1)) - 0.00041128801884293388*std::cos(13*std::acos(2*x2 - 1)) + 0.00036173037230146954*std::cos(14*std::acos(2*x2 - 1)) - 0.00030698029053774099*std::cos(15*std::acos(2*x2 - 1)) + 0.00025115095829012434*std::cos(16*std::acos(2*x2 - 1)) - 0.0001974244291607286*std::cos(17*std::acos(2*x2 - 1)) + 0.00014809868196024513*std::cos(18*std::acos(2*x2 - 1)) - 0.0001046441819414591*std::cos(19*std::acos(2*x2 - 1)) + 6.7804575796859103e-5*std::cos(20*std::acos(2*x2 - 1)) - 3.7737406585236285e-5*std::cos(21*std::acos(2*x2 - 1)) + 1.4169685055377826e-5*std::cos(22*std::acos(2*x2 - 1)) + 3.4573813680152854e-6*std::cos(23*std::acos(2*x2 - 1)) - 1.5869447232832877e-5*std::cos(24*std::acos(2*x2 - 1)) + 2.3869059742907e-5*std::cos(25*std::acos(2*x2 - 1)) - 2.8270704788192883e-5*std::cos(26*std::acos(2*x2 - 1)) + 2.9854948325385093e-5*std::cos(27*std::acos(2*x2 - 1)) - 2.9336857270458338e-5*std::cos(28*std::acos(2*x2 - 1)) + 2.7345700851255296e-5*std::cos(29*std::acos(2*x2 - 1)) - 2.4414125849614669e-5*std::cos(30*std::acos(2*x2 - 1)) + 2.0975342359990047e-5*std::cos(31*std::acos(2*x2 - 1)) - 1.7366718699566026e-5*std::cos(32*std::acos(2*x2 - 1)) + 1.3838004646320072e-5*std::cos(33*std::acos(2*x2 - 1)) - 1.0562420344022154e-5*std::cos(34*std::acos(2*x2 - 1)) + 7.6490811992196749e-6*std::cos(35*std::acos(2*x2 - 1)) - 5.1555789122350304e-6*std::cos(36*std::acos(2*x2 - 1)) + 3.099897370904226e-6*std::cos(37*std::acos(2*x2 - 1)) - 1.4711453061618453e-6*std::cos(38*std::acos(2*x2 - 1)) + 2.3881848083610448e-7*std::cos(39*std::acos(2*x2 - 1)) + 6.3952808264553131e-7*std::cos(40*std::acos(2*x2 - 1)) - 1.2121968864499953e-6*std::cos(41*std::acos(2*x2 - 1)) + 1.5287674546837017e-6*std::cos(42*std::acos(2*x2 - 1)) - 1.6368750482564979e-6*std::cos(43*std::acos(2*x2 - 1)) + 1.580170816606448e-6*std::cos(44*std::acos(2*x2 - 1)) - 1.3972690793098503e-6*std::cos(45*std::acos(2*x2 - 1)) + 1.121487899140177e-6*std::cos(46*std::acos(2*x2 - 1)) - 7.8117969216756329e-7*std::cos(47*std::acos(2*x2 - 1)) + 4.0053993378577471e-7*std::cos(48*std::acos(2*x2 - 1)) - 2.7873600226969817e-18*std::cos(49*std::acos(2*x2 - 1)) - 0.0017241936761803608;


         QDA = -0.0061636768851878005*x2 + 0.0028643144975169857*std::cos(2*std::acos(2*x2 - 1)) - 0.0025189357136866011*std::cos(3*std::acos(2*x2 - 1)) + 0.002069711893731192*std::cos(4*std::acos(2*x2 - 1)) - 0.0015470707660719242*std::cos(5*std::acos(2*x2 - 1)) + 0.00098506796605821981*std::cos(6*std::acos(2*x2 - 1)) - 0.00041847222701149882*std::cos(7*std::acos(2*x2 - 1)) - 0.00012010496168644909*std::cos(8*std::acos(2*x2 - 1)) + 0.00060260862248479648*std::cos(9*std::acos(2*x2 - 1)) - 0.0010073804561400277*std::cos(10*std::acos(2*x2 - 1)) + 0.0013201952007468913*std::cos(11*std::acos(2*x2 - 1)) - 0.0015345685858873157*std::cos(12*std::acos(2*x2 - 1)) + 0.0016513051016688034*std::cos(13*std::acos(2*x2 - 1)) - 0.0016774995156244403*std::cos(14*std::acos(2*x2 - 1)) + 0.0016251024236868632*std::cos(15*std::acos(2*x2 - 1)) - 0.0015093060941546771*std::cos(16*std::acos(2*x2 - 1)) + 0.0013468906348878078*std::cos(17*std::acos(2*x2 - 1)) - 0.0011547303381543304*std::cos(18*std::acos(2*x2 - 1)) + 0.00094854370747171287*std::cos(19*std::acos(2*x2 - 1)) - 0.00074198115327222646*std::cos(20*std::acos(2*x2 - 1)) + 0.00054605030788216135*std::cos(21*std::acos(2*x2 - 1)) - 0.0003688763364239934*std::cos(22*std::acos(2*x2 - 1)) + 0.00021573688931837251*std::cos(23*std::acos(2*x2 - 1)) - 8.9316154011789429e-5*std::cos(24*std::acos(2*x2 - 1)) - 9.9002166939462696e-6*std::cos(25*std::acos(2*x2 - 1)) + 8.3152070110157757e-5*std::cos(26*std::acos(2*x2 - 1)) - 0.00013291594920144509*std::cos(27*std::acos(2*x2 - 1)) + 0.0001624430329431611*std::cos(28*std::acos(2*x2 - 1)) - 0.00017535669566776365*std::cos(29*std::acos(2*x2 - 1)) + 0.00017532560759641645*std::cos(30*std::acos(2*x2 - 1)) - 0.00016581911838231834*std::cos(31*std::acos(2*x2 - 1)) + 0.00014994123528813618*std::cos(32*std::acos(2*x2 - 1)) - 0.00013033468769269095*std::cos(33*std::acos(2*x2 - 1)) + 0.00010914245496420367*std::cos(34*std::acos(2*x2 - 1)) - 8.8013508220932649e-5*std::cos(35*std::acos(2*x2 - 1)) + 6.8139477799246748e-5*std::cos(36*std::acos(2*x2 - 1)) - 5.031078470906537e-5*std::cos(37*std::acos(2*x2 - 1)) + 3.4982337521102774e-5*std::cos(38*std::acos(2*x2 - 1)) - 2.2341646627357292e-5*std::cos(39*std::acos(2*x2 - 1)) + 1.23734648201926e-5*std::cos(40*std::acos(2*x2 - 1)) - 4.918296643964379e-6*std::cos(41*std::acos(2*x2 - 1)) - 2.7858573371995793e-7*std::cos(42*std::acos(2*x2 - 1)) + 3.5259209677190364e-6*std::cos(43*std::acos(2*x2 - 1)) - 5.1570481982314353e-6*std::cos(44*std::acos(2*x2 - 1)) + 5.5040866302193387e-6*std::cos(45*std::acos(2*x2 - 1)) - 4.8875778191506775e-6*std::cos(46*std::acos(2*x2 - 1)) + 3.5962319293774782e-6*std::cos(47*std::acos(2*x2 - 1)) - 1.898033634153107e-6*std::cos(48*std::acos(2*x2 - 1)) - 1.7372769825277583e-18*std::cos(49*std::acos(2*x2 - 1)) + 0.0046598928471686864;
        
        }
        
        else
        {
            QF = 1;
            QG = 0;
            QA = 1;
            QU = 0;
            QW = 0;
            QDF = 0;
            QDG = 0;
            QDA = 0;
        }
        
        amrex::Real u = QU;
        
        state_fab(i, j, k, Idx::u) = u;
        
        amrex::Real psi = 1 + u;
        
        state_fab(i, j, k, Idx::w_U_0) = -y/R*QW;
        state_fab(i, j, k, Idx::w_U_1) = x/R*QW;
        state_fab(i, j, k, Idx::w_U_2) = x/R*QW; //0;
        state_fab(i, j, k, Idx::chi) = 1/(psi*psi*psi*psi);
        
        
        amrex::Real dthetadx = -x*z/(std::sqrt(rdomain*rdomain - z*z)*(R*R + rdomain*rdomain - 2*R*std::sqrt(rdomain*rdomain - z*z)));
        amrex::Real dthetady = -y*z/(std::sqrt(rdomain*rdomain - z*z)*(R*R + rdomain*rdomain - 2*R*std::sqrt(rdomain*rdomain - z*z)));
        
        amrex::Real f1 = -R + std::sqrt(rdomain*rdomain - z*z);
        amrex::Real f2 = 1 + z*z/((R - std::sqrt(x*x+y*y))*(R - std::sqrt(x*x+y*y)));
        
        amrex::Real dthetadz = 1.0/(f1*f2);

        amrex::Real A_L_theta = QA/e;
        
        state_fab(i, j, k, Idx::A_L_0) = dthetadx*A_L_theta;
        state_fab(i, j, k, Idx::A_L_1) = dthetady*A_L_theta;
        state_fab(i, j, k, Idx::A_L_2) = dthetadz*A_L_theta;
        
        amrex::Real costheta = (rdomain*rdomain-R*R-r*r)/(2*R*r);//(std::sqrt(x*x+y*y) - R)/r;
        amrex::Real sintheta = z/r;
        
        state_fab(i, j, k, Idx::phi1) = std::sqrt(2)*eta*QF*costheta;
        state_fab(i, j, k, Idx::phi2) = std::sqrt(2)*eta*QF*sintheta; 
        
        amrex::Real cosphi = x/rdomain;
        amrex::Real sinphi = y/rdomain;
        
        state_fab(i, j, k, Idx::sigma1) = std::sqrt(2)*eta*QG*(cosphi*cosphi - sinphi*sinphi);
        state_fab(i, j, k, Idx::sigma2) = -std::sqrt(2)*eta*QG*2*cosphi*sinphi;
        
        state_fab(i, j, k, Idx::dDtsigma1) = std::sqrt(2)*omega*eta*QG*2*cosphi*sinphi;
        state_fab(i, j, k, Idx::dDtsigma2) = std::sqrt(2)*omega*eta*QG*(cosphi*cosphi - sinphi*sinphi);
        
        amrex::Real V = 0.25*lambda*std::pow(eta,4)*std::pow((QF*QF+QG*QG-1),2) + eta*eta*OMEGA*OMEGA*QG*QG;
        amrex::Real A = eta*eta/(r*r)*(QF*QF*(1-QA)*(1-QA) + r*r*QDF*QDF);
        amrex::Real B = eta*eta*(R*R*omega*omega*QG*QG/(x*x+y*y) + QDG*QDG);
        amrex::Real C = QDA*QDA/(2*e*e*r*r);
        
        state_fab(i, j, k, Idx::V) = V;
        state_fab(i, j, k, Idx::A) = A;
        state_fab(i, j, k, Idx::B) = B;
        state_fab(i, j, k, Idx::C) = C;
        
        state_fab(i, j, k, Idx::S_U_phi) = R*eta*eta*omega*omega*QG*QG/(psi*psi*(x_U_0*x_U_0+x_U_1*x_U_1));
        
}
#endif