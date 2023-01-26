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
        
        amrex::Real x2 = r/200;
        
        amrex::Real QF, QG, QA, QU, QW, QDF, QDG, QDA;
        
        //amrex::Real x2 = std::sqrt(x_U_0*x_U_0 + x_U_1*x_U_1 + x_U_2*x_U_2)/1000;
        
        if(r < 200.0)
        {
        
            QF = 0.63832545161588994*x2 - 0.24281663099137327*std::cos(2*std::acos(2*x2 - 1)) + 0.15483080215801073*std::cos(3*std::acos(2*x2 - 1)) - 0.081425371081413261*std::cos(4*std::acos(2*x2 - 1)) + 0.032865613175971795*std::cos(5*std::acos(2*x2 - 1)) - 0.0071875705931113316*std::cos(6*std::acos(2*x2 - 1)) - 0.0029017640706157609*std::cos(7*std::acos(2*x2 - 1)) + 0.004780745239444998*std::cos(8*std::acos(2*x2 - 1)) - 0.0035491316257708176*std::cos(9*std::acos(2*x2 - 1)) + 0.0018571419660234657*std::cos(10*std::acos(2*x2 - 1)) - 0.00066084558627966042*std::cos(11*std::acos(2*x2 - 1)) + 5.8959930942702188e-5*std::cos(12*std::acos(2*x2 - 1)) + 0.00014153824540662219*std::cos(13*std::acos(2*x2 - 1)) - 0.00014807289084322467*std::cos(14*std::acos(2*x2 - 1)) + 9.4741358412878118e-5*std::cos(15*std::acos(2*x2 - 1)) - 4.3695732001240717e-5*std::cos(16*std::acos(2*x2 - 1)) + 1.2487257971306037e-5*std::cos(17*std::acos(2*x2 - 1)) + 1.1267418540857055e-6*std::cos(18*std::acos(2*x2 - 1)) - 4.4137124135187888e-6*std::cos(19*std::acos(2*x2 - 1)) + 3.4046406570710508e-6*std::cos(20*std::acos(2*x2 - 1)) - 1.5014092017602565e-6*std::cos(21*std::acos(2*x2 - 1)) + 4.9789861306798189e-8*std::cos(22*std::acos(2*x2 - 1)) + 7.1949500441217951e-7*std::cos(23*std::acos(2*x2 - 1)) - 9.9340893749272352e-7*std::cos(24*std::acos(2*x2 - 1)) + 1.0130029197905066e-6*std::cos(25*std::acos(2*x2 - 1)) - 9.4515870973158056e-7*std::cos(26*std::acos(2*x2 - 1)) + 8.7225446104496895e-7*std::cos(27*std::acos(2*x2 - 1)) - 8.2071946835501897e-7*std::cos(28*std::acos(2*x2 - 1)) + 7.9045245171202594e-7*std::cos(29*std::acos(2*x2 - 1)) - 7.7349017757813443e-7*std::cos(30*std::acos(2*x2 - 1)) + 7.624999500480123e-7*std::cos(31*std::acos(2*x2 - 1)) - 7.530723172902154e-7*std::cos(32*std::acos(2*x2 - 1)) + 7.4331310158728958e-7*std::cos(33*std::acos(2*x2 - 1)) - 7.3278392172596019e-7*std::cos(34*std::acos(2*x2 - 1)) + 7.2164393282239669e-7*std::cos(35*std::acos(2*x2 - 1)) - 7.1017145816462679e-7*std::cos(36*std::acos(2*x2 - 1)) + 6.9857703461490006e-7*std::cos(37*std::acos(2*x2 - 1)) - 6.8697289319630261e-7*std::cos(38*std::acos(2*x2 - 1)) + 6.754003452520696e-7*std::cos(39*std::acos(2*x2 - 1)) - 6.63863996226699e-7*std::cos(40*std::acos(2*x2 - 1)) + 6.5235491080936625e-7*std::cos(41*std::acos(2*x2 - 1)) - 6.4086153155409886e-7*std::cos(42*std::acos(2*x2 - 1)) + 6.2937295100887192e-7*std::cos(43*std::acos(2*x2 - 1)) - 6.1787961273818137e-7*std::cos(44*std::acos(2*x2 - 1)) + 6.0637492581630537e-7*std::cos(45*std::acos(2*x2 - 1)) - 5.9485973353048949e-7*std::cos(46*std::acos(2*x2 - 1)) + 5.8334438458526594e-7*std::cos(47*std::acos(2*x2 - 1)) - 5.7197000105646545e-7*std::cos(48*std::acos(2*x2 - 1)) + 2.8031418982006617e-7*std::cos(49*std::acos(2*x2 - 1)) + 0.50576669238477256;


         QG = -0.40002566742499091*x2 + 0.07354723999498175*std::cos(2*std::acos(2*x2 - 1)) - 0.013885391878921453*std::cos(3*std::acos(2*x2 - 1)) - 0.0061672912129416323*std::cos(4*std::acos(2*x2 - 1)) + 0.0082672260890876252*std::cos(5*std::acos(2*x2 - 1)) - 0.005500995192996164*std::cos(6*std::acos(2*x2 - 1)) + 0.0026127350396768904*std::cos(7*std::acos(2*x2 - 1)) - 0.00082029811045166399*std::cos(8*std::acos(2*x2 - 1)) + 1.4317577076259578e-5*std::cos(9*std::acos(2*x2 - 1)) + 0.00021234117742739123*std::cos(10*std::acos(2*x2 - 1)) - 0.00019268605086088207*std::cos(11*std::acos(2*x2 - 1)) + 0.0001126948269211581*std::cos(12*std::acos(2*x2 - 1)) - 4.6847806806569396e-5*std::cos(13*std::acos(2*x2 - 1)) + 1.0320090184773145e-5*std::cos(14*std::acos(2*x2 - 1)) + 3.9468948778253991e-6*std::cos(15*std::acos(2*x2 - 1)) - 6.4763105399310573e-6*std::cos(16*std::acos(2*x2 - 1)) + 4.7214906553043272e-6*std::cos(17*std::acos(2*x2 - 1)) - 2.422044880268252e-6*std::cos(18*std::acos(2*x2 - 1)) + 8.4391792233360973e-7*std::cos(19*std::acos(2*x2 - 1)) - 7.1775357405774803e-8*std::cos(20*std::acos(2*x2 - 1)) - 1.7667143929760903e-7*std::cos(21*std::acos(2*x2 - 1)) + 1.8077203058236056e-7*std::cos(22*std::acos(2*x2 - 1)) - 1.1308458628514921e-7*std::cos(23*std::acos(2*x2 - 1)) + 5.029529421780045e-8*std::cos(24*std::acos(2*x2 - 1)) - 1.2619282709429076e-8*std::cos(25*std::acos(2*x2 - 1)) - 3.5835187816146731e-9*std::cos(26*std::acos(2*x2 - 1)) + 7.4500938739788094e-9*std::cos(27*std::acos(2*x2 - 1)) - 6.2645002644248624e-9*std::cos(28*std::acos(2*x2 - 1)) + 4.0372419696954667e-9*std::cos(29*std::acos(2*x2 - 1)) - 2.3251834499321934e-9*std::cos(30*std::acos(2*x2 - 1)) + 1.3895086481398206e-9*std::cos(31*std::acos(2*x2 - 1)) - 1.0169163627692044e-9*std::cos(32*std::acos(2*x2 - 1)) + 9.3563726739756398e-10*std::cos(33*std::acos(2*x2 - 1)) - 9.582770190485269e-10*std::cos(34*std::acos(2*x2 - 1)) + 9.9193162980237572e-10*std::cos(35*std::acos(2*x2 - 1)) - 1.0058990075076568e-9*std::cos(36*std::acos(2*x2 - 1)) + 9.9906967030754981e-10*std::cos(37*std::acos(2*x2 - 1)) - 9.7936370867632739e-10*std::cos(38*std::acos(2*x2 - 1)) + 9.5430548649953187e-10*std::cos(39*std::acos(2*x2 - 1)) - 9.2833855472775542e-10*std::cos(40*std::acos(2*x2 - 1)) + 9.0311350952126541e-10*std::cos(41*std::acos(2*x2 - 1)) - 8.7859720522765055e-10*std::cos(42*std::acos(2*x2 - 1)) + 8.5419827328156751e-10*std::cos(43*std::acos(2*x2 - 1)) - 8.2981844153499825e-10*std::cos(44*std::acos(2*x2 - 1)) + 8.0701773826883199e-10*std::cos(45*std::acos(2*x2 - 1)) - 7.8984265246339624e-10*std::cos(46*std::acos(2*x2 - 1)) + 7.8437381799921751e-10*std::cos(47*std::acos(2*x2 - 1)) - 7.8871191166228095e-10*std::cos(48*std::acos(2*x2 - 1)) + 3.919597270517672e-10*std::cos(49*std::acos(2*x2 - 1)) + 0.34186183186843488;


         QA = 0.83009278287297616*x2 - 0.27768966638640791*std::cos(2*std::acos(2*x2 - 1)) + 0.12988395790027363*std::cos(3*std::acos(2*x2 - 1)) - 0.025476745094575908*std::cos(4*std::acos(2*x2 - 1)) - 0.022638566648801792*std::cos(5*std::acos(2*x2 - 1)) + 0.03078585265727548*std::cos(6*std::acos(2*x2 - 1)) - 0.021565503051391311*std::cos(7*std::acos(2*x2 - 1)) + 0.010205145917913282*std::cos(8*std::acos(2*x2 - 1)) - 0.0026836626675105883*std::cos(9*std::acos(2*x2 - 1)) - 0.00066093562670769754*std::cos(10*std::acos(2*x2 - 1)) + 0.0014161456803043829*std::cos(11*std::acos(2*x2 - 1)) - 0.0010980053034770343*std::cos(12*std::acos(2*x2 - 1)) + 0.00058303455664891757*std::cos(13*std::acos(2*x2 - 1)) - 0.00020769891405703414*std::cos(14*std::acos(2*x2 - 1)) + 1.7233077009546488e-5*std::cos(15*std::acos(2*x2 - 1)) + 4.6014724167469517e-5*std::cos(16*std::acos(2*x2 - 1)) - 4.7580946058216333e-5*std::cos(17*std::acos(2*x2 - 1)) + 3.0228536556563268e-5*std::cos(18*std::acos(2*x2 - 1)) - 1.3758863929260595e-5*std::cos(19*std::acos(2*x2 - 1)) + 3.7123083962628534e-6*std::cos(20*std::acos(2*x2 - 1)) + 6.7184589790103555e-7*std::cos(21*std::acos(2*x2 - 1)) - 1.7347705447007256e-6*std::cos(22*std::acos(2*x2 - 1)) + 1.4108732114236078e-6*std::cos(23*std::acos(2*x2 - 1)) - 7.9229180882656181e-7*std::cos(24*std::acos(2*x2 - 1)) + 3.1429423871104367e-7*std::cos(25*std::acos(2*x2 - 1)) - 5.4909599252889185e-8*std::cos(26*std::acos(2*x2 - 1)) - 4.3288381972733658e-8*std::cos(27*std::acos(2*x2 - 1)) + 5.6946846863081824e-8*std::cos(28*std::acos(2*x2 - 1)) - 4.0174599057196669e-8*std::cos(29*std::acos(2*x2 - 1)) + 2.039988261301473e-8*std::cos(30*std::acos(2*x2 - 1)) - 6.9665619314984717e-9*std::cos(31*std::acos(2*x2 - 1)) + 3.6086793622769342e-10*std::cos(32*std::acos(2*x2 - 1)) + 1.7700425049282931e-9*std::cos(33*std::acos(2*x2 - 1)) - 1.7711139881082527e-9*std::cos(34*std::acos(2*x2 - 1)) + 1.1270894009762777e-9*std::cos(35*std::acos(2*x2 - 1)) - 5.2086819579688637e-10*std::cos(36*std::acos(2*x2 - 1)) + 1.4428371691854736e-10*std::cos(37*std::acos(2*x2 - 1)) + 2.8368057902738997e-11*std::cos(38*std::acos(2*x2 - 1)) - 7.7519549807059107e-11*std::cos(39*std::acos(2*x2 - 1)) + 6.8499098754282173e-11*std::cos(40*std::acos(2*x2 - 1)) - 3.6414340293111636e-11*std::cos(41*std::acos(2*x2 - 1)) - 4.7110475243083982e-12*std::cos(42*std::acos(2*x2 - 1)) + 4.6724717844837116e-11*std::cos(43*std::acos(2*x2 - 1)) - 7.6402054821622656e-11*std::cos(44*std::acos(2*x2 - 1)) + 7.5338839333749519e-11*std::cos(45*std::acos(2*x2 - 1)) - 3.2255934034886025e-11*std::cos(46*std::acos(2*x2 - 1)) - 4.0940244250275981e-11*std::cos(47*std::acos(2*x2 - 1)) + 1.0701985825312121e-10*std::cos(48*std::acos(2*x2 - 1)) - 6.5443634886253044e-11*std::cos(49*std::acos(2*x2 - 1)) + 0.34901821221024826;


         QU = -0.032462779159564061*x2 + 0.0057686066256195535*std::cos(2*std::acos(2*x2 - 1)) - 0.0013308775612759578*std::cos(3*std::acos(2*x2 - 1)) - 0.0001768086140290021*std::cos(4*std::acos(2*x2 - 1)) + 0.00050554972132949732*std::cos(5*std::acos(2*x2 - 1)) - 0.00043712975898278279*std::cos(6*std::acos(2*x2 - 1)) + 0.00027726869191181799*std::cos(7*std::acos(2*x2 - 1)) - 0.00013712925560159663*std::cos(8*std::acos(2*x2 - 1)) + 4.6509070987108401e-5*std::cos(9*std::acos(2*x2 - 1)) - 1.0944398552793399e-6*std::cos(10*std::acos(2*x2 - 1)) - 1.4364281079548065e-5*std::cos(11*std::acos(2*x2 - 1)) + 1.4613401812068561e-5*std::cos(12*std::acos(2*x2 - 1)) - 9.6730868675527835e-6*std::cos(13*std::acos(2*x2 - 1)) + 4.6863451133043e-6*std::cos(14*std::acos(2*x2 - 1)) - 1.4069554001756526e-6*std::cos(15*std::acos(2*x2 - 1)) - 1.6388026577831143e-7*std::cos(16*std::acos(2*x2 - 1)) + 6.1979159538619745e-7*std::cos(17*std::acos(2*x2 - 1)) - 5.4853156458008888e-7*std::cos(18*std::acos(2*x2 - 1)) + 3.3007626487678715e-7*std::cos(19*std::acos(2*x2 - 1)) - 1.4143839432619532e-7*std::cos(20*std::acos(2*x2 - 1)) + 2.9880715046503419e-8*std::cos(21*std::acos(2*x2 - 1)) + 1.6758589920740784e-8*std::cos(22*std::acos(2*x2 - 1)) - 2.5700605158729742e-8*std::cos(23*std::acos(2*x2 - 1)) + 1.9296057218326294e-8*std::cos(24*std::acos(2*x2 - 1)) - 1.0276351670984562e-8*std::cos(25*std::acos(2*x2 - 1)) + 3.676673327653148e-9*std::cos(26*std::acos(2*x2 - 1)) - 2.3635540205280459e-10*std::cos(27*std::acos(2*x2 - 1)) - 9.5143593269579407e-10*std::cos(28*std::acos(2*x2 - 1)) + 9.9200646973087758e-10*std::cos(29*std::acos(2*x2 - 1)) - 6.5325984266222992e-10*std::cos(30*std::acos(2*x2 - 1)) + 3.1692061585096425e-10*std::cos(31*std::acos(2*x2 - 1)) - 1.0331883628076375e-10*std::cos(32*std::acos(2*x2 - 1)) + 6.1510670757172196e-12*std::cos(33*std::acos(2*x2 - 1)) + 1.8688655203464893e-11*std::cos(34*std::acos(2*x2 - 1)) - 1.1318507355213757e-11*std::cos(35*std::acos(2*x2 - 1)) - 3.6527599552995461e-12*std::cos(36*std::acos(2*x2 - 1)) + 1.4777688020857136e-11*std::cos(37*std::acos(2*x2 - 1)) - 1.797221029824143e-11*std::cos(38*std::acos(2*x2 - 1)) + 1.2168045896494468e-11*std::cos(39*std::acos(2*x2 - 1)) + 2.0761192361000907e-12*std::cos(40*std::acos(2*x2 - 1)) - 2.1228227012827586e-11*std::cos(41*std::acos(2*x2 - 1)) + 3.6017354688002311e-11*std::cos(42*std::acos(2*x2 - 1)) - 2.9708638414006943e-11*std::cos(43*std::acos(2*x2 - 1)) - 1.9219502549831085e-11*std::cos(44*std::acos(2*x2 - 1)) + 1.2735239128244105e-10*std::cos(45*std::acos(2*x2 - 1)) - 2.9130082788975327e-10*std::cos(46*std::acos(2*x2 - 1)) + 4.769561316323453e-10*std::cos(47*std::acos(2*x2 - 1)) - 6.2746052302492222e-10*std::cos(48*std::acos(2*x2 - 1)) + 3.3935029850405851e-10*std::cos(49*std::acos(2*x2 - 1)) + 0.027953900226934739;


         QW = 0.0016947878157233594*x2 - 0.00014354793823375775*std::cos(2*std::acos(2*x2 - 1)) - 2.8680583607548555e-5*std::cos(3*std::acos(2*x2 - 1)) + 3.5278491865887772e-5*std::cos(4*std::acos(2*x2 - 1)) - 1.7353103834272619e-5*std::cos(5*std::acos(2*x2 - 1)) + 5.748499948477433e-6*std::cos(6*std::acos(2*x2 - 1)) - 1.0517962948302559e-6*std::cos(7*std::acos(2*x2 - 1)) - 2.6574956932456781e-7*std::cos(8*std::acos(2*x2 - 1)) + 4.01804814541524e-7*std::cos(9*std::acos(2*x2 - 1)) - 2.5934440928522685e-7*std::cos(10*std::acos(2*x2 - 1)) + 1.2120457787633468e-7*std::cos(11*std::acos(2*x2 - 1)) - 4.0764403466583811e-8*std::cos(12*std::acos(2*x2 - 1)) + 5.7133146813596692e-9*std::cos(13*std::acos(2*x2 - 1)) + 4.8767699116175952e-9*std::cos(14*std::acos(2*x2 - 1)) - 5.5079046442003578e-9*std::cos(15*std::acos(2*x2 - 1)) + 3.5023346557453274e-9*std::cos(16*std::acos(2*x2 - 1)) - 1.6384491771985749e-9*std::cos(17*std::acos(2*x2 - 1)) + 5.3374622712809147e-10*std::cos(18*std::acos(2*x2 - 1)) - 4.6112698240650108e-11*std::cos(19*std::acos(2*x2 - 1)) - 9.7300918956744538e-11*std::cos(20*std::acos(2*x2 - 1)) + 9.7167020831090641e-11*std::cos(21*std::acos(2*x2 - 1)) - 5.9547109559109854e-11*std::cos(22*std::acos(2*x2 - 1)) + 2.6192995983862365e-11*std::cos(23*std::acos(2*x2 - 1)) - 6.5004361721939262e-12*std::cos(24*std::acos(2*x2 - 1)) - 2.1308033072558781e-12*std::cos(25*std::acos(2*x2 - 1)) + 4.5133853039964247e-12*std::cos(26*std::acos(2*x2 - 1)) - 4.2692606647910361e-12*std::cos(27*std::acos(2*x2 - 1)) + 3.3612834291167183e-12*std::cos(28*std::acos(2*x2 - 1)) - 2.5772857213705724e-12*std::cos(29*std::acos(2*x2 - 1)) + 2.0925840515263549e-12*std::cos(30*std::acos(2*x2 - 1)) - 1.8482674392679183e-12*std::cos(31*std::acos(2*x2 - 1)) + 1.7413113514252269e-12*std::cos(32*std::acos(2*x2 - 1)) - 1.693663139917624e-12*std::cos(33*std::acos(2*x2 - 1)) + 1.6625926633651639e-12*std::cos(34*std::acos(2*x2 - 1)) - 1.6308531818010829e-12*std::cos(35*std::acos(2*x2 - 1)) + 1.5946437628328008e-12*std::cos(36*std::acos(2*x2 - 1)) - 1.5552905356818035e-12*std::cos(37*std::acos(2*x2 - 1)) + 1.5150281238429068e-12*std::cos(38*std::acos(2*x2 - 1)) - 1.4754542526129494e-12*std::cos(39*std::acos(2*x2 - 1)) + 1.4373444103145308e-12*std::cos(40*std::acos(2*x2 - 1)) - 1.4009671214333828e-12*std::cos(41*std::acos(2*x2 - 1)) + 1.3664201596244991e-12*std::cos(42*std::acos(2*x2 - 1)) - 1.3338387216066098e-12*std::cos(43*std::acos(2*x2 - 1)) + 1.3033449938940597e-12*std::cos(44*std::acos(2*x2 - 1)) - 1.2746338332745732e-12*std::cos(45*std::acos(2*x2 - 1)) + 1.2466282932896695e-12*std::cos(46*std::acos(2*x2 - 1)) - 1.2193505431076839e-12*std::cos(47*std::acos(2*x2 - 1)) + 1.200275369791135e-12*std::cos(48*std::acos(2*x2 - 1)) - 5.9113541030317835e-13*std::cos(49*std::acos(2*x2 - 1)) - 0.0015451459303222033;


         QDF = -0.031998108479205069*x2 + 0.011459874447863821*std::cos(2*std::acos(2*x2 - 1)) - 0.006286340436378564*std::cos(3*std::acos(2*x2 - 1)) + 0.0021699798689102643*std::cos(4*std::acos(2*x2 - 1)) + 0.00022772182044586851*std::cos(5*std::acos(2*x2 - 1)) - 0.0011165978815756014*std::cos(6*std::acos(2*x2 - 1)) + 0.0010902346041831289*std::cos(7*std::acos(2*x2 - 1)) - 0.00071034888044436539*std::cos(8*std::acos(2*x2 - 1)) + 0.00032531154125657488*std::cos(9*std::acos(2*x2 - 1)) - 7.1501993571265361e-5*std::cos(10*std::acos(2*x2 - 1)) - 4.6118709099402732e-5*std::cos(11*std::acos(2*x2 - 1)) + 7.3884762344074518e-5*std::cos(12*std::acos(2*x2 - 1)) - 6.0269163277888865e-5*std::cos(13*std::acos(2*x2 - 1)) + 3.7084634537760659e-5*std::cos(14*std::acos(2*x2 - 1)) - 1.8808546538719765e-5*std::cos(15*std::acos(2*x2 - 1)) + 8.6620849010684443e-6*std::cos(16*std::acos(2*x2 - 1)) - 4.8258423847794271e-6*std::cos(17*std::acos(2*x2 - 1)) + 4.4163959623549874e-6*std::cos(18*std::acos(2*x2 - 1)) - 5.2314714803789851e-6*std::cos(19*std::acos(2*x2 - 1)) + 6.0936150656411179e-6*std::cos(20*std::acos(2*x2 - 1)) - 6.5933345525413278e-6*std::cos(21*std::acos(2*x2 - 1)) + 6.724210083337082e-6*std::cos(22*std::acos(2*x2 - 1)) - 6.6152422010153075e-6*std::cos(23*std::acos(2*x2 - 1)) + 6.3932407263712035e-6*std::cos(24*std::acos(2*x2 - 1)) - 6.1384035269929266e-6*std::cos(25*std::acos(2*x2 - 1)) + 5.8867367340234412e-6*std::cos(26*std::acos(2*x2 - 1)) - 5.6469185404305104e-6*std::cos(27*std::acos(2*x2 - 1)) + 5.4157169698642596e-6*std::cos(28*std::acos(2*x2 - 1)) - 5.1873133401674855e-6*std::cos(29*std::acos(2*x2 - 1)) + 4.9572522556310337e-6*std::cos(30*std::acos(2*x2 - 1)) - 4.7232169131963217e-6*std::cos(31*std::acos(2*x2 - 1)) + 4.4844999227613458e-6*std::cos(32*std::acos(2*x2 - 1)) - 4.2412482202057334e-6*std::cos(33*std::acos(2*x2 - 1)) + 3.9939108228196904e-6*std::cos(34*std::acos(2*x2 - 1)) - 3.7429526619564989e-6*std::cos(35*std::acos(2*x2 - 1)) + 3.4887575440402829e-6*std::cos(36*std::acos(2*x2 - 1)) - 3.2316266553914981e-6*std::cos(37*std::acos(2*x2 - 1)) + 2.9718079536960616e-6*std::cos(38*std::acos(2*x2 - 1)) - 2.7095246461343232e-6*std::cos(39*std::acos(2*x2 - 1)) + 2.4449930502831486e-6*std::cos(40*std::acos(2*x2 - 1)) - 2.1784307936888635e-6*std::cos(41*std::acos(2*x2 - 1)) + 1.9100593488181098e-6*std::cos(42*std::acos(2*x2 - 1)) - 1.6401044155506199e-6*std::cos(43*std::acos(2*x2 - 1)) + 1.3687959046322676e-6*std::cos(44*std::acos(2*x2 - 1)) - 1.0963676376060921e-6*std::cos(45*std::acos(2*x2 - 1)) + 8.2305574269520544e-7*std::cos(46*std::acos(2*x2 - 1)) - 5.4909394640388121e-7*std::cos(47*std::acos(2*x2 - 1)) + 2.7470927951541757e-7*std::cos(48*std::acos(2*x2 - 1)) + 3.1763556851023951e-17*std::cos(49*std::acos(2*x2 - 1)) + 0.024920634679829995;


         QDG = 0.0034533396356233857*x2 + 0.00031038023468549789*std::cos(2*std::acos(2*x2 - 1)) - 0.0012152344915091187*std::cos(3*std::acos(2*x2 - 1)) + 0.0011435079130591511*std::cos(4*std::acos(2*x2 - 1)) - 0.00072184872754495178*std::cos(5*std::acos(2*x2 - 1)) + 0.00031678117051666803*std::cos(6*std::acos(2*x2 - 1)) - 6.1726003771822222e-5*std::cos(7*std::acos(2*x2 - 1)) - 4.9003563961764426e-5*std::cos(8*std::acos(2*x2 - 1)) + 6.9522350142215477e-5*std::cos(9*std::acos(2*x2 - 1)) - 5.158074072138169e-5*std::cos(10*std::acos(2*x2 - 1)) + 2.7053902314504832e-5*std::cos(11*std::acos(2*x2 - 1)) - 9.1895975762492002e-6*std::cos(12*std::acos(2*x2 - 1)) + 7.0086189491417342e-9*std::cos(13*std::acos(2*x2 - 1)) + 2.9908930958917695e-6*std::cos(14*std::acos(2*x2 - 1)) - 2.8826310809831073e-6*std::cos(15*std::acos(2*x2 - 1)) + 1.8068187121898831e-6*std::cos(16*std::acos(2*x2 - 1)) - 8.1020134604510355e-7*std::cos(17*std::acos(2*x2 - 1)) + 2.0150386281838664e-7*std::cos(18*std::acos(2*x2 - 1)) + 6.173917053126857e-8*std::cos(19*std::acos(2*x2 - 1)) - 1.1918655109520024e-7*std::cos(20*std::acos(2*x2 - 1)) + 9.0449457030373074e-8*std::cos(21*std::acos(2*x2 - 1)) - 4.4984175594103272e-8*std::cos(22*std::acos(2*x2 - 1)) + 1.0909365883611386e-8*std::cos(23*std::acos(2*x2 - 1)) + 7.0349941761233369e-9*std::cos(24*std::acos(2*x2 - 1)) - 1.3232496042008038e-8*std::cos(25*std::acos(2*x2 - 1)) + 1.3344667096603976e-8*std::cos(26*std::acos(2*x2 - 1)) - 1.136905698901599e-8*std::cos(27*std::acos(2*x2 - 1)) + 9.3215962954934251e-9*std::cos(28*std::acos(2*x2 - 1)) - 7.8609192978868405e-9*std::cos(29*std::acos(2*x2 - 1)) + 6.9799842727525716e-9*std::cos(30*std::acos(2*x2 - 1)) - 6.4658022609777294e-9*std::cos(31*std::acos(2*x2 - 1)) + 6.1184845934494521e-9*std::cos(32*std::acos(2*x2 - 1)) - 5.8149725227088034e-9*std::cos(33*std::acos(2*x2 - 1)) + 5.5009609448435255e-9*std::cos(34*std::acos(2*x2 - 1)) - 5.1633408925550063e-9*std::cos(35*std::acos(2*x2 - 1)) + 4.8066053217773026e-9*std::cos(36*std::acos(2*x2 - 1)) - 4.439089978043675e-9*std::cos(37*std::acos(2*x2 - 1)) + 4.0672900401123296e-9*std::cos(38*std::acos(2*x2 - 1)) - 3.6947698678390842e-9*std::cos(39*std::acos(2*x2 - 1)) + 3.3229280613673116e-9*std::cos(40*std::acos(2*x2 - 1)) - 2.952095301667594e-9*std::cos(41*std::acos(2*x2 - 1)) + 2.5823712665124437e-9*std::cos(42*std::acos(2*x2 - 1)) - 2.2140699779408913e-9*std::cos(43*std::acos(2*x2 - 1)) + 1.847757120693774e-9*std::cos(44*std::acos(2*x2 - 1)) - 1.4838261042827717e-9*std::cos(45*std::acos(2*x2 - 1)) + 1.1214375157477977e-9*std::cos(46*std::acos(2*x2 - 1)) - 7.5716721864144663e-10*std::cos(47*std::acos(2*x2 - 1)) + 3.8412245556547844e-10*std::cos(48*std::acos(2*x2 - 1)) + 1.0752236232446089e-17*std::cos(49*std::acos(2*x2 - 1)) - 0.0035716180382855955;


         QDA = -0.016493086615981038*x2 + 0.0024749078506912368*std::cos(2*std::acos(2*x2 - 1)) + 0.0028610988856768035*std::cos(3*std::acos(2*x2 - 1)) - 0.0053181685887073822*std::cos(4*std::acos(2*x2 - 1)) + 0.0048992486839919613*std::cos(5*std::acos(2*x2 - 1)) - 0.0030543006044871928*std::cos(6*std::acos(2*x2 - 1)) + 0.0012049278935149274*std::cos(7*std::acos(2*x2 - 1)) - 3.5115081364772042e-5*std::cos(8*std::acos(2*x2 - 1)) - 0.00042790361750880612*std::cos(9*std::acos(2*x2 - 1)) + 0.00044794661409551187*std::cos(10*std::acos(2*x2 - 1)) - 0.00029571583122835361*std::cos(11*std::acos(2*x2 - 1)) + 0.0001363930066605291*std::cos(12*std::acos(2*x2 - 1)) - 3.2193240780892744e-5*std::cos(13*std::acos(2*x2 - 1)) - 1.5196736016867467e-5*std::cos(14*std::acos(2*x2 - 1)) + 2.5962745934999881e-5*std::cos(15*std::acos(2*x2 - 1)) - 2.0366684969529462e-5*std::cos(16*std::acos(2*x2 - 1)) + 1.123796057746961e-5*std::cos(17*std::acos(2*x2 - 1)) - 4.189082421740394e-6*std::cos(18*std::acos(2*x2 - 1)) + 3.5563300549647519e-7*std::cos(19*std::acos(2*x2 - 1)) + 1.0393120133792711e-6*std::cos(20*std::acos(2*x2 - 1)) - 1.1292977776606256e-6*std::cos(21*std::acos(2*x2 - 1)) + 7.5713532536108162e-7*std::cos(22*std::acos(2*x2 - 1)) - 3.6599492144332791e-7*std::cos(23*std::acos(2*x2 - 1)) + 1.0813040298726132e-7*std::cos(24*std::acos(2*x2 - 1)) + 1.4307048162286413e-8*std::cos(25*std::acos(2*x2 - 1)) - 4.9017501999706585e-8*std::cos(26*std::acos(2*x2 - 1)) + 4.286018262123813e-8*std::cos(27*std::acos(2*x2 - 1)) - 2.5641658976336612e-8*std::cos(28*std::acos(2*x2 - 1)) + 1.0969788863219921e-8*std::cos(29*std::acos(2*x2 - 1)) - 2.3402749487583016e-9*std::cos(30*std::acos(2*x2 - 1)) - 1.2702019380814621e-9*std::cos(31*std::acos(2*x2 - 1)) + 1.9790149644569992e-9*std::cos(32*std::acos(2*x2 - 1)) - 1.5011584982395336e-9*std::cos(33*std::acos(2*x2 - 1)) + 8.1078111847791074e-10*std::cos(34*std::acos(2*x2 - 1)) - 2.9679494635319321e-10*std::cos(35*std::acos(2*x2 - 1)) + 2.1814558116636713e-11*std::cos(36*std::acos(2*x2 - 1)) + 7.8232056259831023e-11*std::cos(37*std::acos(2*x2 - 1)) - 8.495590125738693e-11*std::cos(38*std::acos(2*x2 - 1)) + 5.6672156022444616e-11*std::cos(39*std::acos(2*x2 - 1)) - 2.4490396193625746e-11*std::cos(40*std::acos(2*x2 - 1)) + 1.8726160617196939e-12*std::cos(41*std::acos(2*x2 - 1)) + 5.3695921696430341e-12*std::cos(42*std::acos(2*x2 - 1)) + 5.8299153257250137e-12*std::cos(43*std::acos(2*x2 - 1)) - 3.4813884229831616e-11*std::cos(44*std::acos(2*x2 - 1)) + 7.3064121961697187e-11*std::cos(45*std::acos(2*x2 - 1)) - 1.0261916287506907e-10*std::cos(46*std::acos(2*x2 - 1)) + 1.0273971203610565e-10*std::cos(47*std::acos(2*x2 - 1)) - 6.4135092408268518e-11*std::cos(48*std::acos(2*x2 - 1)) + 4.9991757988315085e-17*std::cos(49*std::acos(2*x2 - 1)) + 0.013634481900124368;
        
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