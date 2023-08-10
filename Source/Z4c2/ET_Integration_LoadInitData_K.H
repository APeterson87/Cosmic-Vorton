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
        amrex::Real e = std::sqrt(lambda/2.0);
        
        amrex::Real rdomain = std::sqrt(x*x + y*y + z*z);
        
        amrex::Real rsq = std::pow((std::sqrt(x*x + y*y)-R),2) + z*z;
        amrex::Real r = std::sqrt(rsq);
        
        amrex::Real x2 = r/500;
        
        
        amrex::Real QF, QA;
        
        if(r < 500.0)
        {
            QF = 0.4027642425104318*x2 - 0.18389804471411586*std::cos(2*std::acos(2*x2 - 1)) + 0.15785299937437691*std::cos(3*std::acos(2*x2 - 1)) - 0.12705591130856209*std::cos(4*std::acos(2*x2 - 1)) + 0.095476923821582724*std::cos(5*std::acos(2*x2 - 1)) - 0.066454401105316849*std::cos(6*std::acos(2*x2 - 1)) + 0.042203947720811624*std::cos(7*std::acos(2*x2 - 1)) - 0.023695995294027064*std::cos(8*std::acos(2*x2 - 1)) + 0.010837852846319043*std::cos(9*std::acos(2*x2 - 1)) - 0.002823817816695386*std::cos(10*std::acos(2*x2 - 1)) - 0.0014872210869433123*std::cos(11*std::acos(2*x2 - 1)) + 0.0032650768876229929*std::cos(12*std::acos(2*x2 - 1)) - 0.00351448500401297*std::cos(13*std::acos(2*x2 - 1)) + 0.0029904875622377026*std::cos(14*std::acos(2*x2 - 1)) - 0.0021942044300177956*std::cos(15*std::acos(2*x2 - 1)) + 0.0014142663155602229*std::cos(16*std::acos(2*x2 - 1)) - 0.00078499575203161941*std::cos(17*std::acos(2*x2 - 1)) + 0.00034255646011249578*std::cos(18*std::acos(2*x2 - 1)) - 6.9827385533750494e-5*std::cos(19*std::acos(2*x2 - 1)) - 7.2439677694600621e-5*std::cos(20*std::acos(2*x2 - 1)) + 0.00012676980842400189*std::cos(21*std::acos(2*x2 - 1)) - 0.00012943688543993505*std::cos(22*std::acos(2*x2 - 1)) + 0.00010710284965785057*std::cos(23*std::acos(2*x2 - 1)) - 7.6990340906982202e-5*std::cos(24*std::acos(2*x2 - 1)) + 4.869405016596183e-5*std::cos(25*std::acos(2*x2 - 1)) - 2.6435364062852781e-5*std::cos(26*std::acos(2*x2 - 1)) + 1.1114226556524845e-5*std::cos(27*std::acos(2*x2 - 1)) - 1.8896726518488616e-6*std::cos(28*std::acos(2*x2 - 1)) - 2.7520118106080099e-6*std::cos(29*std::acos(2*x2 - 1)) + 4.3673727191367168e-6*std::cos(30*std::acos(2*x2 - 1)) - 4.248918325301318e-6*std::cos(31*std::acos(2*x2 - 1)) + 3.3276444457788767e-6*std::cos(32*std::acos(2*x2 - 1)) - 2.1927722615843609e-6*std::cos(33*std::acos(2*x2 - 1)) + 1.1635335836236715e-6*std::cos(34*std::acos(2*x2 - 1)) - 3.725628428142419e-7*std::cos(35*std::acos(2*x2 - 1)) - 1.604239206302871e-7*std::cos(36*std::acos(2*x2 - 1)) + 4.7301246562832244e-7*std::cos(37*std::acos(2*x2 - 1)) - 6.2320288779678379e-7*std::cos(38*std::acos(2*x2 - 1)) + 6.6814971249294641e-7*std::cos(39*std::acos(2*x2 - 1)) - 6.5420500528282077e-7*std::cos(40*std::acos(2*x2 - 1)) + 6.1414874707489686e-7*std::cos(41*std::acos(2*x2 - 1)) - 5.6834075525154237e-7*std::cos(42*std::acos(2*x2 - 1)) + 5.2752763162952632e-7*std::cos(43*std::acos(2*x2 - 1)) - 4.9592603294963677e-7*std::cos(44*std::acos(2*x2 - 1)) + 4.7387666495153635e-7*std::cos(45*std::acos(2*x2 - 1)) - 4.5981173599739245e-7*std::cos(46*std::acos(2*x2 - 1)) + 4.5153822701704538e-7*std::cos(47*std::acos(2*x2 - 1)) - 4.4696240472652526e-7*std::cos(48*std::acos(2*x2 - 1)) + 4.4441700539597784e-7*std::cos(49*std::acos(2*x2 - 1)) - 4.4273988427691313e-7*std::cos(50*std::acos(2*x2 - 1)) + 4.4121997504892319e-7*std::cos(51*std::acos(2*x2 - 1)) - 4.3948993997192098e-7*std::cos(52*std::acos(2*x2 - 1)) + 4.3741280377103398e-7*std::cos(53*std::acos(2*x2 - 1)) - 4.3498620652844544e-7*std::cos(54*std::acos(2*x2 - 1)) + 4.3227246420037002e-7*std::cos(55*std::acos(2*x2 - 1)) - 4.2935355873234549e-7*std::cos(56*std::acos(2*x2 - 1)) + 4.2630604836328689e-7*std::cos(57*std::acos(2*x2 - 1)) - 4.2319007094275551e-7*std::cos(58*std::acos(2*x2 - 1)) + 4.2004712984961001e-7*std::cos(59*std::acos(2*x2 - 1)) - 4.1690234687564998e-7*std::cos(60*std::acos(2*x2 - 1)) + 4.1376846885397833e-7*std::cos(61*std::acos(2*x2 - 1)) - 4.1065005622863138e-7*std::cos(62*std::acos(2*x2 - 1)) + 4.0754692697854422e-7*std::cos(63*std::acos(2*x2 - 1)) - 4.0445664505409096e-7*std::cos(64*std::acos(2*x2 - 1)) + 4.013761398611454e-7*std::cos(65*std::acos(2*x2 - 1)) - 3.9830260255106775e-7*std::cos(66*std::acos(2*x2 - 1)) + 3.9523382724081235e-7*std::cos(67*std::acos(2*x2 - 1)) - 3.9216827570076607e-7*std::cos(68*std::acos(2*x2 - 1)) + 3.8910502678755898e-7*std::cos(69*std::acos(2*x2 - 1)) - 3.8604361168523282e-7*std::cos(70*std::acos(2*x2 - 1)) + 3.8298384341141356e-7*std::cos(71*std::acos(2*x2 - 1)) - 3.7992569851090507e-7*std::cos(72*std::acos(2*x2 - 1)) + 3.7686924220063645e-7*std::cos(73*std::acos(2*x2 - 1)) - 3.7381456945324854e-7*std::cos(74*std::acos(2*x2 - 1)) + 3.7076175154505953e-7*std::cos(75*std::acos(2*x2 - 1)) - 3.6771084349282446e-7*std::cos(76*std::acos(2*x2 - 1)) + 3.6466188321933324e-7*std::cos(77*std::acos(2*x2 - 1)) - 3.6161488620699289e-7*std::cos(78*std::acos(2*x2 - 1)) + 3.585698741866844e-7*std::cos(79*std::acos(2*x2 - 1)) - 3.5552686371113917e-7*std::cos(80*std::acos(2*x2 - 1)) + 3.5248585450280145e-7*std::cos(81*std::acos(2*x2 - 1)) - 3.4944686611027009e-7*std::cos(82*std::acos(2*x2 - 1)) + 3.4640992565074247e-7*std::cos(83*std::acos(2*x2 - 1)) - 3.4337501315408192e-7*std::cos(84*std::acos(2*x2 - 1)) + 3.403420777998295e-7*std::cos(85*std::acos(2*x2 - 1)) - 3.3731102507157662e-7*std::cos(86*std::acos(2*x2 - 1)) + 3.3428164418480844e-7*std::cos(87*std::acos(2*x2 - 1)) - 3.3125361649163909e-7*std::cos(88*std::acos(2*x2 - 1)) + 3.2822654731819167e-7*std::cos(89*std::acos(2*x2 - 1)) - 3.2519996431140674e-7*std::cos(90*std::acos(2*x2 - 1)) + 3.2217339213276575e-7*std::cos(91*std::acos(2*x2 - 1)) - 3.1914651246744974e-7*std::cos(92*std::acos(2*x2 - 1)) + 3.1611928258921951e-7*std::cos(93*std::acos(2*x2 - 1)) - 3.1309200383689006e-7*std::cos(94*std::acos(2*x2 - 1)) + 3.1006537510799848e-7*std::cos(95*std::acos(2*x2 - 1)) - 3.0704084429296341e-7*std::cos(96*std::acos(2*x2 - 1)) + 3.0401916045232069e-7*std::cos(97*std::acos(2*x2 - 1)) - 3.0104901199527745e-7*std::cos(98*std::acos(2*x2 - 1)) + 1.4904154979188816e-7*std::cos(99*std::acos(2*x2 - 1)) + 0.6948447478116353;


        QA = 0.55276503846435909*x2 - 0.23936818884352123*std::cos(2*std::acos(2*x2 - 1)) + 0.18662182001817759*std::cos(3*std::acos(2*x2 - 1)) - 0.12846424571256707*std::cos(4*std::acos(2*x2 - 1)) + 0.074552212778171589*std::cos(5*std::acos(2*x2 - 1)) - 0.031566725668219449*std::cos(6*std::acos(2*x2 - 1)) + 0.0022962348727725099*std::cos(7*std::acos(2*x2 - 1)) + 0.013917876902396703*std::cos(8*std::acos(2*x2 - 1)) - 0.019914492249186726*std::cos(9*std::acos(2*x2 - 1)) + 0.019295304251206263*std::cos(10*std::acos(2*x2 - 1)) - 0.015343330578890839*std::cos(11*std::acos(2*x2 - 1)) + 0.010461431282100952*std::cos(12*std::acos(2*x2 - 1)) - 0.006068246695419487*std::cos(13*std::acos(2*x2 - 1)) + 0.0027734527376973035*std::cos(14*std::acos(2*x2 - 1)) - 0.00065389755578986776*std::cos(15*std::acos(2*x2 - 1)) - 0.00048914074586221302*std::cos(16*std::acos(2*x2 - 1)) + 0.00094483857333296263*std::cos(17*std::acos(2*x2 - 1)) - 0.00098576512449042274*std::cos(18*std::acos(2*x2 - 1)) + 0.00082158787425834316*std::cos(19*std::acos(2*x2 - 1)) - 0.00059082213384450792*std::cos(20*std::acos(2*x2 - 1)) + 0.00037157081218074615*std::cos(21*std::acos(2*x2 - 1)) - 0.00019856800733484246*std::cos(22*std::acos(2*x2 - 1)) + 7.9578274839640917e-5*std::cos(23*std::acos(2*x2 - 1)) - 8.2216888077826095e-6*std::cos(24*std::acos(2*x2 - 1)) - 2.7349960891508152e-5*std::cos(25*std::acos(2*x2 - 1)) + 3.936060005058431e-5*std::cos(26*std::acos(2*x2 - 1)) - 3.7947861057990515e-5*std::cos(27*std::acos(2*x2 - 1)) + 3.039151172942213e-5*std::cos(28*std::acos(2*x2 - 1)) - 2.1278087551670349e-5*std::cos(29*std::acos(2*x2 - 1)) + 1.3075166567800764e-5*std::cos(30*std::acos(2*x2 - 1)) - 6.7926111380362669e-6*std::cos(31*std::acos(2*x2 - 1)) + 2.563100425916412e-6*std::cos(32*std::acos(2*x2 - 1)) - 7.7886171151453354e-8*std::cos(33*std::acos(2*x2 - 1)) - 1.1258144354906602e-6*std::cos(34*std::acos(2*x2 - 1)) + 1.5002500691387727e-6*std::cos(35*std::acos(2*x2 - 1)) - 1.4108203638153738e-6*std::cos(36*std::acos(2*x2 - 1)) + 1.1155264080867089e-6*std::cos(37*std::acos(2*x2 - 1)) - 7.74629689165654e-7*std::cos(38*std::acos(2*x2 - 1)) + 4.7290932321357448e-7*std::cos(39*std::acos(2*x2 - 1)) - 2.438783735358302e-7*std::cos(40*std::acos(2*x2 - 1)) + 9.055613165408205e-8*std::cos(41*std::acos(2*x2 - 1)) - 8.5459673539178826e-10*std::cos(42*std::acos(2*x2 - 1)) - 4.2369955348051103e-8*std::cos(43*std::acos(2*x2 - 1)) + 5.561387308485044e-8*std::cos(44*std::acos(2*x2 - 1)) - 5.2119078651838291e-8*std::cos(45*std::acos(2*x2 - 1)) + 4.1209347299148158e-8*std::cos(46*std::acos(2*x2 - 1)) - 2.8672564131998479e-8*std::cos(47*std::acos(2*x2 - 1)) + 1.7573916044621996e-8*std::cos(48*std::acos(2*x2 - 1)) - 9.1295373882277486e-9*std::cos(49*std::acos(2*x2 - 1)) + 3.4534277295467097e-9*std::cos(50*std::acos(2*x2 - 1)) - 1.1082883222268691e-10*std::cos(51*std::acos(2*x2 - 1)) - 1.5194168319354606e-9*std::cos(52*std::acos(2*x2 - 1)) + 2.0381040848060916e-9*std::cos(53*std::acos(2*x2 - 1)) - 1.9301244196662237e-9*std::cos(54*std::acos(2*x2 - 1)) + 1.5391657618574195e-9*std::cos(55*std::acos(2*x2 - 1)) - 1.0804797560132195e-9*std::cos(56*std::acos(2*x2 - 1)) + 6.6964058777208013e-10*std::cos(57*std::acos(2*x2 - 1)) - 3.5392306620141145e-10*std::cos(58*std::acos(2*x2 - 1)) + 1.3943033012886463e-10*std::cos(59*std::acos(2*x2 - 1)) - 1.1363103166739386e-11*std::cos(60*std::acos(2*x2 - 1)) - 5.2556100096912139e-11*std::cos(61*std::acos(2*x2 - 1)) + 7.4294580904021856e-11*std::cos(62*std::acos(2*x2 - 1)) - 7.1847101185129603e-11*std::cos(63*std::acos(2*x2 - 1)) + 5.8144697390183353e-11*std::cos(64*std::acos(2*x2 - 1)) - 4.1388821189758396e-11*std::cos(65*std::acos(2*x2 - 1)) + 2.606783192082851e-11*std::cos(66*std::acos(2*x2 - 1)) - 1.4112860186044074e-11*std::cos(67*std::acos(2*x2 - 1)) + 5.8657158991115566e-12*std::cos(68*std::acos(2*x2 - 1)) - 8.3035488207539032e-13*std::cos(69*std::acos(2*x2 - 1)) - 1.7962905468627e-12*std::cos(70*std::acos(2*x2 - 1)) + 2.8366042154059912e-12*std::cos(71*std::acos(2*x2 - 1)) - 2.9683512592360017e-12*std::cos(72*std::acos(2*x2 - 1)) + 2.6695798349685163e-12*std::cos(73*std::acos(2*x2 - 1)) - 2.2462032234216167e-12*std::cos(74*std::acos(2*x2 - 1)) + 1.8362152076623062e-12*std::cos(75*std::acos(2*x2 - 1)) - 1.455517997794864e-12*std::cos(76*std::acos(2*x2 - 1)) + 1.0487721802121541e-12*std::cos(77*std::acos(2*x2 - 1)) - 5.2201298839094079e-13*std::cos(78*std::acos(2*x2 - 1)) - 2.1256434112881806e-13*std::cos(79*std::acos(2*x2 - 1)) + 1.2058444520679501e-12*std::cos(80*std::acos(2*x2 - 1)) - 2.4343999038833886e-12*std::cos(81*std::acos(2*x2 - 1)) + 3.7616246922889118e-12*std::cos(82*std::acos(2*x2 - 1)) - 4.9228000148504591e-12*std::cos(83*std::acos(2*x2 - 1)) + 5.5308726348846093e-12*std::cos(84*std::acos(2*x2 - 1)) - 5.1208603329966351e-12*std::cos(85*std::acos(2*x2 - 1)) + 3.2572208819026116e-12*std::cos(86*std::acos(2*x2 - 1)) + 3.3430723467287038e-13*std::cos(87*std::acos(2*x2 - 1)) - 5.5996058484497624e-12*std::cos(88*std::acos(2*x2 - 1)) + 1.1977717429001444e-11*std::cos(89*std::acos(2*x2 - 1)) - 1.8343038893364749e-11*std::cos(90*std::acos(2*x2 - 1)) + 2.310922460335707e-11*std::cos(91*std::acos(2*x2 - 1)) - 2.4527354106074206e-11*std::cos(92*std::acos(2*x2 - 1)) + 2.1207103781395453e-11*std::cos(93*std::acos(2*x2 - 1)) - 1.2703683591186454e-11*std::cos(94*std::acos(2*x2 - 1)) - 5.2582938003808977e-14*std::cos(95*std::acos(2*x2 - 1)) + 1.4714741577992463e-11*std::cos(96*std::acos(2*x2 - 1)) - 2.8015218617172266e-11*std::cos(97*std::acos(2*x2 - 1)) + 3.6693972513268669e-11*std::cos(98*std::acos(2*x2 - 1)) - 1.9650103592894208e-11*std::cos(99*std::acos(2*x2 - 1)) + 0.57875914602311962;
        }
        
        else
        {
            QF = 1.0;
            QA = 1.0;
        }
        
        //QA = 0;
        /*Black Hole binary puncture data from initial data solver */
        
        
        

        amrex::Real psi = 1 + initial_data(i, j, k, InitIdx::u);
        
        amrex::Real alphapsi = initial_data(i, j, k, InitIdx::alphapsi);
        
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

        state_fab(i, j, k, Idx::alpha) = 1; //std::pow(psi,-2); //salphapsi/psi;

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