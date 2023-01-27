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
        
        amrex::Real R = 200;
        
        amrex::Real x = x_U_0;
        amrex::Real y = x_U_1;
        amrex::Real z = x_U_2;
        
        
        amrex::Real rdomain = std::sqrt(x*x + y*y + z*z);
        
        amrex::Real rsq = std::pow((std::sqrt(x*x + y*y)-R),2) + z*z;
        amrex::Real r = std::sqrt(rsq);
        
        amrex::Real x2 = r/200;
        
        amrex::Real V, A, B, u;
        
        //amrex::Real x2 = std::sqrt(x_U_0*x_U_0 + x_U_1*x_U_1 + x_U_2*x_U_2)/1000;
        
        if(r < 200.0)
        {
        
        amrex::Real QF = 0.6246372484922188*x2 - 0.24379950891048621*std::cos(2*std::acos(2*x2 - 1)) + 0.15994651481218308*std::cos(3*std::acos(2*x2 - 1)) - 0.086057833075017576*std::cos(4*std::acos(2*x2 - 1)) + 0.035195121426024351*std::cos(5*std::acos(2*x2 - 1)) - 0.0076468361904048661*std::cos(6*std::acos(2*x2 - 1)) - 0.0032691884006957293*std::cos(7*std::acos(2*x2 - 1)) + 0.0052398502045463945*std::cos(8*std::acos(2*x2 - 1)) - 0.0038243195984793525*std::cos(9*std::acos(2*x2 - 1)) + 0.001947154880655954*std::cos(10*std::acos(2*x2 - 1)) - 0.00065040210352139069*std::cos(11*std::acos(2*x2 - 1)) + 2.0337057774086342e-5*std::cos(12*std::acos(2*x2 - 1)) + 0.00017259412824313647*std::cos(13*std::acos(2*x2 - 1)) - 0.00016325922148225161*std::cos(14*std::acos(2*x2 - 1)) + 9.8157964990224639e-5*std::cos(15*std::acos(2*x2 - 1)) - 4.1689601602296761e-5*std::cos(16*std::acos(2*x2 - 1)) + 9.4471127538159338e-6*std::cos(17*std::acos(2*x2 - 1)) + 3.2940905417694344e-6*std::cos(18*std::acos(2*x2 - 1)) - 5.4137927603080394e-6*std::cos(19*std::acos(2*x2 - 1)) + 3.6026798361891987e-6*std::cos(20*std::acos(2*x2 - 1)) - 1.3377620948530622e-6*std::cos(21*std::acos(2*x2 - 1)) - 1.846817664456768e-7*std::cos(22*std::acos(2*x2 - 1)) + 8.9860906989785994e-7*std::cos(23*std::acos(2*x2 - 1)) - 1.0956007049818528e-6*std::cos(24*std::acos(2*x2 - 1)) + 1.060229314345934e-6*std::cos(25*std::acos(2*x2 - 1)) - 9.6554874844806404e-7*std::cos(26*std::acos(2*x2 - 1)) + 8.8525548714787905e-7*std::cos(27*std::acos(2*x2 - 1)) - 8.3548210750780671e-7*std::cos(28*std::acos(2*x2 - 1)) + 8.0923025504164858e-7*std::cos(29*std::acos(2*x2 - 1)) - 7.9539414412840848e-7*std::cos(30*std::acos(2*x2 - 1)) + 7.8593429860690511e-7*std::cos(31*std::acos(2*x2 - 1)) - 7.7680681148314945e-7*std::cos(32*std::acos(2*x2 - 1)) + 7.6669888270444742e-7*std::cos(33*std::acos(2*x2 - 1)) - 7.556059071350818e-7*std::cos(34*std::acos(2*x2 - 1)) + 7.4390961465539696e-7*std::cos(35*std::acos(2*x2 - 1)) - 7.3195839959960818e-7*std::cos(36*std::acos(2*x2 - 1)) + 7.1995687741036041e-7*std::cos(37*std::acos(2*x2 - 1)) - 7.0798885343573148e-7*std::cos(38*std::acos(2*x2 - 1)) + 6.9607056584625182e-7*std::cos(39*std::acos(2*x2 - 1)) - 6.8419202825520387e-7*std::cos(40*std::acos(2*x2 - 1)) + 6.7233851662115707e-7*std::cos(41*std::acos(2*x2 - 1)) - 6.6049753947791667e-7*std::cos(42*std::acos(2*x2 - 1)) + 6.4865893930135532e-7*std::cos(43*std::acos(2*x2 - 1)) - 6.3681400936863031e-7*std::cos(44*std::acos(2*x2 - 1)) + 6.2495676735399575e-7*std::cos(45*std::acos(2*x2 - 1)) - 6.1308841120544244e-7*std::cos(46*std::acos(2*x2 - 1)) + 6.0121849650499049e-7*std::cos(47*std::acos(2*x2 - 1)) - 5.894917391932597e-7*std::cos(48*std::acos(2*x2 - 1)) + 2.8890016208135694e-7*std::cos(49*std::acos(2*x2 - 1)) + 0.51818629699070007;


        amrex::Real QA = 0.83766708793493172*x2 - 0.27736430846593413*std::cos(2*std::acos(2*x2 - 1)) + 0.1266844112694708*std::cos(3*std::acos(2*x2 - 1)) - 0.02210063254468039*std::cos(4*std::acos(2*x2 - 1)) - 0.02445937632172111*std::cos(5*std::acos(2*x2 - 1)) + 0.030953945823981643*std::cos(6*std::acos(2*x2 - 1)) - 0.020855971826156696*std::cos(7*std::acos(2*x2 - 1)) + 0.0093781794112381238*std::cos(8*std::acos(2*x2 - 1)) - 0.0021168265967512097*std::cos(9*std::acos(2*x2 - 1)) - 0.00091810623918994402*std::cos(10*std::acos(2*x2 - 1)) + 0.0014626457210800796*std::cos(11*std::acos(2*x2 - 1)) - 0.0010490189463799558*std::cos(12*std::acos(2*x2 - 1)) + 0.00051575061607879447*std::cos(13*std::acos(2*x2 - 1)) - 0.00015701825220955368*std::cos(14*std::acos(2*x2 - 1)) - 9.8735907418413815e-6*std::cos(15*std::acos(2*x2 - 1)) + 5.5214217680333438e-5*std::cos(16*std::acos(2*x2 - 1)) - 4.7155808067172317e-5*std::cos(17*std::acos(2*x2 - 1)) + 2.6516676034166881e-5*std::cos(18*std::acos(2*x2 - 1)) - 1.0178124572015473e-5*std::cos(19*std::acos(2*x2 - 1)) + 1.4323175082618655e-6*std::cos(20*std::acos(2*x2 - 1)) + 1.6955373663869733e-6*std::cos(21*std::acos(2*x2 - 1)) - 1.9611055969460944e-6*std::cos(22*std::acos(2*x2 - 1)) + 1.281971400139692e-6*std::cos(23*std::acos(2*x2 - 1)) - 5.8725945431895021e-7*std::cos(24*std::acos(2*x2 - 1)) + 1.5621949992980499e-7*std::cos(25*std::acos(2*x2 - 1)) + 3.0870415793082806e-8*std::cos(26*std::acos(2*x2 - 1)) - 7.4002075921980603e-8*std::cos(27*std::acos(2*x2 - 1)) + 5.8148466441226843e-8*std::cos(28*std::acos(2*x2 - 1)) - 3.1060854045056319e-8*std::cos(29*std::acos(2*x2 - 1)) + 1.1116253733484127e-8*std::cos(30*std::acos(2*x2 - 1)) - 9.6368981891692407e-10*std::cos(31*std::acos(2*x2 - 1)) - 2.3943400768999812e-9*std::cos(32*std::acos(2*x2 - 1)) + 2.4594174197201113e-9*std::cos(33*std::acos(2*x2 - 1)) - 1.5308567954952501e-9*std::cos(34*std::acos(2*x2 - 1)) + 6.6771532400422551e-10*std::cos(35*std::acos(2*x2 - 1)) - 1.5594651264771997e-10*std::cos(36*std::acos(2*x2 - 1)) - 5.601652822151415e-11*std::cos(37*std::acos(2*x2 - 1)) + 9.8795898545445127e-11*std::cos(38*std::acos(2*x2 - 1)) - 7.5200200233549452e-11*std::cos(39*std::acos(2*x2 - 1)) + 3.9221591252580623e-11*std::cos(40*std::acos(2*x2 - 1)) - 8.8815309273737597e-12*std::cos(41*std::acos(2*x2 - 1)) - 1.4892965333190844e-11*std::cos(42*std::acos(2*x2 - 1)) + 3.3803636972917772e-11*std::cos(43*std::acos(2*x2 - 1)) - 4.4671339544111177e-11*std::cos(44*std::acos(2*x2 - 1)) + 3.9547348035240404e-11*std::cos(45*std::acos(2*x2 - 1)) - 1.2936988286194051e-11*std::cos(46*std::acos(2*x2 - 1)) - 2.8795500706113586e-11*std::cos(47*std::acos(2*x2 - 1)) + 6.5579219138411915e-11*std::cos(48*std::acos(2*x2 - 1)) - 3.9322252051721129e-11*std::cos(49*std::acos(2*x2 - 1)) + 0.34234270421444901;


        amrex::Real QU = -0.03368113546039312*x2 + 0.0061944504724723006*std::cos(2*std::acos(2*x2 - 1)) - 0.001532799734729521*std::cos(3*std::acos(2*x2 - 1)) - 0.00013715470617930636*std::cos(4*std::acos(2*x2 - 1)) + 0.00053992176738089045*std::cos(5*std::acos(2*x2 - 1)) - 0.00048409619135376719*std::cos(6*std::acos(2*x2 - 1)) + 0.00031102310572364935*std::cos(7*std::acos(2*x2 - 1)) - 0.00015384794051652472*std::cos(8*std::acos(2*x2 - 1)) + 5.1211509726755624e-5*std::cos(9*std::acos(2*x2 - 1)) + 3.4442552164378736e-8*std::cos(10*std::acos(2*x2 - 1)) - 1.7026253973486968e-5*std::cos(11*std::acos(2*x2 - 1)) + 1.6763747160852535e-5*std::cos(12*std::acos(2*x2 - 1)) - 1.0813551366486713e-5*std::cos(13*std::acos(2*x2 - 1)) + 5.0389888693876071e-6*std::cos(14*std::acos(2*x2 - 1)) - 1.3544096146829475e-6*std::cos(15*std::acos(2*x2 - 1)) - 3.329310401443329e-7*std::cos(16*std::acos(2*x2 - 1)) + 7.615049758025313e-7*std::cos(17*std::acos(2*x2 - 1)) - 6.2480961967057871e-7*std::cos(18*std::acos(2*x2 - 1)) + 3.5379092416969501e-7*std::cos(19*std::acos(2*x2 - 1)) - 1.3773795006852003e-7*std::cos(20*std::acos(2*x2 - 1)) + 1.8181743042126315e-8*std::cos(21*std::acos(2*x2 - 1)) + 2.6645683056620434e-8*std::cos(22*std::acos(2*x2 - 1)) - 3.1109104372452044e-8*std::cos(23*std::acos(2*x2 - 1)) + 2.1053445530385125e-8*std::cos(24*std::acos(2*x2 - 1)) - 1.0086426963701467e-8*std::cos(25*std::acos(2*x2 - 1)) + 2.881237645472809e-9*std::cos(26*std::acos(2*x2 - 1)) + 4.6756461123328953e-10*std::cos(27*std::acos(2*x2 - 1)) - 1.3578810110830388e-9*std::cos(28*std::acos(2*x2 - 1)) + 1.1432459295749061e-9*std::cos(29*std::acos(2*x2 - 1)) - 6.6078602318253708e-10*std::cos(30*std::acos(2*x2 - 1)) + 2.7349463443201573e-10*std::cos(31*std::acos(2*x2 - 1)) - 6.0019354043629995e-11*std::cos(32*std::acos(2*x2 - 1)) - 1.9564258660200917e-11*std::cos(33*std::acos(2*x2 - 1)) + 2.7423647710375045e-11*std::cos(34*std::acos(2*x2 - 1)) - 9.5679867878537187e-12*std::cos(35*std::acos(2*x2 - 1)) - 9.8580039914831062e-12*std::cos(36*std::acos(2*x2 - 1)) + 2.1776579607605749e-11*std::cos(37*std::acos(2*x2 - 1)) - 2.4075605177070365e-11*std::cos(38*std::acos(2*x2 - 1)) + 1.6608479589525537e-11*std::cos(39*std::acos(2*x2 - 1)) - 1.0716504658432912e-13*std::cos(40*std::acos(2*x2 - 1)) - 2.1729073017675932e-11*std::cos(41*std::acos(2*x2 - 1)) + 3.8693604457232443e-11*std::cos(42*std::acos(2*x2 - 1)) - 3.2052747282337136e-11*std::cos(43*std::acos(2*x2 - 1)) - 2.1961768969056699e-11*std::cos(44*std::acos(2*x2 - 1)) + 1.4041134189622853e-10*std::cos(45*std::acos(2*x2 - 1)) - 3.1682942723529622e-10*std::cos(46*std::acos(2*x2 - 1)) + 5.1179245464550164e-10*std::cos(47*std::acos(2*x2 - 1)) - 6.6595870786439573e-10*std::cos(48*std::acos(2*x2 - 1)) + 3.586452579757089e-10*std::cos(49*std::acos(2*x2 - 1)) + 0.028899737031107469;


        amrex::Real QDF = -0.032769941177544157*x2 + 0.011909125489846862*std::cos(2*std::acos(2*x2 - 1)) - 0.0066329414722070068*std::cos(3*std::acos(2*x2 - 1)) + 0.0023122866169214907*std::cos(4*std::acos(2*x2 - 1)) + 0.00025171959709980744*std::cos(5*std::acos(2*x2 - 1)) - 0.0012072431233295965*std::cos(6*std::acos(2*x2 - 1)) + 0.0011693445280730414*std::cos(7*std::acos(2*x2 - 1)) - 0.00074955445878884083*std::cos(8*std::acos(2*x2 - 1)) + 0.00033096430344445652*std::cos(9*std::acos(2*x2 - 1)) - 6.1173489157778795e-5*std::cos(10*std::acos(2*x2 - 1)) - 5.8468619851397793e-5*std::cos(11*std::acos(2*x2 - 1)) + 8.1915689062832351e-5*std::cos(12*std::acos(2*x2 - 1)) - 6.3349538121746307e-5*std::cos(13*std::acos(2*x2 - 1)) + 3.7040991346174948e-5*std::cos(14*std::acos(2*x2 - 1)) - 1.7636727542662859e-5*std::cos(15*std::acos(2*x2 - 1)) + 7.5934546113664019e-6*std::cos(16*std::acos(2*x2 - 1)) - 4.2959883262047035e-6*std::cos(17*std::acos(2*x2 - 1)) + 4.3814202148903163e-6*std::cos(18*std::acos(2*x2 - 1)) - 5.4818668506258547e-6*std::cos(19*std::acos(2*x2 - 1)) + 6.4386717501064637e-6*std::cos(20*std::acos(2*x2 - 1)) - 6.9229459904862905e-6*std::cos(21*std::acos(2*x2 - 1)) + 7.0005346392472914e-6*std::cos(22*std::acos(2*x2 - 1)) - 6.8416856068923777e-6*std::cos(23*std::acos(2*x2 - 1)) + 6.5871724001730209e-6*std::cos(24*std::acos(2*x2 - 1)) - 6.3157946392287932e-6*std::cos(25*std::acos(2*x2 - 1)) + 6.0570550925137214e-6*std::cos(26*std::acos(2*x2 - 1)) - 5.8137067795039838e-6*std::cos(27*std::acos(2*x2 - 1)) + 5.5790147391337995e-6*std::cos(28*std::acos(2*x2 - 1)) - 5.3458344600105165e-6*std::cos(29*std::acos(2*x2 - 1)) + 5.1096588445157928e-6*std::cos(30*std::acos(2*x2 - 1)) - 4.8685955874266637e-6*std::cos(31*std::acos(2*x2 - 1)) + 4.622377142882219e-6*std::cos(32*std::acos(2*x2 - 1)) - 4.3714367421902926e-6*std::cos(33*std::acos(2*x2 - 1)) + 4.116353350228895e-6*std::cos(34*std::acos(2*x2 - 1)) - 3.8576221562672521e-6*std::cos(35*std::acos(2*x2 - 1)) + 3.5956140162457081e-6*std::cos(36*std::acos(2*x2 - 1)) - 3.33060947346692e-6*std::cos(37*std::acos(2*x2 - 1)) + 3.0628432630968142e-6*std::cos(38*std::acos(2*x2 - 1)) - 2.7925352545261829e-6*std::cos(39*std::acos(2*x2 - 1)) + 2.5199055070002257e-6*std::cos(40*std::acos(2*x2 - 1)) - 2.2451788951595411e-6*std::cos(41*std::acos(2*x2 - 1)) + 1.9685851668215511e-6*std::cos(42*std::acos(2*x2 - 1)) - 1.6903581879199288e-6*std::cos(43*std::acos(2*x2 - 1)) + 1.4107356897864229e-6*std::cos(44*std::acos(2*x2 - 1)) - 1.1299590576081186e-6*std::cos(45*std::acos(2*x2 - 1)) + 8.4827178684777625e-7*std::cos(46*std::acos(2*x2 - 1)) - 5.6591489906394156e-7*std::cos(47*std::acos(2*x2 - 1)) + 2.8312357446030206e-7*std::cos(48*std::acos(2*x2 - 1)) + 6.9853230553814623e-17*std::cos(49*std::acos(2*x2 - 1)) + 0.025462735192165901;


        amrex::Real QDA = -0.016200673460423136*x2 + 0.00228861147416075*std::cos(2*std::acos(2*x2 - 1)) + 0.0029942910815648886*std::cos(3*std::acos(2*x2 - 1)) - 0.0053124912075209352*std::cos(4*std::acos(2*x2 - 1)) + 0.0047623505254366268*std::cos(5*std::acos(2*x2 - 1)) - 0.0028665413455994356*std::cos(6*std::acos(2*x2 - 1)) + 0.0010478584540984512*std::cos(7*std::acos(2*x2 - 1)) + 5.3309309315775454e-5*std::cos(8*std::acos(2*x2 - 1)) - 0.00045265775428071701*std::cos(9*std::acos(2*x2 - 1)) + 0.00043434000188440666*std::cos(10*std::acos(2*x2 - 1)) - 0.00026903558833193307*std::cos(11*std::acos(2*x2 - 1)) + 0.00011255633432848434*std::cos(12*std::acos(2*x2 - 1)) - 1.7269782371703762e-5*std::cos(13*std::acos(2*x2 - 1)) - 2.1539496331137988e-5*std::cos(14*std::acos(2*x2 - 1)) + 2.6695548073611686e-5*std::cos(15*std::acos(2*x2 - 1)) - 1.8577404298188185e-5*std::cos(16*std::acos(2*x2 - 1)) + 9.0269100727267289e-6*std::cos(17*std::acos(2*x2 - 1)) - 2.5443493901095262e-6*std::cos(18*std::acos(2*x2 - 1)) - 5.1914102979345557e-7*std::cos(19*std::acos(2*x2 - 1)) + 1.3233572858197992e-6*std::cos(20*std::acos(2*x2 - 1)) - 1.0920708977732423e-6*std::cos(21*std::acos(2*x2 - 1)) + 6.1122803129302428e-7*std::cos(22*std::acos(2*x2 - 1)) - 2.2918012058967663e-7*std::cos(23*std::acos(2*x2 - 1)) + 2.15182385914147e-8*std::cos(24*std::acos(2*x2 - 1)) + 5.270582676546149e-8*std::cos(25*std::acos(2*x2 - 1)) - 5.6591901828755266e-8*std::cos(26*std::acos(2*x2 - 1)) + 3.6653130360110266e-8*std::cos(27*std::acos(2*x2 - 1)) - 1.6630581124470007e-8*std::cos(28*std::acos(2*x2 - 1)) + 4.0898262789215038e-9*std::cos(29*std::acos(2*x2 - 1)) + 1.3848043454288893e-9*std::cos(30*std::acos(2*x2 - 1)) - 2.5799593668101729e-9*std::cos(31*std::acos(2*x2 - 1)) + 1.9822949651479151e-9*std::cos(32*std::acos(2*x2 - 1)) - 1.0475739489498342e-9*std::cos(33*std::acos(2*x2 - 1)) + 3.5907139597012134e-10*std::cos(34*std::acos(2*x2 - 1)) - 6.5861421044858708e-12*std::cos(35*std::acos(2*x2 - 1)) - 1.0833168367671382e-10*std::cos(36*std::acos(2*x2 - 1)) + 1.056959093932795e-10*std::cos(37*std::acos(2*x2 - 1)) - 6.6879223034016513e-11*std::cos(38*std::acos(2*x2 - 1)) + 3.0610620646041326e-11*std::cos(39*std::acos(2*x2 - 1)) - 8.2228171718400669e-12*std::cos(40*std::acos(2*x2 - 1)) - 7.6681190902805512e-13*std::cos(41*std::acos(2*x2 - 1)) - 9.3990667779948038e-13*std::cos(42*std::acos(2*x2 - 1)) + 1.1743339558218337e-11*std::cos(43*std::acos(2*x2 - 1)) - 3.001117274328166e-11*std::cos(44*std::acos(2*x2 - 1)) + 5.1054380560903951e-11*std::cos(45*std::acos(2*x2 - 1)) - 6.5603985409193975e-11*std::cos(46*std::acos(2*x2 - 1)) + 6.295649088578288e-11*std::cos(47*std::acos(2*x2 - 1)) - 3.8536011535257519e-11*std::cos(48*std::acos(2*x2 - 1)) + 4.6391466779887221e-17*std::cos(49*std::acos(2*x2 - 1)) + 0.013432998848748524;
        
        amrex::Real lambda = 2.0;
        amrex::Real eta = 0.05;
        amrex::Real e = 1.0;
        
        V = 0.25*lambda*std::pow(eta,4)*std::pow((QF*QF-1),2);
        A = eta*eta/(r*r)*(QF*QF*(1-QA)*(1-QA) + r*r*QDF*QDF);
        B = QDA*QDA/(2*e*e*r*r);
        u = QU;
        }
        
        else
        {
        V = 0;
        A = 0;
        B = 0;
        u = 0;
        }
        
        state_fab(i, j, k, Idx::V) = V;
        state_fab(i, j, k, Idx::A) = A;
        state_fab(i, j, k, Idx::B) = B;
        
        
        //amrex::Real rho = V+A+B;
        
        
        state_fab(i, j, k, Idx::u) = 0; //u;
        state_fab(i, j, k, Idx::chi) = 1; //std::pow(1+u,-4);
        
}
#endif