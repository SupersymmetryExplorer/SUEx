#include <stdio.h>
#include <math.h>
#include "variables.h"
#include "loop.h"
#include "complex.h"

// ============================================================================================================
// Self energy of W and Z particles
// ============================================================================================================

//
// Precision corrections in the minimal supersymmetric standard model
//		By D.M. Pierce, Jonathan A. Bagger, K.T. Matchev and Ren-jie Zhang
//		arXiv:0709.1075v1  [hep-ph]       ....   (All the equations are given in the appendix of this paper)
// ------------------------------------------------------------------------------------------------------------

//
// Calculates essentially the pizz(q),piww(q) at q=mz and q=0 for calculating s^2_w, g1, g2 (mz) in Drbar scheme
// Where pizz, and piww are self energies of the Z and W bosons
// ------------------------------------------------------------------------------------------------------------

void su_pixx(double sw2, double g, double g1, double tbeta, double *pizz, double *piww, double *piww0, double rmtop)
{

    // May be something wrong in the variable passing.. Things are not matching

    double mb = mbpole;
    double mt = mtpole;

    double sq2, cw2, sw, cw, cwm2;
    double mup, mdo, me, mmu;
    double mcq, eps, eps0, ms, ct, st, cb, sb;
    double cta, sta, cbeta2, cbet, sbet;
    double c2b, sal, cal, s2al, mtsave, s2a;

    double mc1 = dmc1, mc2 = dmc2, mn1 = dmn1, mn2 = dmn2, mn3 = dmn3, mn4 = dmn4, gluino = mgluino;

    //
    // basic parameters and definitions used:
    //----------------------------------------------------------------------------------------------------------------------------

    sq2 = sqrt(2.0);
    cw2 = 1.0 - sw2;
    sw = sqrt(sw2);
    cw = sqrt(cw2);
    cwm2 = 1.0 / cw2;

    //
    // defining other parameters for the higgs mass calculation
    //----------------------------------------------------------------------------------------------------------------------------
    double b;
    b = atan(tbeta);
    beta = b;
    mup = 1.e-2;
    mdo = 1.e-2;
    me = 0.5e-3;
    mmu = 0.106;
    ms = 0.190;
    mcq = 1.42;
    eps = 1.0e-2;
    eps0 = pow(eps, 2);

    gmn[1] = fabs(dxmn[1]); // Probabily the Nutralino and the Chargino Masses
    gmn[2] = fabs(dxmn[2]);
    gmn[3] = fabs(dxmn[3]);
    gmn[4] = fabs(dxmn[4]);
    gmc[1] = mc1;
    gmc[2] = mc2;

    ct = cos(thet);
    st = sin(thet);
    cb = cos(theb);
    sb = sin(theb);
    cta = cos(thel);
    sta = sin(thel);

    cbeta2 = 1.0 / (1.0 + pow(tbeta, 2));
    cbet = sqrt(cbeta2);
    sbet = sqrt(1.0 - cbeta2);
    c2b = 2 * cbeta2 - 1.0;

    sal = sin(alfa);
    cal = cos(alfa);
    s2a = 2 * sal * cal;
    s2al = s2a;

    mtsave = mt;

    if (rmtop != 0.0)
        mt = rmtop;

    //
    // Z boson self-energy at pow(q,2)=pow(mz,2)
    //----------------------------------------------------------------------------------------------------------------------------

    double qsz, pizzf, pizzb, pizzh0, pizzhs, pizzsu, pizzsd;
    double pizzsl, pizzc, pizzn, pizzsusy, pizzs, pizzf0;
    double pizzsm, pizzb0, pizzh00, pizzhs0, pizzsu0;

    qsz = pow(mz, 2);

    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Contribution from standard model particles (3rd from last of Eq. D4).
    // -------------------------------------------------------------------------

    pizzf = 3 * (pow((.50 - 2 * sw2 / 3), 2) + pow((2 * sw2 / 3), 2)) * (creal(su_bh(qsz, mt, mt)) + creal(su_bh(qsz, mcq, mcq)) + creal(su_bh(qsz, mup, mup)));
    pizzf += +3 * (pow((-.50 + sw2 / 3), 2) + pow((-sw2 / 3), 2)) * (creal(su_bh(qsz, mb, mb)) + creal(su_bh(qsz, ms, ms)) + creal(su_bh(qsz, mdo, mdo)));
    pizzf += +(pow((-.50 + sw2), 2) + pow((-sw2), 2)) * (creal(su_bh(qsz, me, me)) + creal(su_bh(qsz, mmu, mmu)) + creal(su_bh(qsz, mtau, mtau)));
    pizzf += +pow(.50, 2) * 3 * creal(su_bh(qsz, eps, eps));
    pizzf += -12 * (.50 - 2 * sw2 / 3) * (2 * sw2 / 3) * (pow(mt, 2) * creal(su_b0(qsz, mt, mt)) + pow(mcq, 2) * creal(su_b0(qsz, mcq, mcq)) + pow(mup, 2) * creal(su_b0(qsz, mup, mup)));
    pizzf += -12 * (-.50 + sw2 / 3) * (-sw2 / 3) * (pow(mb, 2) * creal(su_b0(qsz, mb, mb)) + pow(ms, 2) * creal(su_b0(qsz, ms, ms)) + pow(mdo, 2) * creal(su_b0(qsz, mdo, mdo)));
    pizzf += -4 * (-.50 + sw2) * (-sw2) * (pow(me, 2) * creal(su_b0(qsz, me, me)) + pow(mmu, 2) * creal(su_b0(qsz, mmu, mmu)) + pow(mtau, 2) * creal(su_b0(qsz, mtau, mtau)));


    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Contribution from W and Z (3rd line and the first part of 4th line of Eq. D4).
    // -------------------------------------------------------------------------    
    pizzb = -2 * pow(cw, 4) * (2 * qsz + pow(mw, 2) - pow(mz, 2) * pow(sw, 4) / pow(cw, 2)) * creal(su_b0(qsz, mw, mw)) - (8 * pow(cw, 4) + pow((cw2 - sw2), 2)) * creal(su_bt22(qsz, mw, mw));

    //
    // This line is added and substracted later in SUSPECT code. I am not completely sure about this line. 
    // So I am keeping it. But I don't understand it properly.  
    // -----------------------------------------------------------------------------------------------------
    pizzh0 = -creal(su_bt22(qsz, mz, ml)) - pow(mz, 2) * creal(su_b0(qsz, mz, ml));   

    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Contribution from sfermions (1st, 2nd and the second part of 4th line of Eq. D4).
    // -------------------------------------------------------------------------

    double ml = hml, mh = hmh, ma = hma, mch = hmch;
    pizzhs = -pow(sin(beta - alfa), 2) * (creal(su_bt22(qsz, ma, mh)) + creal(su_bt22(qsz, mz, ml)) - pow(mz, 2) * creal(su_b0(qsz, mz, ml)));
    pizzhs += -pow(cos(beta - alfa), 2) * (creal(su_bt22(qsz, mz, mh)) + creal(su_bt22(qsz, ma, ml)) - pow(mz, 2) * creal(su_b0(qsz, mz, mh)));
    pizzhs += -pow((pow(cw, 2) - pow(sw, 2)), 2) * creal(su_bt22(qsz, mch, mch)) - pizzh0;

    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  gf = If3 - ef sw^2        (Check A7 and A8)
    //  Contribution from sfermions (5th line of Eq. D4).
    // -------------------------------------------------------------------------
    mst1 = mst1sfbp;
    mst2 = mst2sfbp;

    pizzsu = -12 * pow(((.50 - 2 * sw2 / 3) * pow(cos(thet), 2) - (2 * sw2 / 3) * pow(sin(thet), 2)), 2) * creal(su_bt22(qsz, mst1, mst1));
    pizzsu += -12 * pow((-(.50 - 2 * sw2 / 3) * pow(sin(thet), 2) + (2 * sw2 / 3) * pow(cos(thet), 2)), 2) * creal(su_bt22(qsz, mst2, mst2));
    pizzsu += -24 * pow(((.50) * sin(thet) * cos(thet)), 2) * creal(su_bt22(qsz, mst1, mst2));
    pizzsu += -24 * pow((.50 - 2 * sw2 / 3), 2) * creal(su_bt22(qsz, msu1, msu1)) - 24 * pow((+2 * sw2 / 3), 2) * creal(su_bt22(qsz, msu2, msu2));

    pizzsd = -12 * pow(((-.50 + sw2 / 3) * pow(cos(theb), 2) - (-sw2 / 3) * pow(sin(theb), 2)), 2) * creal(su_bt22(qsz, msb1, msb1));
    pizzsd += -12 * pow((-(-.50 + sw2 / 3) * pow(sin(theb), 2) + (-sw2 / 3) * pow(cos(theb), 2)), 2) * creal(su_bt22(qsz, msb2, msb2));
    pizzsd += -24 * pow(((-0.50) * sin(theb) * cos(theb)), 2) * creal(su_bt22(qsz, msb1, msb2));
    pizzsd += -24 * pow((-.50 + sw2 / 3), 2) * creal(su_bt22(qsz, msd1, msd1)) - 24 * pow((-sw2 / 3), 2) * creal(su_bt22(qsz, msd2, msd2));

    pizzsl = -4 * pow(((-.50 + sw2) * pow(cos(thel), 2) - (-sw2) * pow(sin(thel), 2)), 2) * creal(su_bt22(qsz, msta1, msta1));
    pizzsl += -4 * pow((-(-.50 + sw2) * pow(sin(thel), 2) + (-sw2) * pow(cos(thel), 2)), 2) * creal(su_bt22(qsz, msta2, msta2));
    pizzsl += -8 * pow(((-.50) * sin(thel) * cos(thel)), 2) * creal(su_bt22(qsz, msta1, msta2));
    pizzsl += -8 * pow((-.50 + sw2), 2) * creal(su_bt22(qsz, mse1, mse1));
    pizzsl += -8 * pow((-sw2), 2) * creal(su_bt22(qsz, mse2, mse2));

    pizzsl += -8 * pow((.50), 2) * creal(su_bt22(qsz, msn1, msn1)) - 4 * pow((.50), 2) * creal(su_bt22(qsz, msntau, msntau));

    double mampizzs;
    mampizzs = pizzsl + pizzsd + pizzsu;

    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Correction from Neutralinos (2nd last line of Eq. D4).             
    // -------------------------------------------------------------------------
    
    pizzn = 0.0;

    for (int i = 1; i <= 4; i++)
        for (int j = 1; j <= 4; j++)
            pizzn = pizzn + 1.0 / 4 * pow((z[i][3] * z[j][3] - z[i][4] * z[j][4]), 2) * (creal(su_bh(qsz, gmn[i], gmn[j])) - 2 * dxmn[i] * dxmn[j] * creal(su_b0(qsz, gmn[i], gmn[j])));

    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Correction from Charginos (last line of Eq. D4).             
    // -------------------------------------------------------------------------

    pizzc = 0.0;
    double temp;

    for (int i = 1; i <= 2; i++)
        for (int j = 1; j <= 2; j++)
        {
            temp = 0.0;
            temp += (pow((2 * cw2 * v[i][1] * v[j][1] + (cw2 - sw2) * v[i][2] * v[j][2]), 2) + pow((2 * cw2 * u[i][1] * u[j][1] + (cw2 - sw2) * u[i][2] * u[j][2]), 2)) * creal(su_bh(qsz, gmc[i], gmc[j]));
            temp += +4 * (2 * cw2 * v[i][1] * v[j][1] + (cw2 - sw2) * v[i][2] * v[j][2]) * (2 * cw2 * u[i][1] * u[j][1] + (cw2 - sw2) * u[i][2] * u[j][2]) * gmc[i] * gmc[j] * creal(su_b0(qsz, gmc[i], gmc[j]));
            pizzc = pizzc + 1.0 / 4 * temp;
        }

    //
    // Sum of the susy contributions for pizz and final pizz(pow(mz,2))
    //----------------------------------------------------------------------------------------------------------------------------

    double alph = alpha;
    pizzsm = alph / 4 / pi / sw2 / cw2 * (pizzf + pizzb + pizzh0);
    pizzsusy = alph / 4 / pi / sw2 / cw2 * (pizzhs + mampizzs + pizzn + pizzc);
    *pizz = pizzsm + pizzsusy;


    /*

    //
    // Z boson self-energy at pow(q,2)=0
    //----------------------------------------------------------------------------------------------------------------------------

    qsz = eps;

    pizzf0 = 3 * (pow((.50 - 2 * sw2 / 3), 2) + pow((2 * sw2 / 3), 2)) * (creal(su_bh(qsz, mt, mt)) + creal(su_bh(qsz, mcq, mcq)) + creal(su_bh(qsz, mup, mup)));
    pizzf0 += 3 * (pow((-.50 + sw2 / 3), 2) + pow((-sw2 / 3), 2)) * (creal(su_bh(qsz, mb, mb)) + creal(su_bh(qsz, ms, ms)) + creal(su_bh(qsz, mdo, mdo)));
    pizzf0 += +(pow((-.50 + sw2), 2) + pow((-sw2), 2)) * (creal(su_bh(qsz, me, me)) + creal(su_bh(qsz, mmu, mmu)) + creal(su_bh(qsz, mtau, mtau)));
    pizzf0 += +pow(.50, 2) * 3 * creal(su_bh(qsz, eps, eps)) - 12 * (.50 - 2 * sw2 / 3) * (2 * sw2 / 3) * (pow(mt, 2) * creal(su_b0(qsz, mt, mt)) + pow(mcq, 2) * creal(su_b0(qsz, mcq, mcq)) + pow(mup, 2) * creal(su_b0(qsz, mup, mup)));
    pizzf0 += -12 * (-.50 + sw2 / 3) * (-sw2 / 3) * (pow(mb, 2) * creal(su_b0(qsz, mb, mb)) + pow(ms, 2) * creal(su_b0(qsz, ms, ms)) + pow(mdo, 2) * creal(su_b0(qsz, mdo, mdo)));
    pizzf0 += -4 * (-.50 + sw2) * (-sw2) * (pow(me, 2) * creal(su_b0(qsz, me, me)) + pow(mmu, 2) * creal(su_b0(qsz, mmu, mmu)) + pow(mtau, 2) * creal(su_b0(qsz, mtau, mtau)));

    pizzb0 = -2 * pow(cw, 4) * (2 * qsz + pow(mw, 2) - pow(mz, 2) * pow(sw, 4) / pow(cw, 2)) * creal(su_b0(qsz, mw, mw)) - (8 * pow(cw, 4) + pow((cw2 - sw2), 2)) * creal(su_bt22(qsz, mw, mw));

    pizzh00 = -creal(su_bt22(qsz, mz, ml)) - pow(mz, 2) * creal(su_b0(qsz, mz, ml));

    pizzhs0 = -pow(sin(beta - alfa), 2) * (creal(su_bt22(qsz, ma, mh)) + creal(su_bt22(qsz, mz, ml)) - pow(mz, 2) * creal(su_b0(qsz, mz, ml)));
    pizzhs0 += -pow(cos(beta - alfa), 2) * (creal(su_bt22(qsz, mz, mh)) + creal(su_bt22(qsz, ma, ml)) - pow(mz, 2) * creal(su_b0(qsz, mz, mh)));
    pizzhs0 += -pow((pow(cw, 2) - pow(sw, 2)), 2) * creal(su_bt22(qsz, mch, mch)) - pizzh00;

    pizzsu0 = -12 * pow(((.50 - 2 * sw2 / 3) * pow(cos(thet), 2) - (2 * sw2 / 3) * pow(sin(thet), 2)), 2) * creal(su_bt22(qsz, mst1, mst1));
    pizzsu0 += -12 * pow((-(.50 - 2 * sw2 / 3) * pow(sin(thet), 2) + (2 * sw2 / 3) * pow(cos(thet), 2)), 2) * creal(su_bt22(qsz, mst2, mst2));
    pizzsu0 += -24 * pow(((.50) * sin(thet) * cos(thet)), 2) * creal(su_bt22(qsz, mst1, mst2));
    pizzsu0 += -24 * pow((.50 - 2 * sw2 / 3), 2) * creal(su_bt22(qsz, msu1, msu1)) - 24 * pow((+2 * sw2 / 3), 2) * creal(su_bt22(qsz, msu2, msu2));

    double pizzsd0 = -12 * pow(((-.50 + sw2 / 3) * pow(cos(theb), 2) - (-sw2 / 3) * pow(sin(theb), 2)), 2) * creal(su_bt22(qsz, msb1, msb1));
    pizzsd0 += -12 * pow((-(-.50 + sw2 / 3) * pow(sin(theb), 2) + (-sw2 / 3) * pow(cos(theb), 2)), 2) * creal(su_bt22(qsz, msb2, msb2));
    pizzsd0 += -24 * pow(((-0.50) * sin(theb) * cos(theb)), 2) * creal(su_bt22(qsz, msb1, msb2));
    pizzsd0 += -24 * pow((-.50 + sw2 / 3), 2) * creal(su_bt22(qsz, msd1, msd1)) - 24 * pow((-sw2 / 3), 2) * creal(su_bt22(qsz, msd2, msd2));

    double pizzsl0 = -4 * pow(((-.50 + sw2) * pow(cos(thel), 2) - (-sw2) * pow(sin(thel), 2)), 2) * creal(su_bt22(qsz, msta1, msta1));
    pizzsl0 += -4 * pow((-(-.50 + sw2) * pow(sin(thel), 2) + (-sw2) * pow(cos(thel), 2)), 2) * creal(su_bt22(qsz, msta2, msta2));
    pizzsl0 += -8 * pow(((-.50) * sin(thel) * cos(thel)), 2) * creal(su_bt22(qsz, msta1, msta2)) - 8 * pow((-.50 + sw2), 2) * creal(su_bt22(qsz, mse1, mse1));
    pizzsl0 += -8 * pow((-sw2), 2) * creal(su_bt22(qsz, mse2, mse2));

    // correction msn1-> msntau (jlk)
    pizzsl0 += -8 * pow((.50), 2) * creal(su_bt22(qsz, msn1, msn1)) - 4 * pow((.50), 2) * creal(su_bt22(qsz, msntau, msntau));

    // Chek what this line is doing here
    double pizzs0 = pizzsl0 + pizzsd0 + pizzsu0;

    double pizzn0 = 0.0;

    for (int i = 1; i <= 4; i++)
        for (int j = 1; j <= 4; j++)
            pizzn0 = pizzn0 + 1.0 / 4 * pow((z[i][3] * z[j][3] - z[i][4] * z[j][4]), 2) * (creal(su_bh(qsz, gmn[i], gmn[j])) - 2 * dxmn[i] * dxmn[j] * creal(su_b0(qsz, gmn[i], gmn[j])));

    double pizzc0 = 0.0;

    for (int i = 1; i <= 2; i++)
        for (int j = 1; j <= 2; j++)
        {
            temp = (pow((2 * cw2 * v[i][1] * v[j][1] + (cw2 - sw2) * v[i][2] * v[j][2]), 2) + pow((2 * cw2 * u[i][1] * u[j][1] + (cw2 - sw2) * u[i][2] * u[j][2]), 2)) * creal(su_bh(qsz, gmc[i], gmc[j]));
            temp += +4 * (2 * cw2 * v[i][1] * v[j][1] + (cw2 - sw2) * v[i][2] * v[j][2]) * (2 * cw2 * u[i][1] * u[j][1] + (cw2 - sw2) * u[i][2] * u[j][2]) * gmc[i] * gmc[j] * creal(su_b0(qsz, gmc[i], gmc[j]));
            pizzc0 = pizzc0 + 1.0 / 4 * temp;
        }

    //
    // Sum of the susy contributions for pizz and final pizz(pow(mz,2))
    //----------------------------------------------------------------------------------------------------------------------------

    double pizzsm0 = alph / 4 / pi / sw2 / cw2 * (pizzf0 + pizzb0 + pizzh00);
    double pizzsusy0 = alph / 4 / pi / sw2 / cw2 * (pizzhs0 + pizzs0 + pizzn0 + pizzc0); // Check something is wrong
    double pizz0 = pizzsm0 + pizzsusy0;

    */    


    //
    // W boson self-energy at pow(q,2)=pow(mw,2)
    //----------------------------------------------------------------------------------------------------------------------------

    double qsw = pow(mw, 2); // Check something wrong here..

    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Correction from SM particles (2nd last line first part of Eq D9).             
    // -------------------------------------------------------------------------

    temp = 3.0 / 2 * (creal(su_bh(qsw, mt, mb)) + creal(su_bh(qsw, mcq, ms)));
    temp += 3.0 / 2 * (creal(su_bh(qsw, mup, mdo)));
    temp += 0.50 * (creal(su_bh(qsw, me, eps)) + creal(su_bh(qsw, mmu, eps)));
    temp += 0.50 * creal(su_bh(qsw, mtau, eps));
    double piwwf = temp;

    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Correction from W and Z (2nd part of 3rd line and 4th aand 5th line of Eq D9).             
    // -------------------------------------------------------------------------
    double piwwb = -(1.0 + 8 * pow(cw, 2)) * creal(su_bt22(qsw, mz, mw));
    piwwb += -pow(sw, 2) * (8 * creal(su_bt22(qsw, mw, eps)) + 4 * qsw * creal(su_b0(qsw, mw, eps)));
    piwwb += -((4 * qsw + pow(mz, 2) + pow(mw, 2)) * pow(cw, 2) - pow(mz, 2) * pow(sw, 4)) * creal(su_b0(qsw, mz, mw));


    // This line is added in substracted later in SUSPECT code. I am not completely sure about this line. 
    // So I am keeping it. But I don't understand it properly.  
    // -----------------------------------------------------------------------------------------------------
    double piwwh0 = -creal(su_bt22(qsw, ml, mw)) - pow(mw, 2) * creal(su_b0(qsw, ml, mw));


    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Correction from Hhiggs etc (1st and 2nd line and first part of 3rd line of Eq. D9).             
    // -------------------------------------------------------------------------

    double piwwhs = -pow(sin(beta - alfa), 2) * (creal(su_bt22(qsw, mh, mch)) + creal(su_bt22(qsw, ml, mw)) - pow(mw, 2) * creal(su_b0(qsw, ml, mw)));
    piwwhs += -pow(cos(beta - alfa), 2) * (creal(su_bt22(qsw, ml, mch)) + creal(su_bt22(qsw, mh, mw)) - pow(mw, 2) * creal(su_b0(qsw, mh, mw)));
    piwwhs += -creal(su_bt22(qsw, ma, mch)) - piwwh0;


    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Correction from SUSY particle (2nd last line of Eq. D9).             
    // -------------------------------------------------------------------------

    temp = 2 * creal(su_bt22(qsw, msu1, msd1)) + pow(cos(thet), 2) * pow(cos(theb), 2) * creal(su_bt22(qsw, mst1, msb1));
    temp += +pow(cos(thet), 2) * pow(sin(theb), 2) * creal(su_bt22(qsw, mst1, msb2)) + pow(sin(thet), 2) * pow(cos(theb), 2) * creal(su_bt22(qsw, mst2, msb1));
    temp += +pow(sin(thet), 2) * pow(sin(theb), 2) * creal(su_bt22(qsw, mst2, msb2));
    double piwws = -2 * 3 * temp;
    piwws = piwws - 2 * (2 * creal(su_bt22(qsw, msn1, mse1)) + pow(cos(thel), 2) * creal(su_bt22(qsw, msntau, msta1)) + pow(sin(thel), 2) * creal(su_bt22(qsw, msntau, msta2)));

    
    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Correction from Charginos and Neutralinos (last line of Eq. D9).             
    // -------------------------------------------------------------------------

    double piwwnc = 0.0;
    for (int i = 1; i <= 4; i++)
        for (int j = 1; j <= 2; j++)
        {
            piwwnc = piwwnc + (pow((-z[i][2] * v[j][1] + z[i][4] * v[j][2] / sq2), 2) + pow((-z[i][2] * u[j][1] - z[i][3] * u[j][2] / sq2), 2)) * creal(su_bh(qsw, gmn[i], gmc[j]));
            piwwnc = piwwnc + 4 * (-z[i][2] * v[j][1] + z[i][4] * v[j][2] / sq2) * (-z[i][2] * u[j][1] - z[i][3] * u[j][2] / sq2) * dxmn[i] * gmc[j] * creal(su_b0(qsw, gmn[i], gmc[j]));
        }

    //
    // Sum of the susy contributions for piww and final piww(pow(mw,2))
    //----------------------------------------------------------------------------------------------------------------------------

    double piwwsm = alph / 4 / pi / sw2 * (piwwf + piwwb + piwwh0);
    double piwwsusy = alph / 4 / pi / sw2 * (piwwhs + piwws + piwwnc);
    *piww = piwwsm + piwwsusy;

    //
    // W boson self-energy at pow(q,2)=0
    //----------------------------------------------------------------------------------------------------------------------------

    qsw = eps;

    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Correction from SM particles (2nd last line first part of Eq D9).             
    // -------------------------------------------------------------------------

    double piwwf0 = 3.0 / 2 * (creal(su_bh(qsw, mt, mb)) + creal(su_bh(qsw, mcq, ms)) + creal(su_bh(qsw, mup, mdo)));
    piwwf0 += +0.50 * (creal(su_bh(qsw, me, eps)) + creal(su_bh(qsw, mmu, eps)) + creal(su_bh(qsw, mtau, eps)));


    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Correction from W and Z (2nd part of 3rd line and 4th aand 5th line of Eq D9).             
    // -------------------------------------------------------------------------    

    double piwwb0 = -(1.0 + 8 * pow(cw, 2)) * creal(su_bt22(qsw, mz, mw)) - pow(sw, 2) * (8 * creal(su_bt22(qsw, mw, eps)) + 4 * qsw * creal(su_b0(qsw, mw, eps)));
    piwwb0 += -((4 * qsw + pow(mz, 2) + pow(mw, 2)) * pow(cw, 2) - pow(mz, 2) * pow(sw, 4)) * creal(su_b0(qsw, mz, mw));


    // This line is added in substracted later in SUSPECT code. I am not completely sure about this line. 
    // So I am keeping it. But I don't understand it properly.  
    // -----------------------------------------------------------------------------------------------------

    double piwwh00 = -creal(su_bt22(qsw, ml, mw)) - pow(mw, 2) * creal(su_b0(qsw, ml, mw));

    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Correction from Hhiggs etc (1st and 2nd line and first part of 3rd line of Eq. D9).             
    // -------------------------------------------------------------------------

    double piwwhs0 = -pow(sin(beta - alfa), 2) * (creal(su_bt22(qsw, mh, mch)) + creal(su_bt22(qsw, ml, mw)) - pow(mw, 2) * creal(su_b0(qsw, ml, mw)));
    piwwhs0 += -pow(cos(beta - alfa), 2) * (creal(su_bt22(qsw, ml, mch)) + creal(su_bt22(qsw, mh, mw)) - pow(mw, 2) * creal(su_b0(qsw, mh, mw)));
    piwwhs0 += -creal(su_bt22(qsw, ma, mch)) - piwwh00;


    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Correction from SUSY particle (2nd last line of Eq. D9).             
    // -------------------------------------------------------------------------

    double piwws0 = 2 * creal(su_bt22(qsw, msu1, msd1));
    piwws0 += +pow(cos(thet), 2) * pow(cos(theb), 2) * creal(su_bt22(qsw, mst1, msb1));
    piwws0 += +pow(cos(thet), 2) * pow(sin(theb), 2) * creal(su_bt22(qsw, mst1, msb2));
    piwws0 += +pow(sin(thet), 2) * pow(cos(theb), 2) * creal(su_bt22(qsw, mst2, msb1));
    piwws0 += +pow(sin(thet), 2) * pow(sin(theb), 2) * creal(su_bt22(qsw, mst2, msb2));
    piwws0 = -2 * 3 * piwws0;
    piwws0 += -2 * (2 * creal(su_bt22(qsw, msn1, mse1)) + pow(cos(thel), 2) * creal(su_bt22(qsw, msntau, msta1)) + pow(sin(thel), 2) * creal(su_bt22(qsw, msntau, msta2)));


    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Correction from Charginos and Neutralinos (last line of Eq. D9).             
    // -------------------------------------------------------------------------

    double piwwnc0 = 0.0;
    for (int i = 1; i <= 4; i++)
        for (int j = 1; j <= 2; j++)
        {
            piwwnc0 = piwwnc0 + (pow((-z[i][2] * v[j][1] + z[i][4] * v[j][2] / sq2), 2) + pow((-z[i][2] * u[j][1] - z[i][3] * u[j][2] / sq2), 2)) * creal(su_bh(qsw, gmn[i], gmc[j]));
            piwwnc0 = piwwnc0 + +4 * (-z[i][2] * v[j][1] + z[i][4] * v[j][2] / sq2) * (-z[i][2] * u[j][1] - z[i][3] * u[j][2] / sq2) * dxmn[i] * gmc[j] * creal(su_b0(qsw, gmn[i], gmc[j]));
        }

    //
    // Sum of the susy contributions for piww and final piww(0)
    //----------------------------------------------------------------------------------------------------------------------------

    double piwwsm0 = alph / 4 / pi / sw2 * (piwwf0 + piwwb0 + piwwh00);
    double piwwsusy0 = alph / 4 / pi / sw2 * (piwwhs0 + piwws0 + piwwnc0);
    *piww0 = piwwsm0 + piwwsusy0;

    mt = mtsave;
}


//
//  Gluino : 
//  calculates the radiative correction to the gluino mass, delgino.
//  the input parameters at EWSB scale are:
//  alphas,m3: the strong coupling constant and the su[3] gaugino mass,
//  msu1,msu2,msd1,msd2,msb1,msb2,mst1,mst2: the squark masses.
//-------------------------------------------------------------------------------------------------------------------------------------

void su_ginocr(double alphas, double m3, double mb, double mt, double *delgino)
{

    if (kpole == 1)
        scale = sqrt(mst1 * mst2);

    double mu = 0.005;
    double md = 0.015;
    double ms = 0.190;
    double mc = 1.42;
    double msc1 = msu1;
    double msc2 = msu2;
    double mss1 = msd1;
    double mss2 = msd2;

    double mst1 = mst1bp;
    double mst2 = mst2bp;

    double sumb1 = 0.0;

    //            
    //  PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
    //  Implementation of the Eq. D44
    // -------------------------------------------------------------------------

    sumb1 += +creal(su_b1(pow(m3, 2), mu, msu1)) + creal(su_b1(pow(m3, 2), mu, msu2));
    sumb1 += +creal(su_b1(pow(m3, 2), md, msd1)) + creal(su_b1(pow(m3, 2), md, msd2));
    sumb1 += +creal(su_b1(pow(m3, 2), ms, mss1)) + creal(su_b1(pow(m3, 2), ms, mss2));
    sumb1 += +creal(su_b1(pow(m3, 2), mc, msc1)) + creal(su_b1(pow(m3, 2), mc, msc2));
    sumb1 += +creal(su_b1(pow(m3, 2), mb, msb1)) + creal(su_b1(pow(m3, 2), mb, msb2));
    sumb1 += +creal(su_b1(pow(m3, 2), mt, mst1)) + creal(su_b1(pow(m3, 2), mt, mst2));

    double sumb0 = mb * sin(2 * thebbp) * (creal(su_b0(pow(m3, 2), mb, msb1)) - creal(su_b0(pow(m3, 2), mb, msb2)));
    sumb0 += +mt * sin(2 * thetbp) * (creal(su_b0(pow(m3, 2), mt, mst1)) - creal(su_b0(pow(m3, 2), mt, mst2)));

    *delgino = 3 * alphas / pi / 4 * m3 * (5.e0 - 3 * log(pow(m3, 2) / pow(scale, 2))) - alphas / pi / 4 * m3 * sumb1 - alphas / pi / 4 * sumb0;
}


double su_frho(double x, double y)
{
    if (x == y)
        return (0.0);
    else
        return (x + y - 2 * x * y / (x - y) * log(x / y));
}

//
// Define first the functions to be used in the calculation. These functions are used by su_runningcp()
// Ref : Eq(C7)
// ------------------------------------------------------------------------------------------------------------------------------------------

double su_frhol(double r)
{
    double frhol1;
    frhol1 = 19.0 - 33.0 / 2 * r + 43.0 / 12 * pow(r, 2) + 7.0 / 120 * pow(r, 3) - pi * sqrt(r) * (4.0 - 3.0 / 2 * r + 3.0 / 32 * pow(r, 2) + pow(r, 3) / 256);
    frhol1 = -pow(pi, 2) * (2.0 - 2 * r + pow(r, 2) / 2) - log(r) * (3 * r - pow(r, 2) / 2);
    return (frhol1);
}

//
// Define first the functions to be used in the calculation. These functions are used by su_runningcp()
// Ref : Eq(C8)
// ------------------------------------------------------------------------------------------------------------------------------------------

double su_frhoh(double r)
{
    double frhoh1;
    frhoh1 = pow(log(r), 2) * (3.0 / 2 - 9.0 / r - 15.e0 / pow(r, 2) - 48.0 / pow(r, 3) - 168.0 / pow(r, 4) - 612.0 / pow(r, 5));
    frhoh1 = -log(r) * (27.0 / 2 + 4.0 / r - 125.0 / 4 / pow(r, 2) - 558.0 / 5 / pow(r, 3) - 8307.0 / 20 / pow(r, 4) - 109321.0 / 70 / pow(r, 5));
    frhoh1 = +pow(pi, 2) * (1.0 - 4.0 / r - 5.0 / pow(r, 2) - 16.0 / pow(r, 3) - 56.0 / pow(r, 4) - 204.0 / pow(r, 5));
    frhoh1 = +49.0 / 4 + 2.0 / 3 / r + 1613.0 / 48 / pow(r, 2) + 87.570 / pow(r, 3) + 341959.0 / 1200 / pow(r, 4) + 9737663.0 / 9800 / pow(r, 5);
    return (frhoh1);
}

//
// Here come the subroutines for the radiative corrections to electroweak parameters (s^2_w etc) in drbar scheme:
// This routine is for the running gauge couplings and sin^2theta_w a la bpmz
// -----------------------------------------------------------------------------------------------------------------------------------------------------------

void su_runningcp(double alphas, double mt, double rmt, double rmg, double tbeta, double pizz, double piww, double piww0, double *alphadr, double *alphasdr, double *sw2dr)
{
    //
    // First thing to calculate: the running of alpha, including only the contributions of susy particles and the one of the top quark.
    // The commons are as those of the hloop routine.
    // Note that alphamz and alphasmz are the input values read by suspect  i.e. alpha(mz) and alphas(mz) in the msbar scheme.
    // ------------------------------------------------------------------------------------------------------------------------------------------

    double mc1 = dmc1;
    double mc2 = dmc2;
    double alph = alpha;

    double beta = atan(tbeta);

    double ct = cos(thet);
    double st = sin(thet);
    double cb = cos(theb);
    double sb = sin(theb);
    double msb1 = msb1sfbp;
    double msb2 = msb2sfbp;
    mch = mchrunz;

    //
    // The DR electromagnetic coupling, Eq(C1)
    // ---------------------------------------------------------------------------------------------------------------

    double temp = 0.0;
    temp += -7*log(mw/mz);											// I think this part is missing in the SUSPECT code.  But check it properly
    temp += +16.0 / 9 * log(mt / mz) + log(mch / mz) / 3 + 4.0 / 9 * (2 * log(msu1 / mz) + 2 * log(msu2 / mz) + log(mst1 / mz) + log(mst2 / mz));
    temp += +1.0 / 9 * (2 * log(msd1 / mz) + 2 * log(msd2 / mz) + log(msb1 / mz) + log(msb2 / mz));
    temp += +1.0 / 3 * (2 * log(mse1 / mz) + 2 * log(mse2 / mz) + log(msta1 / mz) + log(msta2 / mz)) + 4.0 / 3 * (log(mc1 / mz) + log(mc2 / mz));
    double dalpha = 0.0682 - alph / (2 * pi) * temp;                // I think 0.682 is missing SUSPECT code.
    *alphadr = alph / (1.0 - dalpha); // I think this part is missing in the SUSYSOFT code.  But check it properly

    //
    // The DR electromagnetic coupling, Eq(2), Eq(3)
    // Related to g_3 couplings g_3^2(Mz)/ (4\pi) = \alpha_s(Mz) / (1 - \Delta \alpha_s)  
    // ---------------------------------------------------------------------------------------------------------------

    double dalphas;
    dalphas = 1.0 / 2 - 2.0 / 3 * log(rmt / mz) - 2 * log(fabs(rmg) / mz);
    dalphas += -1.0 / 6 * (2 * log(msu1 / mz) + 2 * log(msu2 / mz) + log(mst1 / mz) + log(mst2 / mz) + 2 * log(msd1 / mz) + 2 * log(msd2 / mz) + log(msb1 / mz) + log(msb2 / mz));
    dalphas = alphas / (2 * pi) * dalphas;
    *alphasdr = alphas / (1.0 - dalphas);

    //
    // Now we calculate sin^2theta_w in the DR-scheme. In fact it turns out that we don't need input value measured at lep
    // or for the w mass in principle. However, to define the tree level quantities we need the value of m_w.
    // ---------------------------------------------------------------------------------------------------------------

    double cw2 = 1.0 - sw2;

    //
    // We first start with the calculation of the rho parameter at two-loop
    // Two--loop qcd corrections to the t/b contributions.
    // Eq(C6)
    // ---------------------------------------------------------------------------------------------------------------

    double drho2f = *alphadr / 4 / pi / sw2 * alphas / pi * (-2.145 * pow(mt, 2) / pow(mw, 2) + 1.262 * log(mt / mz) - 2.240 - 0.85 * pow(mz, 2) / pow(mt, 2));

    //
    //	Two loop higgs contribution:
    // Eq(C9), Eq(C10)
    // ---------------------------------------------------------------------------------------------------------------

    double drhol = su_frhol(ml / mt);
    double drhoh, drhoa;

    if (ma / mt <= 1.90)
    {
        drhoh = su_frhol(mh / mt);
        drhoa = su_frhol(mh / mt);
    }
    else
    {
        drhoh = su_frhoh(mh / mt);
        drhoa = su_frhoh(mh / mt);
    }

    //
    // The corresponding G^2 m_t^4 higher-order contributions from the heavy Higgs bosons
    // However these contributions are negligable
    // Eq. C10 
    //----------------------------------------------------------------------------------------------------

    temp = 0.0;
    temp += pow(cos(alfa), 2) / pow(sin(beta), 2) * drhol;
    temp += +pow(sin(alfa), 2) / pow(sin(beta), 2) * drhoh;
    temp += -1.0 / pow(tbeta, 2) * drhoa;
    double drho2h = pow((3 * gf * pow(mt, 2) / 8 / pow(pi, 2) / sqrt(2.0)), 2) / 3 * temp;

    //
    // Two loop qcd corrections to the stop/sbottom contributions
    // ---------------------------------------------------------------------------------------------------------------

    double test = 0.0;
    test = pow(ct, 2) * (pow(cb, 2) * su_frho(mst1, msb1) + pow(sb, 2) * su_frho(mst1, msb2));
    test += +pow(st, 2) * (pow(cb, 2) * su_frho(mst2, msb1) + pow(sb, 2) * su_frho(mst2, msb2));
    test += -pow(ct, 2) * pow(st, 2) * su_frho(mst1, mst2) - pow(cb, 2) * pow(sb, 2) * su_frho(msb1, msb2);
    double drho2s = 3 * gf / (8 * pow(pi, 2) * sqrt(2.0)) * 2 * alphas / 3 / pi * (1.0 + pow(pi, 2) / 3) * test;



    double drho2 = drho2h + drho2f + drho2s;
    double drho = (pizz / pow(mz, 2) - piww / pow(mw, 2)) / (1.0 + pizz / pow(mz, 2)) + drho2;
    double rhohat = 1.0 / (1.0 - drho);

    //
    // Then we calculate deltar, first evaluating the higher orders
    // Ref : Eq(C5), Eq(C12)
    // ---------------------------------------------------------------------------------------------------------------

    double dr2f = (*alphadr) / 4 / pi / sw2 / cw2 * alphas / pi * (2.1450 * pow(mt, 2) / pow(mz, 2) + 0.5750 * log(mt / mz) - 0.2240 - 0.1440 * pow(mz, 2) / pow(mt, 2));
    double drvb = rhohat * (*alphadr) / 4 / pi / sw2 * (6.0 + log(cw2) / sw2 * (7.0 / 2 - 5.0 / 2 * sw2 - sw2 * (5.0 - 3.0 / 2)));

    //
    // (drvb is the non-universal pure sm contribution: its contribution is about 0.01 thus non-negligible in principle)
    // Ref : Eq(C6), Eq(C3)
    // ---------------------------------------------------------------------------------------------------------------

    double dr1loop = rhohat * piww0 / pow(mw, 2) - pizz / pow(mz, 2) + drvb;
    double dr2h = -drho2h * rhohat * (1.0 - dr1loop);
    double deltar = dr1loop + dr2h + dr2f;

    //
    // Now calculate sin^2(theta_w) by solving the usual equation:
    // Ref : Eq(C3)
    // ---------------------------------------------------------------------------------------------------------------

    double deter = (*alphadr) * pi / sqrt(2.0) / pow(mz, 2) / gf / (1.0 - deltar);
    *sw2dr = (1.0 - sqrt(1.0 - 4 * deter)) / 2;
}


//
// The main subroutine for the EWSB and calculates the tadpole corrections to the higgs mass terms squared.
//
//	The input are:
// 	q2: 				the scale at which ewsb is supposed to happen,
// 	mu: 				the higgsino parameter mu ar ewsb scale
// 	at,ab,al: 		the third generation trilinear couplings at ewsb scale
// 	ytau, yt, yb: 	the yukawa couplings (at ewsb scale)
// 	msta1,msta2,msb1,msb2,mst1,mst2,..,thet,theb,thel: masses and mixing of tau,b,top,.. etc sfermions at EWSB scale
//							(Input via global variables)
// 	Other important input parameters, such as the higgs, chargino, neutralino masses and couplings as well as SM
//		parameters are called via commons.
//
// The output are :
//	 	dvdvd2, dvdvu2, which are (up to some appropriate overall constants) the derivatives of the full one-loop scalar
//		potential including the contributions of all SM and susy particles a la pbmz (hep-ph/9606211).
//		Another output is pizz which allow to calculated the rc to pow(mz,2).
//
// The consistency of the EWSB mechanism is performed by the main program
// -------------------------------------------------------------------------------------------------------------------------------------

void su_vloop2(double q2, double mu, double at, double ab, double al, double *dvdvd2, double *dvdvu2, double *pizz1)
{
    double pizz = *pizz1;
    scale = sqrt(q2);                                                                   // Basic parameters and definitions used:
    double g = g2ewsb;
    double gstrong = sqrt(4.0 * pi * alsewsb);
    mz = zm;
    mw = wm;

    //
    // Defining s^2_w at EWSB scale:
    //----------------------------------------------------------------------------------------------------------------------------

    double cw = 1.0 / sqrt(1.0 + pow((g1ewsb / g2ewsb), 2));
    double sw = g1ewsb / g2ewsb * cw;
    double sw2 = pow(sw, 2);
    double cw2 = pow(cw, 2);
    double cwm2 = 1.0 / cw2;
    double alph = pow((g * sw), 2) / (4 * pi);

    //
    // Defining mtau,mb,mt running masses at EWSB scale:
    //----------------------------------------------------------------------------------------------------------------------------

    if (isnan(pizz) || pow(mz, 2) + pizz <= 0.0) // Protection from imagionary masses
    {
        pizz = 0.0;

        if (irge == irgmax)
            inonpert = -1;
    }

    double vd2 = 2 * (pow(mz, 2) + pizz) / (pow(g1ewsb, 2) + pow(g2ewsb, 2)) / (1.0 + pow(tbeta, 2));
    double rmz = sqrt(pow(mz, 2) + pizz);
    double vu2 = vd2 * pow(tbeta, 2);
    double vev2 = 2.0 * (vu2 + vd2);
    double vd = sqrt(vd2);
    double vu = sqrt(vu2);

    double rmtau = ytauewsb * vd;
    double rmb = ybewsb * vd;
    double rmt = ytewsb * vu;

    double rmw = rmz * cw;
    double mzsave = mz;
    double mwsave = mw;
    mz = rmz;
    mw = rmw;

    //
    // Left-Right mixing : Top, Bottom, Tau
    //----------------------------------------------------------------------------------------------------------------------------
    // thet = thetbp;
    // theb = thebbp;
    // thel = thelbp;

    double ct = cos(thetbp);
    double st = sin(thetbp);
    double cb = cos(thebbp);
    double sb = sin(thebbp);
    double cta = cos(thelbp);
    double sta = sin(thelbp);

    beta = atan(tbeta);
    double cbeta2 = 1.0 / (1.0 + pow(tbeta, 2));
    double cbet = sqrt(cbeta2);
    double sbet = sqrt(1.0 - cbeta2);
    double c2b = 2 * cbeta2 - 1.0;
    double ma = ama0;
    double alfasave = alfa;                                                                     //  	 alfa running
    alfa = 0.5 * atan(tan(2.0 * beta) * (pow(ma, 2) + pow(mz, 2)) / (pow(ma, 2) - pow(mz, 2))); //	... Eq(A15)


    if (cos(2.0 * beta) * (pow(ma, 2) - pow(mz, 2)) > 0)
        alfa = alfa - pi / 2.0;

    double sal = sin(alfa);
    double cal = cos(alfa);
    double s2a = 2 * sal * cal;
    double tm = rmt;
    double bm = rmb;
    double taum = rmtau;

    //
    // yt, yb, ytau at EWSB scale are taken from global variables
    // Here the relevant sfermion couplings contributions (Feyman rules) are calculated:
    //----------------------------------------------------------------------------------------------------------------------------

    double sq2 = sqrt(2.0);
    double yb = ybewsb;
    double yt = ytewsb;
    double ytau = ytauewsb;

    // 
    // gdl = -1/2   -  -1/3 sw2
    // gdr =  0     -   1/3 sw2
    // gul =  1/2   -   2/3 sw2
    // gur =  0     -  -2/3 sw2
    // gel = -1/2   -  -1   sw2
    // ger =  0     -   1   sw2
    // gnl =  1/2   -   0   sw2
    //
    // Eq(A7), Eq(A8)
    //-------------------------------------------------------------------------------------------------------------------------------

    double s1blbl = g * mz / cw * (-.50 + sw2 / 3) * cbet + sq2 * yb * bm; // Eq(D53)
    double s1brbr = g * mz / cw * (-sw2 / 3) * cbet + sq2 * yb * bm;       // Eq(D53)
    double s1blbr = yb / sq2 * ab;                                         // Eq(D53)
    double s2blbl = -g * mz / cw * (-.50 + sw2 / 3) * sbet;                // Eq(D53)
    double s2brbr = -g * mz / cw * (-sw2 / 3) * sbet;                      // Eq(D53)
    double s2blbr = -yb / sq2 * mu;                                        // Eq(D53)

    double gs1d1d1 = g * mz / cw * (-.50 + sw2 / 3) * cbet;  // Eq(D53)
    double gs1d2d2 = g * mz / cw * (-sw2 / 3) * cbet;        // Eq(D53)
    double gs2d1d1 = -g * mz / cw * (-.50 + sw2 / 3) * sbet; // Eq(D53)
    double gs2d2d2 = -g * mz / cw * (-sw2 / 3) * sbet;       // Eq(D53)

    double s1taull = g * mz / cw * (-.50 + sw2) * cbet + sq2 * ytau * taum; // Eq(D53)
    double s1taurr = g * mz / cw * (-sw2) * cbet + sq2 * ytau * taum;       // Eq(D53)
    double s1taulr = ytau / sq2 * al;                                       // Eq(D53)
    double s2taull = -g * mz / cw * (-.50 + sw2) * sbet;                    // Eq(D53)
    double s2taurr = -g * mz / cw * (-sw2) * sbet;                          // Eq(D53)
    double s2taulr = -ytau / sq2 * mu;                                      // Eq(D53)

    double gs1e1e1 = g * mz / cw * (-.50 + sw2) * cbet;  // Eq(D53)
    double gs1e2e2 = g * mz / cw * (-sw2) * cbet;        // Eq(D53)
    double gs2e1e1 = -g * mz / cw * (-.50 + sw2) * sbet; // Eq(D53)
    double gs2e2e2 = -g * mz / cw * (-sw2) * sbet;       // Eq(D53)

    double s2tltl = -g * mz / cw * (.50 - 2 * sw2 / 3) * sbet + sq2 * yt * tm; // Eq(D53)
    double s2trtr = -g * mz / cw * (2 * sw2 / 3) * sbet + sq2 * yt * tm;       // Eq(D53)
    double s2tltr = yt / sq2 * at;                                             // Eq(D53)
    double s1tltl = g * mz / cw * (.50 - 2 * sw2 / 3) * cbet;                  // Eq(D53)
    double s1trtr = g * mz / cw * (2 * sw2 / 3) * cbet;                        // Eq(D53)
    double s1tltr = -yt / sq2 * mu;                                            // Eq(D53)

    double gs2u1u1 = -g * mz / cw * (.50 - 2 * sw2 / 3) * sbet; // Eq(D53)
    double gs2u2u2 = -g * mz / cw * (2 * sw2 / 3) * sbet;       // Eq(D53)
    double gs1u1u1 = g * mz / cw * (.50 - 2 * sw2 / 3) * cbet;  // Eq(D53)
    double gs1u2u2 = g * mz / cw * (2 * sw2 / 3) * cbet;        // Eq(D53)

    double gs2n1n1 = -g * mz / cw * (.50) * sbet; // Eq(D53)
    double gs1n1n1 = g * mz / cw * (.50) * cbet;  // Eq(D53)

    //
    // Change of Basis         .. Eq(D54)
    //----------------------------------------------------------------------------------------------------------------------------

    double gs1b1b1 = pow(cb, 2) * s1blbl + 2 * cb * sb * s1blbr + pow(sb, 2) * s1brbr;
    double gs1b2b2 = pow(sb, 2) * s1blbl - 2 * cb * sb * s1blbr + pow(cb, 2) * s1brbr;
    double gs2b1b1 = pow(cb, 2) * s2blbl + 2 * cb * sb * s2blbr + pow(sb, 2) * s2brbr;
    double gs2b2b2 = pow(sb, 2) * s2blbl - 2 * cb * sb * s2blbr + pow(cb, 2) * s2brbr;

    double gs1t1t1 = pow(ct, 2) * s1tltl + 2 * ct * st * s1tltr + pow(st, 2) * s1trtr;
    double gs1t2t2 = pow(st, 2) * s1tltl - 2 * ct * st * s1tltr + pow(ct, 2) * s1trtr;
    double gs2t1t1 = pow(ct, 2) * s2tltl + 2 * ct * st * s2tltr + pow(st, 2) * s2trtr;
    double gs2t2t2 = pow(st, 2) * s2tltl - 2 * ct * st * s2tltr + pow(ct, 2) * s2trtr;

    double gs1tau11 = pow(cta, 2) * s1taull + 2 * cta * sta * s1taulr + pow(sta, 2) * s1taurr;
    double gs1tau22 = pow(sta, 2) * s1taull - 2 * cta * sta * s1taulr + pow(cta, 2) * s1taurr;
    double gs2tau11 = pow(cta, 2) * s2taull + 2 * cta * sta * s2taulr + pow(sta, 2) * s2taurr;
    double gs2tau22 = pow(sta, 2) * s2taull - 2 * cta * sta * s2taulr + pow(cta, 2) * s2taurr;

    //
    //	Ncf = color factor = 3 for (s)quarks and
    //							 = 1 for (s)leptons
    // Sum_f  = Sum over all the quarks and leptons
    // Sum_fd = Sum over all the down type quarks and leptons d, s, b, e, mu and tau
    // Sum_fu = Sum over all the down type quarks and leptons u, c, t, nu_e, nu_mu, nu_tau
    //----------------------------------------------------------------------------------------------------------------------------

    // Down fermion contributions:
    //
    // 1st part of the first line of Eq. E4  (https://arxiv.org/pdf/hep-ph/9606211)
    // -----------------------------------------------------------------------------

    double confd = -2 * ytau * rmtau / vd * su_a(rmtau) - 2 * 3 * yb * rmb / vd * su_a(rmb);

    // Up fermion contributions:
    //
    // 1st part of the first line of Eq. E5  (https://arxiv.org/pdf/hep-ph/9606211)
    // -----------------------------------------------------------------------------
    double confu = -2 * 3 * yt * rmt / vu * su_a(rmt);

    // Vd sfermion contributions:

    //
    // 2nd part of the first line of Eq. E4  (https://arxiv.org/pdf/hep-ph/9606211)
    // -----------------------------------------------------------------------------
    double dvdsf = 0.0;
    dvdsf += 3 * gs1t1t1 * su_a(mst1bp) + 3 * gs1t2t2 * su_a(mst2bp);
    dvdsf += 3 * gs1b1b1 * su_a(msb1) + 3 * gs1b2b2 * su_a(msb2);
    dvdsf += gs1tau11 * su_a(msta1bp) + gs1tau22 * su_a(msta2bp);
    dvdsf += 2 * (3 * gs1u1u1 * su_a(msu1) + 3 * gs1u2u2 * su_a(msu2) + 3 * gs1d1d1 * su_a(msd1) + 3 * gs1d2d2 * su_a(msd2) + gs1e1e1 * su_a(mse1) + gs1e2e2 * su_a(mse2));
    dvdsf += 2 * gs1n1n1 * su_a(msn1bp) + gs1n1n1 * su_a(msntau);
    dvdsf = 1.0 / sq2 / vd * dvdsf;

    // Vu sfermion contributions:

    //
    // 2nd part of the first line of Eq. E5  (https://arxiv.org/pdf/hep-ph/9606211)
    // -----------------------------------------------------------------------------
    double dvusf = 0.0;
    dvusf += 3 * gs2t1t1 * su_a(mst1bp) + 3 * gs2t2t2 * su_a(mst2bp);
    dvusf += +3 * gs2b1b1 * su_a(msb1) + 3 * gs2b2b2 * su_a(msb2);
    dvusf += gs2tau11 * su_a(msta1bp) + gs2tau22 * su_a(msta2bp);
    dvusf += 2 * (3 * gs2u1u1 * su_a(msu1) + 3 * gs2u2u2 * su_a(msu2) + 3 * gs2d1d1 * su_a(msd1) + 3 * gs2d2d2 * su_a(msd2) + gs2e1e1 * su_a(mse1) + gs2e2e2 * su_a(mse2));
    dvusf += 2 * gs2n1n1 * su_a(msn1bp) + gs2n1n1 * su_a(msntau);
    dvusf = 1.0 / sq2 / vu * dvusf;

    //
    // 2nd and 3rd line of Eq. E4  (https://arxiv.org/pdf/hep-ph/9606211)
    // -----------------------------------------------------------------------------

    // Vd higgs contributions:
    double dvdh = 0.0;
    dvdh += -pow(g, 2) * c2b / pow(cw, 2) / 8 * (su_a(ma) + 2 * su_a(mch)) + pow(g, 2) / 2 * su_a(mch);
    dvdh += +pow(g, 2) / pow(cw, 2) / 8 * (3 * pow(sal, 2) - pow(cal, 2) + s2a * tbeta) * su_a(ml);
    dvdh += +pow(g, 2) / pow(cw, 2) / 8 * (3 * pow(cal, 2) - pow(sal, 2) - s2a * tbeta) * su_a(mh);

    // Vd gauge contributions:

    //
    // Last line of Eq. E4  (https://arxiv.org/pdf/hep-ph/9606211)
    // -----------------------------------------------------------------------------

    double dvdwz = 3 * pow(g, 2) / 4 * (2 * su_a(mw) + su_a(mz) / pow(cw, 2)) + pow(g, 2) * c2b / pow(cw, 2) / 8 * (2 * su_a(mw) + su_a(mz));

    // Vd gaugino contributions:
    double dvdino = 0.0;
    double mn1 = dmn1;
    double mn2 = dmn2;
    double mn3 = dmn3;
    double mn4 = dmn4;
    double mc1 = dmc1;
    double mc2 = dmc2;


    //
    // Line number 4 of Eq. E4  (https://arxiv.org/pdf/hep-ph/9606211)
    // -----------------------------------------------------------------------------

    dvdino += +dxmn[1] * z[1][3] * (z[1][2] - z[1][1] * sw / cw) * su_a(mn1);
    dvdino += +dxmn[2] * z[2][3] * (z[2][2] - z[2][1] * sw / cw) * su_a(mn2);
    dvdino += +dxmn[3] * z[3][3] * (z[3][2] - z[3][1] * sw / cw) * su_a(mn3);
    dvdino += +dxmn[4] * z[4][3] * (z[4][2] - z[4][1] * sw / cw) * su_a(mn4);
    dvdino = -pow(g, 2) / mw / cbet * dvdino;

    //
    // Line number 5 of Eq. E4  (https://arxiv.org/pdf/hep-ph/9606211)
    // -----------------------------------------------------------------------------
    dvdino += -sqrt(2.0) * pow(g, 2) / mw / cbet * (mc1 * v[1][1] * u[1][2] * su_a(mc1) + mc2 * v[2][1] * u[2][2] * su_a(mc2));

    //
    // Vu higgs contributions:
    
    //
    // 2nd and 3rd line of Eq. E5  (https://arxiv.org/pdf/hep-ph/9606211)
    // -----------------------------------------------------------------------------
    double dvuh = 0.0;
    dvuh += +pow(g, 2) * c2b / pow(cw, 2) / 8 * (su_a(ma) + 2 * su_a(mch)) + pow(g, 2) / 2 * su_a(mch);
    dvuh += +pow(g, 2) / pow(cw, 2) / 8 * (3 * pow(cal, 2) - pow(sal, 2) + s2a / tbeta) * su_a(ml);
    dvuh += +pow(g, 2) / pow(cw, 2) / 8 * (3 * pow(sal, 2) - pow(cal, 2) - s2a / tbeta) * su_a(mh);

    // Vu gauge contributions:

    //
    // Last line of Eq. E5  (https://arxiv.org/pdf/hep-ph/9606211)
    // -----------------------------------------------------------------------------
    double dvuwz = 3 * pow(g, 2) / 4 * (2 * su_a(mw) + su_a(mz) / pow(cw, 2));
    dvuwz += -pow(g, 2) * c2b / pow(cw, 2) / 8 * (2 * su_a(mw) + su_a(mz));

    // Vu neutralino contributions:

    //
    // Line number 4 of Eq. E5  (https://arxiv.org/pdf/hep-ph/9606211)
    // -----------------------------------------------------------------------------

    double dvuino = 0.0;
    dvuino += +dxmn[1] * z[1][4] * (z[1][2] - z[1][1] * sw / cw) * su_a(mn1);
    dvuino += +dxmn[2] * z[2][4] * (z[2][2] - z[2][1] * sw / cw) * su_a(mn2);
    dvuino += +dxmn[3] * z[3][4] * (z[3][2] - z[3][1] * sw / cw) * su_a(mn3);
    dvuino += +dxmn[4] * z[4][4] * (z[4][2] - z[4][1] * sw / cw) * su_a(mn4);
    dvuino = pow(g, 2) / mw / sbet * dvuino;

    //
    // Line number 5 of Eq. E5  (https://arxiv.org/pdf/hep-ph/9606211)
    // -----------------------------------------------------------------------------
    dvuino += -sqrt(2.0) * pow(g, 2) / mw / sbet * (mc1 * v[1][2] * u[1][1] * su_a(mc1) + mc2 * v[2][2] * u[2][1] * su_a(mc2));


    // 2-loop tadpole parts:
    double tad1st, tad2st, tad1sb, tad2sb, tad1w, tad2w, tad1bl, tad2bl, tad1l, tad2l;
    int imodel = ihflag;

    double am3 = m3;
    if (imodel >= 2)
    {
        // Two loop higgs part are not added in this code till now.. So we are doing the one loop corrections only..

        printf("\nThis program don't support two loop higgs contribution.\n So we are doing one loop corrections only.\n");
        tad1st = 0.0;
        tad2st = 0.0;
        tad1sb = 0.0;
        tad2sb = 0.0;
        tad1w = 0.0;
        tad2w = 0.0;
        tad1bl = 0.0;
        tad2bl = 0.0;
        tad1l = 0.0;
        tad2l = 0.0;
    }

    else
    {
        tad1st = 0.0;
        tad2st = 0.0;
        tad1sb = 0.0;
        tad2sb = 0.0;
        tad1w = 0.0;
        tad2w = 0.0;
        tad1bl = 0.0;
        tad2bl = 0.0;
        tad1l = 0.0;
        tad2l = 0.0;
    }

    //
    // Final contributions including 2-loop corrections:

    double tad1 = -cpi * (confd + dvdsf + dvdh + dvdwz + dvdino);
    double tad1loop = tad1st + tad1sb + tad1w + tad1l + tad1bl;
    *dvdvd2 = tad1 + tad1loop;

    double tad2 = -cpi * (confu + dvusf + dvuh + dvuwz + dvuino);
    double tad2loop = tad2st + tad2sb + tad2w + tad2l + tad2bl;
    *dvdvu2 = tad2 + tad2loop;

    //
    // Z boson self-energy at pow(q,2) = pow(mz,2)
    //	Eq(D1), Eq(D2), Eq(D3), Eq(D4), Eq(D5)
    // ------------------------------------------------------------------------------------------------------------------------------

    double qsz = pow(mzsave, 2);
    double mup = 1.0e-2;
    double mdo = 1.0e-2;
    double me = 0.5e-3;
    double mmu = 0.106;
    double ms = 0.1900;
    double mcq = 1.420;
    double mt = rmt;
    double mb = rmb;

    double eps = 1.0e-2;
    double eps0 = pow(eps, 2);
    gmn[1] = fabs(dxmn[1]);
    gmn[2] = fabs(dxmn[2]);
    gmn[3] = fabs(dxmn[3]);
    gmn[4] = fabs(dxmn[4]);
    gmc[1] = mc1;
    gmc[2] = mc2;

    //
    // Contribution from standard model particles (3rd from last of Eq. D4). (https://arxiv.org/pdf/hep-ph/9606211)
    // ---------------------------------------------------------------------------------------------------------------------------

    double pizzf; // Warning : May be some mistakes here
    pizzf = 0.0;
    pizzf += 3 * (pow((.50 - 2 * sw2 / 3), 2) + pow((2 * sw2 / 3), 2)) * (creal(su_bh(qsz, mt, mt)) + creal(su_bh(qsz, mcq, mcq)) + creal(su_bh(qsz, mup, mup)));
    pizzf += 3 * (pow((-.50 + sw2 / 3), 2) + pow((-sw2 / 3), 2)) * (creal(su_bh(qsz, mb, mb)) + creal(su_bh(qsz, ms, ms)) + creal(su_bh(qsz, mdo, mdo)));
    pizzf += (pow((-.50 + sw2), 2) + pow((-sw2), 2)) * (creal(su_bh(qsz, me, me)) + creal(su_bh(qsz, mmu, mmu)) + creal(su_bh(qsz, mtau, mtau)));
    pizzf += pow(.50, 2) * 3 * creal(su_bh(qsz, eps, eps));
    pizzf += -12 * (.50 - 2 * sw2 / 3) * (2 * sw2 / 3) * (pow(mt, 2) * creal(su_b0(qsz, mt, mt)) + pow(mcq, 2) * creal(su_b0(qsz, mcq, mcq)) + pow(mup, 2) * creal(su_b0(qsz, mup, mup)));
    pizzf += -12 * (-.50 + sw2 / 3) * (-sw2 / 3) * (pow(mb, 2) * creal(su_b0(qsz, mb, mb)) + pow(ms, 2) * creal(su_b0(qsz, ms, ms)) + pow(mdo, 2) * creal(su_b0(qsz, mdo, mdo)));
    pizzf += -4 * (-.50 + sw2) * (-sw2) * (pow(me, 2) * creal(su_b0(qsz, me, me)) + pow(mmu, 2) * creal(su_b0(qsz, mmu, mmu)) + pow(mtau, 2) * creal(su_b0(qsz, mtau, mtau)));


    double pizzhs, pizzs, pizzb, pizzh0, pizzsu, pizzd, pizzsl, pizzsd;

    //
    // Contribution from W and Z (3rd line and the first part of 4th line of Eq. D4). (https://arxiv.org/pdf/hep-ph/9606211)
    // ---------------------------------------------------------------------------------------------------------------------------

    pizzb = -2 * pow(cw, 4) * (2 * qsz + pow(mw, 2) - pow(mz, 2) * pow(sw, 4) / pow(cw, 2)) * creal(su_b0(qsz, mw, mw));
    pizzb += -(8 * pow(cw, 4) + pow((cw2 - sw2), 2)) * creal(su_bt22(qsz, mw, mw));

    //
    // This part is added and subtracted in SUSY. No idea why. 
    // ---------------------------------------------------------------------------------------------------------------------------    
    pizzh0 = -creal(su_bt22(qsz, mz, ml)) - pow(mz, 2) * creal(su_b0(qsz, mz, ml));

    //
    // 1st, 2nd and second part of the 4th line of D4 (https://arxiv.org/pdf/hep-ph/9606211)
    // ---------------------------------------------------------------------------------------------------------------------------

    pizzhs = -pow(sin(beta - alfa), 2) * (creal(su_bt22(qsz, ma, mh)) + creal(su_bt22(qsz, mz, ml)) - pow(mz, 2) * creal(su_b0(qsz, mz, ml)));
    pizzhs += -pow(cos(beta - alfa), 2) * (creal(su_bt22(qsz, mz, mh)) + creal(su_bt22(qsz, ma, ml)) - pow(mz, 2) * creal(su_b0(qsz, mz, mh)));
    pizzhs += -pow((pow(cw, 2) - pow(sw, 2)), 2) * creal(su_bt22(qsz, mch, mch)) - pizzh0;

    //
    // 5th line of D4 (https://arxiv.org/pdf/hep-ph/9606211)
    // ---------------------------------------------------------------------------------------------------------------------------

    pizzsu = -12 * pow(((.50 - 2 * sw2 / 3) * pow(cos(thetbp), 2) - (2 * sw2 / 3) * pow(sin(thetbp), 2)), 2) * creal(su_bt22(qsz, mst1bp, mst1bp));
    pizzsu += -12 * pow((-(.50 - 2 * sw2 / 3) * pow(sin(thetbp), 2) + (2 * sw2 / 3) * pow(cos(thetbp), 2)), 2) * creal(su_bt22(qsz, mst2bp, mst2bp));
    pizzsu += -24 * pow(((.50) * sin(thetbp) * cos(thetbp)), 2) * creal(su_bt22(qsz, mst1bp, mst2bp));
    pizzsu += -24 * pow((.50 - 2 * sw2 / 3), 2) * creal(su_bt22(qsz, msu1, msu1)) - 24 * pow((+2 * sw2 / 3), 2) * creal(su_bt22(qsz, msu2, msu2));

    pizzsd = -12 * pow(((-.50 + sw2 / 3) * pow(cos(thebbp), 2) - (-sw2 / 3) * pow(sin(thebbp), 2)), 2) * creal(su_bt22(qsz, msb1, msb1));
    pizzsd += -12 * pow((-(-.50 + sw2 / 3) * pow(sin(thebbp), 2) + (-sw2 / 3) * pow(cos(thebbp), 2)), 2) * creal(su_bt22(qsz, msb2, msb2));
    pizzsd += -24 * pow(((-0.50) * sin(thebbp) * cos(thebbp)), 2) * creal(su_bt22(qsz, msb1, msb2));
    pizzsd += -24 * pow((-.50 + sw2 / 3), 2) * creal(su_bt22(qsz, msd1, msd1)) - 24 * pow((-sw2 / 3), 2) * creal(su_bt22(qsz, msd2, msd2));

    pizzsl = -4 * pow(((-.50 + sw2) * pow(cos(thelbp), 2) - (-sw2) * pow(sin(thelbp), 2)), 2) * creal(su_bt22(qsz, msta1bp, msta1bp));
    pizzsl += -4 * pow((-(-.50 + sw2) * pow(sin(thelbp), 2) + (-sw2) * pow(cos(thelbp), 2)), 2) * creal(su_bt22(qsz, msta2bp, msta2bp));
    pizzsl += -8 * pow((-.50 * sin(thelbp) * cos(thelbp)), 2) * creal(su_bt22(qsz, msta1bp, msta2bp));
    pizzsl += -8 * pow((-.50 + sw2), 2) * creal(su_bt22(qsz, mse1, mse1)) - 8 * pow(-sw2, 2) * creal(su_bt22(qsz, mse2, mse2));
    pizzsl += -8 * pow(.50, 2) * creal(su_bt22(qsz, msn1bp, msn1bp)) - 4 * pow(.50, 2) * creal(su_bt22(qsz, msntau, msntau));

    pizzs = pizzsl + pizzsd + pizzsu;

    //
    // 2nd last line of D4 (https://arxiv.org/pdf/hep-ph/9606211)
    // ---------------------------------------------------------------------------------------------------------------------------

    double pizzn = 0.0;
    for (int i = 1; i <= 4; i++)
        for (int j = 1; j <= 4; j++)
            pizzn = pizzn + 1.0 / 4 * pow((z[i][3] * z[j][3] - z[i][4] * z[j][4]), 2) * (creal(su_bh(qsz, gmn[i], gmn[j])) - 2 * dxmn[i] * dxmn[j] * creal(su_b0(qsz, gmn[i], gmn[j])));

    //
    // Last line of D4 (https://arxiv.org/pdf/hep-ph/9606211)
    // ---------------------------------------------------------------------------------------------------------------------------

    double pizzc = 0.0;
    double pizzctemp = 0.0;
    for (int i = 1; i <= 2; i++)
        for (int j = 1; j <= 2; j++)
        {
            pizzctemp = 0.0;
            pizzctemp += (pow((2 * cw2 * v[i][1] * v[j][1] + (cw2 - sw2) * v[i][2] * v[j][2]), 2) + pow((2 * cw2 * u[i][1] * u[j][1] + (cw2 - sw2) * u[i][2] * u[j][2]), 2)) * creal(su_bh(qsz, gmc[i], gmc[j]));
            pizzctemp += +4 * (2 * cw2 * v[i][1] * v[j][1] + (cw2 - sw2) * v[i][2] * v[j][2]) * (2 * cw2 * u[i][1] * u[j][1] + (cw2 - sw2) * u[i][2] * u[j][2]) * gmc[i] * gmc[j] * creal(su_b0(qsz, gmc[i], gmc[j]));
            pizzc = pizzc + 1.0 / 4 * pizzctemp;
        }

    //
    // Sum of the susy contributions for pizz and final pizz(pow(mz,2))
    // ------------------------------------------------------------------------------------------------------------------------------
    double pizzsm = alph / 4.0 / pi / sw2 / cw2 * (pizzf + pizzb + pizzh0);
    double pizzsusy = alph / 4.0 / pi / sw2 / cw2 * (pizzhs + pizzs + pizzn + pizzc);
    pizz = pizzsm + pizzsusy;


    *pizz1 = pizz; // Wha have to pass this variable to the main program
    mz = mzsave;
    mw = mwsave;
    alfa = alfasave;

}

//
// Calculates radiative corrections to the two stop masses, including standard and susy-qcd corrections and the yukawa corrections.
//
// The input at the EWSB scale are, respectively: the strong coupling constant,
//																  the gluino mass,
//																  mu parameter,
//																  pseudoscalar higgs boson mass,
// 															  and trilinear couplings of the top and the bottom quarks.
// The outputs are the radiative corrections to the ll,lr,rr entries of the stop mass matrix.
//----------------------------------------------------------------------------------------------------------------------------

//
//	Precision Corrections in the minimal supersymmetric standard model
//		Damien M. Pierce, Jonathan A. Bagger, Konstaintin T. Matchev and Ren-Jie Zhang
//	arXiv:hep-ph/9606211v3
// ---------------------------------------------------------------------------------------------------------------------------

int su_stopcr(double pscale, double mu, double at, double ab, double m3, double *crll, double *crlr, double *crrr)
{
    double yt = ytewsb;
    double yb = ybewsb;
    double gmc[3], gmn[5], antr[5], antl[5];
    double bntl[5], bntr[5], actl[3], actr[3], bctl[3], bctr[3], fttll[5], gttll[5];
    double fttlr[5], gttlr[5], fttrr[5], gttrr[5], fbtll[3], gbtll[3], fbtlr[3];
    double gbtlr[3], fbtrr[3], gbtrr[3];

    double pizz = pizzp;

    double mst1 = mst1bp;
    double mst2 = mst2bp;
    double msta1 = msta1bp;
    double msta2 = msta2bp;
    double msn1 = msn1bp;

    // fix ren. scale (used in b0 functions):
    if (kpole == 1)
        scale = sqrt(mst1 * mst2);

    if (mst1 == 0.0 || msb1 == 0.0)
        return 0;

    //   (nb protection: means mst1 or msb1 undefined yet)
    double sq2 = sqrt(2.0);
    double g = g2ewsb;
    double g1 = g1ewsb;
    double alphas = alsewsb;
    double cw = 1.0 / sqrt(1.0 + pow((g1ewsb / g2ewsb), 2));
    double sw = g1ewsb / g2ewsb * cw;
    double sw2 = pow(sw, 2);

    double vd2 = 2 * (pow(mz, 2) + pizz) / (pow(g1ewsb, 2) + pow(g2ewsb, 2)) / (1.0 + pow(tbeta, 2));
    double vu2 = vd2 * pow(tbeta, 2);
    vd = sqrt(vd2);
    vu = sqrt(vu2);
    double rmb = yb * vd;
    double rmt = yt * vu;

    double mc1 = dmc1;
    double mc2 = dmc2;
    double zero = 1.e-2;
    gmn[1] = fabs(dxmn[1]);
    gmn[2] = fabs(dxmn[2]);
    gmn[3] = fabs(dxmn[3]);
    gmn[4] = fabs(dxmn[4]);
    gmc[1] = mc1;
    gmc[2] = mc2;
    double alfa = alfarunpwsb;
    double mh = mhrunpwsb;
    double ml = mlrunpwsb;
    double ma = marunpwsb;
    double mch = mchrunpwsb;
    double mz = zm;
    double mw = wm;

    //   printf("\nAlpha : %e .................................................... ",alfa);
    double thet = thetbp;
    double theb = thebbp;
    double thel = thelbp;

    double b = atan(tbeta);
    double cbeta2 = 1.0 / (1.0 + pow(tbeta, 2));
    double cbet = sqrt(cbeta2);
    double sbet = sqrt(1.0 - cbeta2);
    double sal = sin(alfa);
    double cal = cos(alfa);

    double ct = cos(thet);
    double st = sin(thet);
    double cb = cos(theb);
    double sb = sin(theb);
    double cta = cos(thel);
    double sta = sin(thel);

    //
    // Higgs couplings
    //-----------------------------------------------------------------------------------------------------------------

    //
    // Feynman rules associated with the CP-even-Higgs-sfermion-sfermion vertices
    //-----------------------------------------------------------------------------------------------------------------

    double s2tltl = -g * mz / cw * (.50 - 2 * sw2 / 3) * sbet + sq2 * yt * rmt; // D53 S2_uL_uL
    double s2trtr = -g * mz / cw * (2 * sw2 / 3) * sbet + sq2 * yt * rmt;       // D53 S2_uR_uR
    double s2tltr = yt / sq2 * at;                                              // D53 S2_uL_uR
    double s1tltl = g * mz / cw * (.50 - 2 * sw2 / 3) * cbet;                   // D53 S1_uL_uL
    double s1trtr = g * mz / cw * (2 * sw2 / 3) * cbet;                         // D53 S1_uR_uR
    double s1tltr = -yt / sq2 * mu;                                             // D53 S1_uL_uR

    double s2tlt1 = ct * s2tltl + st * s2tltr;  // D54
    double s2tlt2 = -st * s2tltl + ct * s2tltr; // D54
    double s2trt1 = ct * s2tltr + st * s2trtr;  // D54
    double s2trt2 = -st * s2tltr + ct * s2trtr; // D54
    double s1tlt1 = ct * s1tltl + st * s1tltr;  // D54
    double s1tlt2 = -st * s1tltl + ct * s1tltr; // D54
    double s1trt1 = ct * s1tltr + st * s1trtr;  // D54
    double s1trt2 = -st * s1tltr + ct * s1trtr; // D54

    double ghtlt1 = cal * s1tlt1 + sal * s2tlt1;  // D55
    double gltlt1 = -sal * s1tlt1 + cal * s2tlt1; // D55
    double ghtlt2 = cal * s1tlt2 + sal * s2tlt2;  // D55
    double gltlt2 = -sal * s1tlt2 + cal * s2tlt2; // D55	Rotating the basis
    double ghtrt1 = cal * s1trt1 + sal * s2trt1;  // D55
    double gltrt1 = -sal * s1trt1 + cal * s2trt1; // D55
    double ghtrt2 = cal * s1trt2 + sal * s2trt2;  // D55
    double gltrt2 = -sal * s1trt2 + cal * s2trt2; // D55

    double atltr = -yt / sq2 * (-mu * sbet - at * cbet); // D56
    double gtltr = +yt / sq2 * (-mu * cbet + at * sbet); // D56

    double gatlt1 = st * atltr;  //
    double gatlt2 = ct * atltr;  // These parts will  be required later
    double gatrt1 = -ct * atltr; //	for calculating top quark mass corrections
    double gatrt2 = st * atltr;  //

    double ggtlt1 = st * gtltr;  //
    double ggtlt2 = ct * gtltr;  // These parts will  be required later
    double ggtrt1 = -ct * gtltr; //	for calculating top quark mass corrections
    double ggtrt2 = st * gtltr;  //

    //
    // Feyman rules for the charged-Higgs-sfermion-sfermion vertices
    //-----------------------------------------------------------------------------------------------------------------

    double gctlbl = g * mw / sq2 * sin(2 * b) - yt * rmt * cbet - yb * rmb * sbet; // D57
    double gctrbr = -yt * rmb * cbet - yb * rmt * sbet;                            // D57
    double gctlbr = yb * (-mu * cbet - ab * sbet);                                 // D57
    double gctrbl = yt * (-mu * sbet - at * cbet);                                 // D57

    double ggtlbl = -g * mw / sq2 * cos(2 * b) - yt * rmt * sbet + yb * rmb * cbet; // D57
    double ggtrbr = 0.0;                                                            // D57
    double ggtlbr = yb * (-mu * sbet + ab * cbet);                                  // D57
    double ggtrbl = -yt * (-mu * cbet + at * sbet);                                 // D57

    double gctlb1 = cb * gctlbl + sb * gctlbr;  //	D58
    double gctlb2 = -sb * gctlbl + cb * gctlbr; //	D58	Left right rotation for the top and
    double gctrb1 = cb * gctrbl + sb * gctrbr;  //	D58	the bottom quark
    double gctrb2 = -sb * gctrbl + cb * gctrbr; //	D58

    double ggtlb1 = cb * ggtlbl + sb * ggtlbr;  //	D58	Left right rotation for the top and
    double ggtlb2 = -sb * ggtlbl + cb * ggtlbr; //	D58	the bottom quark
    double ggtrb1 = cb * ggtrbl + sb * ggtrbr;  //	D58
    double ggtrb2 = -sb * ggtrbl + cb * ggtrbr; //	D58

    //
    // Neutralino/chargino couplings:
    // Initiallize : D20, A8
    //-----------------------------------------------------------------------------------------------------------------

    double ap1tl = 0.0;
    double bp1tl = g1 / sqrt(2.0) * (1.0 / 3.0);
    double ap1tr = g1 / sqrt(2.0) * (-4.0 / 3.0);
    double bp1tr = 0.0;
    double ap2tl = 0.0;
    double bp2tl = sqrt(2.0) * g * (.50);
    double ap2tr = 0.0;
    double bp2tr = 0.0;
    double ap3tl = 0.0;
    double ap3tr = 0.0;
    double bp3tl = 0.0;
    double bp3tr = 0.0;
    double ap4tl = yt;
    double ap4tr = 0.0;
    double bp4tl = 0.0;
    double bp4tr = yt;
    double aw1tl = g;
    double bw1tl = 0.0;
    double aw1tr = 0.0;
    double bw1tr = 0.0;
    double aw2tl = 0.0;
    double bw2tl = -yb;
    double aw2tr = -yt;
    double bw2tr = 0.0;

    //
    // Calculate the coupling to the mass eigen states
    // Eq(D21), Eq(D22)
    //-----------------------------------------------------------------------------------------------------------------

    for (int i = 1; i <= 4; i++)
    {
        antr[i] = z[i][1] * ap1tr + z[i][2] * ap2tr + z[i][3] * ap3tr + z[i][4] * ap4tr;
        bntr[i] = z[i][1] * bp1tr + z[i][2] * bp2tr + z[i][3] * bp3tr + z[i][4] * bp4tr;
        antl[i] = z[i][1] * ap1tl + z[i][2] * ap2tl + z[i][3] * ap3tl + z[i][4] * ap4tl;
        bntl[i] = z[i][1] * bp1tl + z[i][2] * bp2tl + z[i][3] * bp3tl + z[i][4] * bp4tl;
    }

    for (int i = 1; i <= 2; i++)
    {
        actr[i] = v[i][1] * aw1tr + v[i][2] * aw2tr;
        bctr[i] = u[i][1] * bw1tr + u[i][2] * bw2tr;
        actl[i] = v[i][1] * aw1tl + v[i][2] * aw2tl;
        bctl[i] = u[i][1] * bw1tl + u[i][2] * bw2tl;
    }

    for (int i = 1; i <= 4; i++)
    {
        fttll[i] = antl[i] * antl[i] + bntl[i] * bntl[i];
        gttll[i] = bntl[i] * antl[i] + antl[i] * bntl[i];
        fttlr[i] = antl[i] * antr[i] + bntl[i] * bntr[i];
        gttlr[i] = bntl[i] * antr[i] + antl[i] * bntr[i];
        fttrr[i] = antr[i] * antr[i] + bntr[i] * bntr[i];
        gttrr[i] = bntr[i] * antr[i] + antr[i] * bntr[i];
    }

    for (int i = 1; i <= 2; i++)
    {
        fbtll[i] = actl[i] * actl[i] + bctl[i] * bctl[i];
        gbtll[i] = bctl[i] * actl[i] + actl[i] * bctl[i];
        fbtlr[i] = actl[i] * actr[i] + bctl[i] * bctr[i];
        gbtlr[i] = bctl[i] * actr[i] + actl[i] * bctr[i];
        fbtrr[i] = actr[i] * actr[i] + bctr[i] * bctr[i];
        gbtrr[i] = bctr[i] * actr[i] + actr[i] * bctr[i];
    }

    //
    // Chargino/Neutralino :
    //	Page : 43
    //--------------------------------------------------------------------------------------------------------------------

    //
    // Change in LL contribution:
    //--------------------------------------------------------------------------------------------------------------------

    //
    // First line of Eq D47 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------

    double crllqcd = 0.0;
    crllqcd = crllqcd + pow(st, 2) * (creal(su_bf(pow(pscale, 2), mst2, zero)) + su_a(mst2));
    crllqcd = crllqcd + pow(ct, 2) * (creal(su_bf(pow(pscale, 2), mst1, zero)) + su_a(mst1));
    crllqcd = crllqcd + 2 * creal(su_bg(pow(pscale, 2), m3, rmt));
    crllqcd = 16 * pi * alphas / 3 * (crllqcd);


    //
    // Second line of Eq D47 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------    
    double crllyuk = 0.0;
    crllyuk += +pow(yt, 2) * (pow(st, 2) * su_a(mst1) + pow(ct, 2) * su_a(mst2));
    crllyuk += +pow(yb, 2) * (pow(sb, 2) * su_a(msb1) + pow(cb, 2) * su_a(msb2));


    crllyuk += +0.50 * (pow(yt, 2) * pow(sal, 2) - pow(g, 2) * (.50 - 2 * sw2 / 3) / 2.0 / pow(cw, 2) * (-cos(2 * alfa))) * su_a(mh);
    crllyuk += +0.50 * (pow(yt, 2) * pow(cal, 2) - pow(g, 2) * (.50 - 2 * sw2 / 3) / 2.0 / pow(cw, 2) * (+cos(2 * alfa))) * su_a(ml);
    crllyuk += +0.50 * (pow(yt, 2) * pow(sbet, 2) - pow(g, 2) * (.50 - 2 * sw2 / 3) / 2.0 / pow(cw, 2) * (-cos(2 * b))) * su_a(mz);
    crllyuk += +0.50 * (pow(yt, 2) * pow(cbet, 2) - pow(g, 2) * (.50 - 2 * sw2 / 3) / 2.0 / pow(cw, 2) * (cos(2 * b))) * su_a(ma);

    crllyuk += +(pow(yb, 2) * pow(sbet, 2) + pow(g, 2) * ((.50 - 2 * sw2 / 3) / 2.0 / pow(cw, 2) - .50) * (-cos(2 * b))) * su_a(mch);
    crllyuk += +(pow(yb, 2) * pow(cbet, 2) + pow(g, 2) * ((.50 - 2 * sw2 / 3) / 2.0 / pow(cw, 2) - .50) * (cos(2 * b))) * su_a(mw);

    //
    // 4th line of Eq D47 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------

    crllyuk += +pow(ghtlt1, 2) * creal(su_b0(pow(pscale, 2), mh, mst1));
    crllyuk += +pow(ghtlt2, 2) * creal(su_b0(pow(pscale, 2), mh, mst2));
    crllyuk += +pow(gltlt1, 2) * creal(su_b0(pow(pscale, 2), ml, mst1));
    crllyuk += +pow(gltlt2, 2) * creal(su_b0(pow(pscale, 2), ml, mst2));
    crllyuk += +pow(ggtlt1, 2) * creal(su_b0(pow(pscale, 2), mz, mst1));
    crllyuk += +pow(ggtlt2, 2) * creal(su_b0(pow(pscale, 2), mz, mst2));
    crllyuk += +pow(gatlt1, 2) * creal(su_b0(pow(pscale, 2), ma, mst1));
    crllyuk += +pow(gatlt2, 2) * creal(su_b0(pow(pscale, 2), ma, mst2));

    crllyuk += +pow(gctlb1, 2) * creal(su_b0(pow(pscale, 2), mch, msb1));
    crllyuk += +pow(gctlb2, 2) * creal(su_b0(pow(pscale, 2), mch, msb2));
    crllyuk += +pow(ggtlb1, 2) * creal(su_b0(pow(pscale, 2), mw, msb1));
    crllyuk += +pow(ggtlb2, 2) * creal(su_b0(pow(pscale, 2), mw, msb2));



    //
    // 5th, 6th and 7th line of Eq D47 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------

    double crllgau = 0.0;
    crllgau += 4 * pow(g, 2) / pow(cw, 2) * pow((.50 - 2 * sw2 / 3), 2) * su_a(mz) + 2 * pow(g, 2) * su_a(mw);
    crllgau += +pow((2 * g1 * cw / 3), 2) * (pow(ct, 2) * creal(su_bf(pow(pscale, 2), mst1, zero)) + pow(st, 2) * creal(su_bf(pow(pscale, 2), mst2, zero)));
    crllgau += +pow(g, 2) / pow(cw, 2) * pow((.50 - 2 * sw2 / 3), 2) * (pow(ct, 2) * creal(su_bf(pow(pscale, 2), mst1, mz)) + pow(st, 2) * creal(su_bf(pow(pscale, 2), mst2, mz)));
    crllgau += +0.5 * pow(g, 2) * (pow(cb, 2) * creal(su_bf(pow(pscale, 2), msb1, mw)) + pow(sb, 2) * creal(su_bf(pow(pscale, 2), msb2, mw)));  // Need to check this like. I think a factor of 2 missing
    crllgau += +pow(g, 2) / 4 * (pow(ct, 2) * su_a(mst1) + pow(st, 2) * su_a(mst2) + 2 * (pow(cb, 2) * su_a(msb1) + pow(sb, 2) * su_a(msb2)));

    //
    // 8th line of Eq D47 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------

    double temp1;
    double crllhyp = 0.0;
    temp1 = 0.0;
    temp1 += +3.0 * (+.50) * (pow(ct, 2) * su_a(mst1) + pow(st, 2) * su_a(mst2));
    temp1 += +3.0 * (-.50) * (pow(cb, 2) * su_a(msb1) + pow(sb, 2) * su_a(msb2));
    temp1 += +(-.50) * (pow(cta, 2) * su_a(msta1) + pow(sta, 2) * su_a(msta2));
    temp1 += +6.0 * (+.50) * su_a(msu1) + 6.0 * (-.50) * su_a(msd1);
    temp1 += +2.0 * (-.50) * su_a(mse1);
    temp1 += +2.0 * (.50) * su_a(msn1) + (.50) * su_a(msntau);

    crllhyp += pow(g, 2) * 0.50 * temp1;
    crllhyp += +pow(g1, 2) / 4 * pow((1.0 / 3.0), 2) * (pow(ct, 2) * su_a(mst1) + pow(st, 2) * su_a(mst2));

    //
    // 9th line of Eq D47 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------    

    temp1 = 0.0;
    temp1 += 3.0 * (1.0 / 3.0) * (pow(ct, 2) * su_a(mst1) + pow(st, 2) * su_a(mst2));
    temp1 += +3.0 * (1.0 / 3.0) * (pow(cb, 2) * su_a(msb1) + pow(sb, 2) * su_a(msb2));
    temp1 += +(-1.0) * (pow(cta, 2) * su_a(msta1) + pow(sta, 2) * su_a(msta2));
    temp1 += +6.0 * (1.0 / 3.0) * su_a(msu1) + 6.0 * (1.0 / 3.0) * su_a(msd1);
    temp1 += +2.0 * (-1.0) * su_a(mse1);
    temp1 += +2.0 * (-1.0) * su_a(msn1) + (-1.0) * su_a(msntau);
    crllhyp += +pow(g1, 2) / 4 * (1.0 / 3.0) * temp1;
    

    temp1 = 0.0;
    temp1 += 3.0 * (-4.0 / 3.0) * (pow(st, 2) * su_a(mst1) + pow(ct, 2) * su_a(mst2));
    temp1 += +3.0 * (2.0 / 3.0) * (pow(sb, 2) * su_a(msb1) + pow(cb, 2) * su_a(msb2));
    temp1 += +(2.0) * (pow(sta, 2) * su_a(msta1) + pow(cta, 2) * su_a(msta2));
    temp1 += +6.0 * (-4.0 / 3.0) * su_a(msu2) + 6.0 * (2.0 / 3.0) * su_a(msd2);
    temp1 += +2.0 * (2.0) * su_a(mse2);
    crllhyp += +pow(g1, 2) / 4 * (1.0 / 3.0) * temp1;

    //
    // Second Last line of Eq D47 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------

    double crllnino = 0.0;
    for (int i = 1; i <= 4; i++)
        crllnino = crllnino + fttll[i] * creal(su_bg(pow(pscale, 2), gmn[i], rmt)) - 2.0 * rmt * dxmn[i] * gttll[i] * creal(su_b0(pow(pscale, 2), gmn[i], rmt));


    //
    // Last line of Eq D47 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------

    double crllcino = 0.0;
    for (int i = 1; i <= 2; i++)
        crllcino = crllcino + fbtll[i] * creal(su_bg(pow(pscale, 2), gmc[i], rmb)) - 2.0 * rmb * gmc[i] * gttll[i] * creal(su_b0(pow(pscale, 2), gmc[i], rmb));

    *crll = -cpi * (crllqcd + crllyuk + crllgau + crllhyp + crllnino + crllcino);

    //
    // Change in RR contribution:
    //--------------------------------------------------------------------------------------------------------------------

    //
    // First line of Eq D48 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------
    
    double crrrqcd = 0.0;
    temp1 = 0.0;
    temp1 += 2 * creal(su_bg(pow(pscale, 2), m3, rmt));
    temp1 += +pow(st, 2) * (creal(su_bf(pow(pscale, 2), mst1, zero)) + su_a(mst1));
    temp1 += +pow(ct, 2) * (creal(su_bf(pow(pscale, 2), mst2, zero)) + su_a(mst2));
    crrrqcd += 16 * pi * alphas / 3 * temp1;


    //
    // Second line of Eq D48 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------
    double crrryuk = 0.0;
    crrryuk += +pow(yt, 2) * (pow(ct, 2) * su_a(mst1) + pow(st, 2) * su_a(mst2));
    crrryuk += +pow(yt, 2) * (pow(cb, 2) * su_a(msb1) + pow(sb, 2) * su_a(msb2));

    //
    // 3rd line of Eq D48 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------
    crrryuk += +0.50 * (pow(yt, 2) * pow(sal, 2) - pow(g, 2) * (2 * sw2 / 3) / 2.0 / pow(cw, 2) * (-cos(2 * alfa))) * su_a(mh);
    crrryuk += +0.50 * (pow(yt, 2) * pow(cal, 2) - pow(g, 2) * (2 * sw2 / 3) / 2.0 / pow(cw, 2) * (+cos(2 * alfa))) * su_a(ml);
    crrryuk += +0.50 * (pow(yt, 2) * pow(sbet, 2) - pow(g, 2) * (2 * sw2 / 3) / 2.0 / pow(cw, 2) * (-cos(2 * b))) * su_a(mz);
    crrryuk += +0.50 * (pow(yt, 2) * pow(cbet, 2) - pow(g, 2) * (2 * sw2 / 3) / 2.0 / pow(cw, 2) * (cos(2 * b))) * su_a(ma);
    crrryuk += +(pow(yt, 2) * pow(cbet, 2) + pow(g, 2) * ((2 * sw2 / 3) / 2.0 / pow(cw, 2)) * (-cos(2 * b))) * su_a(mch);
    crrryuk += +(pow(yt, 2) * pow(sbet, 2) + pow(g, 2) * ((2 * sw2 / 3) / 2.0 / pow(cw, 2)) * (cos(2 * b))) * su_a(mw);

    //
    // 4th line of Eq D48 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------
    crrryuk += +pow(ghtrt1, 2) * creal(su_b0(pow(pscale, 2), mh, mst1));
    crrryuk += +pow(ghtrt2, 2) * creal(su_b0(pow(pscale, 2), mh, mst2));
    crrryuk += +pow(gltrt1, 2) * creal(su_b0(pow(pscale, 2), ml, mst1));
    crrryuk += +pow(gltrt2, 2) * creal(su_b0(pow(pscale, 2), ml, mst2));
    crrryuk += +pow(ggtrt1, 2) * creal(su_b0(pow(pscale, 2), mz, mst1));
    crrryuk += +pow(ggtrt2, 2) * creal(su_b0(pow(pscale, 2), mz, mst2));
    crrryuk += +pow(gatrt1, 2) * creal(su_b0(pow(pscale, 2), ma, mst1));
    crrryuk += +pow(gatrt2, 2) * creal(su_b0(pow(pscale, 2), ma, mst2));
    crrryuk += +pow(gctrb1, 2) * creal(su_b0(pow(pscale, 2), mch, msb1));
    crrryuk += +pow(gctrb2, 2) * creal(su_b0(pow(pscale, 2), mch, msb2));
    crrryuk += +pow(ggtrb1, 2) * creal(su_b0(pow(pscale, 2), mw, msb1));
    crrryuk += +pow(ggtrb2, 2) * creal(su_b0(pow(pscale, 2), mw, msb2));

    //
    // 5th and 6th line of Eq D48 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------
    double crrrgau = 0;
    crrrgau += 4 * pow(g, 2) / pow(cw, 2) * pow((2 * sw2 / 3), 2) * su_a(mz);
    crrrgau += +pow((2 * g1 * cw / 3), 2) * (pow(st, 2) * creal(su_bf(pow(pscale, 2), mst1, zero)) + pow(ct, 2) * creal(su_bf(pow(pscale, 2), mst2, zero)));
    crrrgau += +pow(g, 2) / pow(cw, 2) * pow((2 * sw2 / 3), 2) * (pow(st, 2) * creal(su_bf(pow(pscale, 2), mst1, mz)) + pow(ct, 2) * creal(su_bf(pow(pscale, 2), mst2, mz)));

    //
    // 7th line of Eq D48 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------    
    double crrrhyp = 0.0;
    crrrhyp += pow(g1, 2) / 4 * pow(-4.0 / 3.0, 2) * (pow(st, 2) * su_a(mst1) + pow(ct, 2) * su_a(mst2));

    //
    // 8th line of Eq D48 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------       
    temp1 = 0.0;
    temp1 += +3.0 * (1.0 / 3.0) * (pow(ct, 2) * su_a(mst1) + pow(st, 2) * su_a(mst2));
    temp1 += +3.0 * (1.0 / 3.0) * (pow(cb, 2) * su_a(msb1) + pow(sb, 2) * su_a(msb2));
    temp1 += +(-1.0) * (pow(cta, 2) * su_a(msta1) + pow(sta, 2) * su_a(msta2));
    temp1 += +6.0 * (1.0 / 3.0) * su_a(msu1) + 6.0 * (1.0 / 3.0) * su_a(msd1);
    temp1 += +2.0 * (-1.0) * su_a(mse1);
    temp1 += +2.0 * (-1.0) * su_a(msn1) + (-1.0) * su_a(msntau);
    crrrhyp += +pow(g1, 2) / 4 * (-4.0 / 3.0) * temp1;
    temp1 = 0.0;
    temp1 += +3.0 * (-4.0 / 3.0) * (pow(st, 2) * su_a(mst1) + pow(ct, 2) * su_a(mst2));
    temp1 += +3.0 * (2.0 / 3.0) * (pow(sb, 2) * su_a(msb1) + pow(cb, 2) * su_a(msb2));
    temp1 += +(2.0) * (pow(sta, 2) * su_a(msta1) + pow(cta, 2) * su_a(msta2));
    temp1 += +6.e0 * (-4.0 / 3.0) * su_a(msu2) + 6.e0 * (2.0 / 3.0) * su_a(msd2);
    temp1 += +2.0 * (2.0) * su_a(mse2);
    crrrhyp += +pow(g1, 2) / 4 * (1.0 / 3.0) * temp1;

    //
    // Second Last line of Eq D48 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------    
    double crrrnino = 0.0;
    for (int i = 1; i <= 4; i++)
        crrrnino = crrrnino + fttrr[i] * creal(su_bg(pow(pscale, 2), gmn[i], rmt)) - 2.0 * rmt * dxmn[i] * gttrr[i] * creal(su_b0(pow(pscale, 2), gmn[i], rmt));

    
    //
    // Last line of Eq D48 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------    
    double crrrcino = 0.0;
    for (int i = 1; i <= 2; i++)
        crrrcino = crrrcino + fbtrr[i] * creal(su_bg(pow(pscale, 2), gmc[i], rmb)) - 2.0 * rmb * gmc[i] * gttrr[i] * creal(su_b0(pow(pscale, 2), gmc[i], rmb));

    *crrr = -cpi * (crrrqcd + crrryuk + crrrgau + crrrhyp + crrrnino + crrrcino);

    //
    // Change in LR contribution:
    //--------------------------------------------------------------------------------------------------------------------

    //
    // First line of Eq D49 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------       
    double crlrqcd = 0.0;
    temp1 = 0.0;
    temp1 += 4 * rmt * m3 * creal(su_b0(pow(pscale, 2), m3, rmt));
    temp1 += +ct * st * (creal(su_bf(pow(pscale, 2), mst1, zero)) - su_a(mst1) - creal(su_bf(pow(pscale, 2), mst2, zero)) + su_a(mst2));
    crlrqcd += 16 * pi * alphas / 3 * temp1;

    //
    // Second line of Eq D49 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------       
    double crlryuk = 0.0;
    crlryuk += 3 * pow(yt, 2) * st * ct * (su_a(mst1) - su_a(mst2));
    crlryuk += +ghtlt1 * ghtrt1 * creal(su_b0(pow(pscale, 2), mh, mst1));
    crlryuk += +ghtlt2 * ghtrt2 * creal(su_b0(pow(pscale, 2), mh, mst2));
    crlryuk += +gltlt1 * gltrt1 * creal(su_b0(pow(pscale, 2), ml, mst1));
    crlryuk += +gltlt2 * gltrt2 * creal(su_b0(pow(pscale, 2), ml, mst2));
    crlryuk += +ggtlt1 * ggtrt1 * creal(su_b0(pow(pscale, 2), mz, mst1));
    crlryuk += +ggtlt2 * ggtrt2 * creal(su_b0(pow(pscale, 2), mz, mst2));
    crlryuk += +gatlt1 * gatrt1 * creal(su_b0(pow(pscale, 2), ma, mst1));
    crlryuk += +gatlt2 * gatrt2 * creal(su_b0(pow(pscale, 2), ma, mst2));
    crlryuk += +gctlb1 * gctrb1 * creal(su_b0(pow(pscale, 2), mch, msb1));
    crlryuk += +gctlb2 * gctrb2 * creal(su_b0(pow(pscale, 2), mch, msb2));
    crlryuk += +ggtlb1 * ggtrb1 * creal(su_b0(pow(pscale, 2), mw, msb1));
    crlryuk += +ggtlb2 * ggtrb2 * creal(su_b0(pow(pscale, 2), mw, msb2));

    double crlrgau = 0.0;
    crlrgau += pow((2 * g1 * cw / 3), 2) * ct * st * (creal(su_bf(pow(pscale, 2), mst1, zero)) - creal(su_bf(pow(pscale, 2), mst2, zero)));
    crlrgau += -pow(g, 2) / pow(cw, 2) * (.50 - 2 * sw2 / 3) * (2 * sw2 / 3) * st * ct * (creal(su_bf(pow(pscale, 2), mst1, mz)) - creal(su_bf(pow(pscale, 2), mst2, mz)));

    double crlrhyp = pow(g1, 2) / 4 * (1.0 / 3.0) * (-4.0 / 3.0) * st * ct * (su_a(mst1) - su_a(mst2));

    //
    // Second Last line of Eq D49 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------      
    double crlrnino = 0.0;
    for (int i = 1; i <= 4; i++)
        crlrnino = crlrnino + fttlr[i] * creal(su_bg(pow(pscale, 2), gmn[i], rmt)) - 2.0 * rmt * dxmn[i] * gttlr[i] * creal(su_b0(pow(pscale, 2), gmn[i], rmt));

    //
    // Last line of Eq D49 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------------------------------          
    double crlrcino = 0.0;
    for (int i = 1; i <= 2; i++)
        crlrcino = crlrcino + fbtlr[i] * creal(su_bg(pow(pscale, 2), gmc[i], rmb)) - 2.0 * rmb * gmc[i] * gttlr[i] * creal(su_b0(pow(pscale, 2), gmc[i], rmb));

    *crlr = -cpi * (crlrqcd + crlryuk + crlrgau + crlrhyp + crlrnino + crlrcino);

    return 0;
}

//
//	Calculates the radiative corrections to the top quark mass including the standard and susy qcd corrections
//	(the standard corrections are also calculable with runm) and the electroweak corrections including the
// contributions of gauge bosons, higgs bosons, Charginos and Neutralinos.
//
// The input are respectively: the strong coupling constant, the pole masses, running masses and yukawa couplings
//			of the top and bottom quarks, tan(beta), the 3d generation squark mass terms and trilinear couplings and mu.
//
// The output delmtop is the radiative correction to the top quark mass.
//-------------------------------------------------------------------------------------------------------------------------------------

void su_topmscr(double alphas, double mt, double mb, double rmt, double rmb, double yt, double yb, double tbeta, double mgl, double mql, double mur, double mdr, double at, double ab, double mu, double *delmtop)
{

    //
    // Basic parameters and definitions used
    //-----------------------------------------------------------------------------------------------------------------------------------

    //	tbeta = beta;
    double b = atan(tbeta);
    double sw = sqrt(sw2);
    double cw = sqrt(1.0 - sw2);
    double g2 = sqrt(g22);
    double e = g2 * sw;
    double cbeta2 = 1.0 / (1.0 + pow(tbeta, 2));
    double cbet = sqrt(cbeta2);
    double sbet = sqrt(1.0 - cbeta2);

    double sa = sin(alfa);
    double ca = cos(alfa);

    //
    // fix ren. scale (used in b0, b1): mz here.
    //----------------------------------------------------------------------------------------------------------------------------------

    double scalsave = scale;
    scale = mz;
    double ct = cos(thet);
    double st = sin(thet);
    double cb = cos(theb);
    double sb = sin(theb);

    //
    // defining couplings:
    //-----------------------------------------------------------------------------------------------------------------------------------

    double gtl = .50 - (2.0 / 3) * pow(sw, 2);
    double gtr = (2.0 / 3) * pow(sw, 2);

    //
    // Various contributions:
    //-----------------------------------------------------------------------------------------------------------------------------------
    double mst1 = mst1sfbp;
    double mst2 = mst2sfbp;


    //
    // 1st and 2nd line of D18 (https://arxiv.org/pdf/hep-ph/9606211). A quantity of (5+ 3ln(m/Q) is missing which is added later in this routine)
    // ------------------------------------------------------------------------------------------

    double temp = 0.0;
    temp += creal(su_b1(pow(mt, 2), mgl, mst1)) + creal(su_b1(pow(mt, 2), mgl, mst2));
    temp += -mgl / mt * sin(2 * thet) * (creal(su_b0(pow(mt, 2), mgl, mst1)) - creal(su_b0(pow(mt, 2), mgl, mst2)));
    double dqcd = 16 * pi * alphas / 3 * temp;

    //
    // Mgluino/mt*sin(2*thet)*
    //-----------------------------------------------------------------------------------------------------------------------------------

    //
    // 3rd and 4th line of D18 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------

    mh = mhrunz;
    ml = mlrunz;

    temp = 0.0;
    temp += pow(sa, 2) * (creal(su_b1(pow(mt, 2), rmt, mh)) + creal(su_b0(pow(mt, 2), rmt, mh)));
    temp += +pow(ca, 2) * (creal(su_b1(pow(mt, 2), rmt, ml)) + creal(su_b0(pow(mt, 2), rmt, ml)));
    temp += +pow(cbet, 2) * (creal(su_b1(pow(mt, 2), rmt, ma)) - creal(su_b0(pow(mt, 2), rmt, ma)));
    temp += +pow(sbet, 2) * (creal(su_b1(pow(mt, 2), rmt, mz)) - creal(su_b0(pow(mt, 2), rmt, mz)));
    double dyuk = pow(yt, 2) / 2 * temp;


    //
    // 5th line and first part of 6th line of D18 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------

    double mhp = mchrunz;
    temp = 0.0;
    temp += (pow(yb, 2) * pow(sbet, 2) + pow(yt, 2) * pow(cbet, 2)) * creal(su_b1(pow(mt, 2), mb, mhp));
    temp += (pow(g2, 2) + pow(yb, 2) * pow(cbet, 2) + pow(yt, 2) * pow(sbet, 2)) * creal(su_b1(pow(mt, 2), mb, mw));
    dyuk += +.5 * temp;
    dyuk += +pow(yb, 2) * pow(cbet, 2) * (creal(su_b0(pow(mt, 2), mb, mhp)) - creal(su_b0(pow(mt, 2), mb, mw)));


    //
    // Second part of 6th line and 7th line of D18 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------

    double dgauge = -pow(e, 2) * pow((2.0 / 3), 2) * (5.0 + 3 * log(pow(mz, 2) / pow(mt, 2)));
    dgauge += +pow(g2, 2) / pow(cw, 2) * ((pow(gtl, 2) + pow(gtr, 2)) * creal(su_b1(pow(mt, 2), mt, mz)) + 4 * gtl * gtr * creal(su_b0(pow(mt, 2), mt, mz)));


    double sq2 = sqrt(2.0);
    double ytr = -4.0 / 3;
    double ytl = 1.0 / 3;
    double ybr = 2.0 / 3;
    double ybl = 1.0 / 3;

    double ap1tl = 0.0;
    double bp1tl = e / cw / sq2 * ytl;
    double ap1tr = e / cw / sq2 * ytr;
    double bp1tr = 0.0;
    double ap2tl = 0.0;
    double bp2tl = sq2 * g2 * (0.50);
    double ap2tr = 0.0;
    double bp2tr = 0.0;
    double ap3tl = 0.0;
    double ap3tr = 0.0;
    double bp3tl = 0.0;
    double bp3tr = 0.0;
    double ap4tl = yt;
    double ap4tr = 0.0;
    double bp4tl = 0.0;
    double bp4tr = yt;


    double antr[5], bntr[5], antl[5], bntl[5];
    double ant1[5], bnt1[5], ant2[5], bnt2[5];

    for (int i = 1; i <= 4; i++)
    {
        antr[i] = z[i][1] * ap1tr + z[i][2] * ap2tr + z[i][3] * ap3tr + z[i][4] * ap4tr;
        bntr[i] = z[i][1] * bp1tr + z[i][2] * bp2tr + z[i][3] * bp3tr + z[i][4] * bp4tr;
        antl[i] = z[i][1] * ap1tl + z[i][2] * ap2tl + z[i][3] * ap3tl + z[i][4] * ap4tl;
        bntl[i] = z[i][1] * bp1tl + z[i][2] * bp2tl + z[i][3] * bp3tl + z[i][4] * bp4tl;
    }

    for (int i = 1; i <= 4; i++)
    {
        ant1[i] = ct * antl[i] + st * antr[i];
        bnt1[i] = ct * bntl[i] + st * bntr[i];
        ant2[i] = -st * antl[i] + ct * antr[i];
        bnt2[i] = -st * bntl[i] + ct * bntr[i];
    }

    double ax1tbl = 0.0;
    double bx1tbl = g2;
    double ax1tbr = 0.0;
    double bx1tbr = 0.0;
    double ax2tbl = -yt;
    double bx2tbl = 0.0;
    double ax2tbr = 0.0;
    double bx2tbr = -yb;

    double actbl[3], actbr[3], bctbl[3], bctbr[3];
    double actb1[3], actb2[3], bctb1[3], bctb2[3];
    for (int i = 1; i <= 2; i++)
    {
        actbl[i] = v[i][1] * ax1tbl + v[i][2] * ax2tbl;
        bctbl[i] = u[i][1] * bx1tbl + u[i][2] * bx2tbl;
        actbr[i] = v[i][1] * ax1tbr + v[i][2] * ax2tbr;
        bctbr[i] = u[i][1] * bx1tbr + u[i][2] * bx2tbr;
    }

    for (int i = 1; i <= 2; i++)
    {
        actb1[i] = cb * actbl[i] + sb * actbr[i];
        bctb1[i] = cb * bctbl[i] + sb * bctbr[i];
        actb2[i] = -sb * actbl[i] + cb * actbr[i];
        bctb2[i] = -sb * bctbl[i] + cb * bctbr[i];
    }

    //
    // Second last line of D18 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------
        
    double dnino = 0.0;
    for (int i = 1; i <= 4; i++)
    {
        dnino = dnino + .5 * ((pow(ant1[i], 2) + pow(bnt1[i], 2)) * creal(su_b1(pow(mt, 2), dxmn[i], mst1)) + 2 * ant1[i] * bnt1[i] * dxmn[i] / mt * creal(su_b0(pow(mt, 2), dxmn[i], mst1)));
        dnino = dnino + .5 * ((pow(ant2[i], 2) + pow(bnt2[i], 2)) * creal(su_b1(pow(mt, 2), dxmn[i], mst2)) + 2 * ant2[i] * bnt2[i] * dxmn[i] / mt * creal(su_b0(pow(mt, 2), dxmn[i], mst2)));
    }


    //
    // Last line of D18 (https://arxiv.org/pdf/hep-ph/9606211)
    // ------------------------------------------------------------------------------------------

    double mc[3];
    mc[1] = dmc1;
    mc[2] = dmc2;
    double dcino = 0.0;
    for (int i = 1; i <= 2; i++)
    {
        dcino = dcino + .5 * ((pow(actb1[i], 2) + pow(bctb1[i], 2)) * creal(su_b1(pow(mt, 2), mc[i], msb1)) + 2 * actb1[i] * bctb1[i] * mc[i] / mt * creal(su_b0(pow(mt, 2), mc[i], msb1)));
        dcino = dcino + .5 * ((pow(actb2[i], 2) + pow(bctb2[i], 2)) * creal(su_b1(pow(mt, 2), mc[i], msb2)) + 2 * actb2[i] * bctb2[i] * mc[i] / mt * creal(su_b0(pow(mt, 2), mc[i], msb2)));
    }


    //
    // Pure QCD correction (including logs):
    //-----------------------------------------------------------------------------------------------------------------------------------

    double mtlog = log(pow((rmt / mz), 2));
    double delmt = alphas / pi * (5.0 / 3 - mtlog) + pow(alphas, 2) * (0.5380 - 0.1815 * mtlog + 0.038 * pow(mtlog, 2));

    //
    // SUSY contributions added:
    //-----------------------------------------------------------------------------------------------------------------------------------

    *delmtop = -delmt + (dqcd + dyuk + dgauge + dnino + dcino) / (16 * pow(pi, 2));
    scale = scalsave;

}

//
// Calculates the radiative corrections to the gaugino and mu mass parameters. input parameters at ewsb scale are:
// 	ml1,mq1,mq3,mu3,md3: sfermion mass parameters of 1st and 3d generations
// 	ma,tb: pseudoscalar higgs boson mass and tan(beta)
// 	yt,yb: top and bottom yukawa couplings,
// 	m1,m2,mu: bare gaugino and higgs mass parameters.
// The outputs are the radiative corrections  rcm1,rcm2,rcmu to m1,m2, mu
// -----------------------------------------------------------------------------------------------------------------------

//
//	Precision corrections in the minimal supersymmetric standard model
//		Damien M. Pierce,
//		Jonathan A. Bagger,
//		Konstantin T. Matchev
//		Ren-Jie Zhang
//	Ref : arXiv:hep-ph/9606211v3
// -----------------------------------------------------------------------------------------------------------------------

void su_radcino(double ml1, double mq1, double mq3, double mu3, double md3, double ma, double yt, double yb, double m1, double m2, double mu, double tb, double *rcm1, double *rcm2, double *rcmu)
{

    if (kpole == 1)
        scale = sqrt(msttr1 * msttr2);

    double cw = g2ewsb / sqrt(pow(g1ewsb, 2) + pow(g2ewsb, 2));
    double sw = sqrt(1.0 - pow(cw, 2));
    double b = atan(tb);
    double ep = 1.0e-5;
    double amu = fabs(mu);
    double alphewsb = pow((g2ewsb * sw), 2) / (4 * pi);

    //
    // The dominant correction to M1 comes from quark/squark, chargino/charged-Higgs and neutralino/neutral-Higgs loops
    // Eq (25)
    // --------------------------------------------------------------------------------------------------------------------

    double rm1 = 0.0;
    rm1 += 11.0 * creal(su_b1(pow(m1, 2), ep, mq1));
    rm1 += +9.0 * creal(su_b1(pow(m1, 2), ep, ml1));
    rm1 += +mu / m1 * sin(2 * b) * (creal(su_b0(pow(m1, 2), amu, ma)) - creal(su_b0(pow(m1, 2), amu, mz)));
    rm1 += +creal(su_b1(pow(m1, 2), amu, ma));
    rm1 += +creal(su_b1(pow(m1, 2), amu, mz));

    *rcm1 = -alphewsb / (4.0 * pi * pow(cw, 2)) * rm1 * m1;

    //
    // The corrections to M2 from quark/squark and Higgs loops and the gauge boson loops
    // 	Eq (27), Eq (29)
    // --------------------------------------------------------------------------------------------------------------------

    double rm2 = 0.0;
    rm2 += 9.0 * creal(su_b1(pow(m2, 2), ep, mq1));
    rm2 += +3.0 * creal(su_b1(pow(m2, 2), ep, ml1));
    rm2 += +mu / m2 * sin(2 * b) * (creal(su_b0(pow(m2, 2), amu, ma)) - creal(su_b0(pow(m2, 2), amu, mz)));
    rm2 += +creal(su_b1(pow(m2, 2), amu, ma));
    rm2 += +creal(su_b1(pow(m2, 2), amu, mz));
    rm2 += -4.0 * (2.0 * creal(su_b0(pow(m2, 2), m2, mw)) - creal(su_b1(pow(m2, 2), m2, mw)));

    *rcm2 = -alphewsb / (4.0 * pi * pow(sw, 2)) * rm2 * m2;

    //
    //	The corrections to  are obtained in a similar manner
    // 	Eq (31)
    // --------------------------------------------------------------------------------------------------------------------

    double rmu1 = 0.0;
    rmu1 += (pow(ybewsb, 2) + pow(ytewsb, 2)) * creal(su_b1(pow(amu, 2), ep, mq3));
    rmu1 += +pow(ytewsb, 2) * creal(su_b1(pow(amu, 2), ep, mu3));
    rmu1 += +pow(ybewsb, 2) * creal(su_b1(pow(amu, 2), ep, md3));

    double rmu2 = 0.0;
    rmu2 += creal(su_b1(pow(amu, 2), m2, ma));
    rmu2 += +creal(su_b1(pow(amu, 2), m2, mz));
    rmu2 += +2.0 * creal(su_b1(pow(amu, 2), amu, mz));
    rmu2 += -4.0 * creal(su_b0(pow(amu, 2), amu, mz));

    *rcmu = (-3.0 / (32 * pow(pi, 2)) * rmu1 - 3 * alphewsb / (16 * pi * pow(sw, 2)) * rmu2) * mu;
}

//
// Eq 16 from (https://arxiv.org/pdf/hep-ph/9606211)
// --------------------------------------------------------------------------
double su_taumscr(double tgbeta, double mu, double m2, double mnstau)
{

    double scalsave = scale;
    scale = mz;
    mtau = 1.7771;
    pi = 4 * atan(1.0);

    double cinostau = g22 / 16 / pow(pi, 2) * mu * m2 * tgbeta / (pow(mu, 2) - pow(m2, 2)) * (creal(su_b0(pow(mtau, 2), m2, mnstau)) - creal(su_b0(pow(mtau, 2), mu, mnstau)));

    scale = scalsave;
    return cinostau;
}


//
// Calculates the qcd (standard+susy) correction to Squark (except stop) masses.
//
//	The input are: the strong coupling constant alphas,
//						the gluino and tree-level squark masses
//	The output is: the correction to the squark mass dmsquark.
//
//	p.s.squark mixing and yukawa's are neglected.
// ---------------------------------------------------------------------------------------------------------------------------

//
//	Precision Corrections in the minimal supersymmetric standard model
//		Damien M. Pierce, Jonathan A. Bagger, Konstaintin T. Matchev and Ren-Jie Zhang
//	arXiv:hep-ph/9606211v3        ..  Eq(33) & Eq(34)
// ---------------------------------------------------------------------------------------------------------------------------

void su_sqcr(double alphas, double mgluino, double msquark, double *dmsquark)
{

    double scale;
    if (kpole == 1)
        scale = sqrt(msttr1 * msttr2);

    double corr;
    double x = pow(mgluino / msquark, 2);
    double corr2 = 2 * alphas / pi / 3 * (1.0 + 3 * x + pow((x - 1.0), 2) * log(fabs(x - 1.0)) - pow(x, 2) * log(x) + 4 * x * log(scale / msquark));
    if (corr2 > -1.0)
    {
        corr = sqrt(1.0 + corr2) - 1.0;
        *dmsquark = msquark * corr;
    }
    else
    {
        *dmsquark = 0.0;
        tachsqrc = -1.0; // Tachionic state
    }
}


//
// Calculates the SUSY radiative corrections to the bottom mass including the susy QCD corrections (the standard ones are
//	calculated with runm) and the dominant electroweak corrections due to the yukawa couplings.
//
// The input are respectively:
//		The strong coupling constant, pole b mass,
// 	The running top and bottom masses, the top yukawa coupling, tan(beta),
// The su[2] gaugino mass, the gluino mass, the 3d generation squark mass  terms, the 3d generation trilinear couplings
//	and the parameter mu.
// The output delmtop is the susy radiative correction to the bottom mass.
// These corrections are then re-summed in the main routine.
//-------------------------------------------------------------------------------------------------------------------------------------

//
//	Referance:
//	 Precision corrections in the minimal supersymmetric standard model
//	 D.M.Pierce, J.A.Bagger, K.T.Matchev and R.J.Zhang
//  arXiv:hep-ph/9606211
//-------------------------------------------------------------------------------------------------------------------------------------

#define b11(x) (.5 * (.5 + 1.0 / (1.0 - x) + log(x) / pow((1.0 - x), 2)))
#define b12(x) (.5 * (.5 + 1.0 / (1.0 - x) + log(x) / pow((1.0 - x), 2) - log(x)))

void su_bmsusycr(double alphas, double mb, double rmt, double rmb, double yt, double tbeta, double m2, double mgluino, double mql, double mur, double mdr, double at, double ab, double mu, double *delmb)
{

    double scalsave = scale;
    scale = mz;
    double msb1 = msb1sfbp;
    double msb2 = msb2sfbp;
    double mst1 = mst1sfbp;
    double mst2 = mst2sfbp;

    double b = atan(tbeta);
    double ct2 = pow(cos(thet), 2);
    double st2 = pow(sin(thet), 2);
    double x1 = pow((msb1 / mgluino), 2);
    double x2 = pow((msb2 / mgluino), 2);
    double mm1 = (msb1> mgluino) ? msb1: mgluino;
    double mm2 = (msb2> mgluino) ?msb2: mgluino;
    double r1, r2;

    if (x1 >= 1.0)
        r1 = b11(x1) - log(pow(mm1, 2) / pow(scale, 2)) / 2;
    else
        r1 = b12(x1) - log(pow(mm1, 2) / pow(scale, 2)) / 2;

    if (x2 >= 1.0)
        r2 = b11(x2) - log(pow(mm2, 2) / pow(scale, 2)) / 2;
    else
        r2 = b12(x2) - log(pow(mm2, 2) / pow(scale, 2)) / 2;

    double ginosq = -alphas / pi / 3 * (r1 + r2 - mgluino / rmb * sin(2 * theb) * (creal(su_b0(pow(mb, 2), mgluino, msb1)) - creal(su_b0(pow(mb, 2), mgluino, msb2))));

    double temp = 1.0;
    temp *= -g22 / 16 / pow(pi, 2) * mu * m2 * tbeta / (pow(mu, 2) - pow(m2, 2));
    temp *= (ct2 * creal(su_b0(pow(mb, 2), m2, mst1)) + st2 * creal(su_b0(pow(mb, 2), m2, mst2)) - ct2 * creal(su_b0(pow(mb, 2), mu, mst1)) - st2 * creal(su_b0(pow(mb, 2), mu, mst2)));

    double cinost = -pow(yt, 2) / pow(pi, 2) / 16 * mu * tbeta / 2.0 * sin(2 * thet) / rmt * (creal(su_b0(pow(mb, 2), mu, mst1)) - creal(su_b0(pow(mb, 2), mu, mst2))) + temp;

    *delmb = ginosq + cinost;
    scale = scalsave;
}

#undef b11
#undef b12


//
// Calculates leading one-loop susy delta_rho contributions of 3rd gen sfermions (plus leading two-loop qcd contributions)
//  input: mt, gmst[2], gmsb[2],gmstau[2],msn: top,stop,sbottom,
//			  stau, stau neutrino masses and stop, sbottom, stau mixing angles
//  output: drho = rho-1
//-------------------------------------------------------------------------------------------------------------------------------------

//
//   PRECISION CORRECTIONS IN THE MINIMAL SUPERSYMMETRIC STANDARD MODEL
//   Equation C4. (https://arxiv.org/abs/hep-ph/9606211)
// -----------------------------------------------------------------------------------------------------

//
// We may also check 
// Leading QCD Corrections to Scalar Quark Contributions to Electroweak Precision Observables
// A. Djouadi, P. Gambino, S. Heinemeyer, W. Hollik, C. Junger, G. Weiglein
// (https://arxiv.org/pdf/hep-ph/9710438) Check Eq. 20
// ------------------------------------------------------------------------------------------------------


#define su_fr(x,y) x+y-2*x*y/(x-y)*log(x/y)

void su_delrho(double mt,double gmst[],double gmsb[],double gmstau[],double msn,double thetat,double thetab,double thel,double drho)
{
    double ct=cos(thetat);
    double st=sin(thetat);
    double cb=cos(thetab);
    double sb=sin(thetab);
    double ctau =cos(thel);
    double stau =sin(thel);
    double cta2=pow(ctau,2);
    double sta2=pow(stau,2);
    double ct2=pow(ct,2);
    double st2=pow(st,2);
    double cb2=pow(cb,2);
    double sb2=pow(sb,2);
    double mt1=pow(gmst[1],2);
    double mt2=pow(gmst[2],2);
    double mb1=pow(gmsb[1],2);
    double mb2=pow(gmsb[2],2);
    double mta1=pow(gmstau[1],2);
    double mta2=pow(gmstau[2],2);

    double drhotb= (ct2*(cb2*su_fr(mt1,mb1)+sb2*su_fr(mt1,mb2)) + st2*(cb2*su_fr(mt2,mb1)+sb2*su_fr(mt2,mb2)) - ct2*st2*su_fr(mt1,mt2)-cb2*sb2*su_fr(mb1,mb2));
    double drhotau= -cta2*sta2*su_fr(mta1,mta2)+cta2*su_fr(mta1,msn*msn) + sta2*su_fr(mta2,msn*msn);
    drho = 3*drhotb*(1.0 +2*0.12/3/pi*(1.0+pow(pi,2)/3))+drhotau;
    drho = gf/(8* pow(pi,2)* sqrt(2.0))*drho;
    }

#undef su_fr
