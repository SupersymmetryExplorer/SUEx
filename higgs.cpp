#include <stdio.h>
#include <math.h>
#include "variables.h"
#include "loop.h"
#include "functions.h"

//
// The following routine is for the evaluation of the Higgs boson mass
// -------------------------------------------------------------------------------------------------------------

//
// Calculates the mssm higgs bosons masses and the angle alpha including radiative corrections for a given
//	input value of the parameter tan(beta). The other input parameters (soft-susy breaking parameters, sparticle
//  masses and mixing angles, sm parameters, are called via common blocks.
//  It returns the pole masses of the
//				1. cp-odd (ama),
//				2. lighter cp-even (aml),
//				3. heavier cp-even (amh),
//				4. charged higgs boson (amch) as well as
//  the running cp-odd (amar) higgs masses.
//
//  it gives also the couplings of the angle beta at the EWSB scale, the mixing
//  alpha and the higgs boson couplings to standard particles in:
//        common/coup_hcoup/gat,gab,glt,glb,ght,ghb,glvv,ghvv,b,a
//  it returns also the couplings of the higgs bosons to sfermions
//        common/su_cplhsf/gcen,gctb,glee,gltt,glbb,ghee,ghtt,ghbb
//      .                   gatt,gabb,gaee
//  and the higgs couplings to charginos and neutralinos:
//         common/su_cplhino/ac1,ac2,ac3,an1,an2,an3,acnl,acnr
//  for the radiative correction of higgs masses, there is imodel=0
//  option where the calculation is made in an approximation based on the
//  of work heinemeyer, hollik, weiglein (hep-ph/0002213), which is fast
//  but approximate, or
//  ichoice[10](=imodel)=1:  full one-loop higgs masses  a la pbmz
//                      =2:  full one-loop pbmz + two-loop bdsz corrections
// -------------------------------------------------------------------------------------------------------------

//
// FeynHiggsFast: a program for a fast calculation of masses and mixing angles in the Higgs Sector of the MSSM
// 	S. Heinemeyer, W.Hollik, G.Heiglein
//		arXiv:hep-ph/0002213v1
// -------------------------------------------------------------------------------------------------------------

int su_susycp(double tgbet)
{
    double aml = aaml;
    double amh = aamh;
    double amch = aamch;
    double ama = aama;

    double qqn[5][5], ssn[5][5], qqc[3][3], ssc[3][3];
    double vlocal = 1.0 / sqrt(sqrt(2.0) * gf);
    tbeta = vuewsb / vdewsb;
    double bet = atan(tbeta);
    double amt = mtpole;
    double am3 = m3;
    double am2 = m2;
    double amz = mz;
    double mst[3], msb[3], msllocal[3];
    double b = beta;

    if (b == 0.0)
        b = bet;

    double sb = sin(bet);
    double cb = cos(bet);
    double amar = ama;
    double marsave = ama;

    // nb at this stage ama is in fact running ma

    double mt = runm(amt, 6);
    double mb = runm(amt, 5);
    double als = falphas(amt, 2);

    sw2 = pow(g1ewsb, 2) / (pow(g1ewsb, 2) + pow(g2ewsb, 2));

    //
    // FeynHiggsFast: a program for a fast calculation of masses and mixing angles in the Higgs Sector of the MSSM
    // 	S. Heinemeyer, W.Hollik, G.Heiglein
    // arXiv:hep-ph/0002213
    // ----------------------------------------------------------------------------------------------------------

    double vev2, madr2;
    double amsq = msq;
    double amur = mtr;
    double amu = mu;
    int imodel = ihflag;
    double amdr = mbr;

    if (imodel == 0)
    {
        vev2 = 2.0 * (pow(vuewsb, 2) + pow(vdewsb, 2));
        madr2 = pow(marsave, 2);
    }

    double amglu = am3;
    double ams2 = sqrt(pow(amsq, 2) * pow(amur, 2) + pow(mt, 2) * (pow(amsq, 2) + pow(amur, 2)) + pow(mt, 4));      // Eq(16)
    double xlam = 1.0 / 8.0 - sw2 / 3.0 + 4 * pow(sw2, 2) / 9.0;                                                    // Eq(13)
    double xt = au - amu / tbeta;                                                                                   // Eq(2)
    double xr = pow(mt, 2) / ams2;                                                                                  // Xr = mt^2 / Ms^2
    double xfac = gf * sqrt(2.0) / pow(pi, 2);                                                                      // Used in different eqn
    double s11 = xfac * pow(amz, 4) * xlam * pow(cb, 2) * log(xr);                                                  // Eq(10)
    double s12 = -xfac * pow(amz, 2) / tbeta * (-3 * pow(mt, 2) / 8.0 + pow(amz, 2) * xlam * pow(sb, 2)) * log(xr); // Eq(11)


    double temp = 0.0;
    temp += -2 * pow(amz, 2) / pow(mt, 2) + 11 * pow(amz, 4) / 10.0 / pow(mt, 4);
    temp += +(12.e0 - 6 * pow(amz, 2) / pow(mt, 2) * pow(sb, 2) + 8 * pow(amz, 4) / pow(mt, 4) * xlam * pow(sb, 4)) * log(xr);
    temp += +pow(xt, 2) / ams2 * (-12.0 + 4.0 * pow(amz, 2) / pow(mt, 2) + 6.0 * xr);
    temp += +pow(xt, 4) / pow(ams2, 2) * (1.0 - 4 * xr + 3 * pow(xr, 2));
    temp += +pow(xt, 6) / pow(ams2, 3) * (3 * xr / 5.0 - 12 * pow(xr, 2) / 5.0 + 2 * pow(xr, 3));
    temp += +pow(xt, 8) / pow(ams2, 4) * (3 * pow(xr, 2) / 7.0 - 12 * pow(xr, 3) / 7.0 + 3 * pow(xr, 4) / 2.0);
    double s22one = xfac * pow(mt, 4) / 8.0 / pow(sb, 2) * temp; // Eq(12)

    temp = 0.0;
    temp += 4.0 + 3 * pow(log(xr), 2) + 2 * log(xr) - 6 * xt / sqrt(ams2);
    temp += -pow(xt, 2) / ams2 * (3 * log(xr) + 8.0) + 17 * pow(xt, 4) / 12.0 / pow(ams2, 2);
    double s22qcd = xfac * als / pi * pow(mt, 4) / pow(sb, 2) * temp; // Eq(14)


    double s22ew = (pow(log(xr), 2) + 2 * pow(xt, 2) / ams2 * log(xr) + pow(xt, 4) / 6.e0 / pow(ams2, 2) * log(xr));
    s22ew = -9 * pow(xfac, 2) / 32.0 / pow(sb, 2) * pow(mt, 6) * s22ew;

    double s22 = s22one + s22qcd + s22ew;

    double xm11 = pow(ama, 2) * pow(sb, 2) + pow(amz, 2) * pow(cb, 2) - s11; // Eq(5) and Eq(6)
    double xm12 = -(pow(ama, 2) + pow(amz, 2)) * sb * cb - s12;              // Eq(5) and Eq(6)
    double xm22 = pow(ama, 2) * pow(cb, 2) + pow(amz, 2) * pow(sb, 2) - s22; // Eq(5) and Eq(6)

    double xml2 = 0.50 * (xm11 + xm22 - sqrt(pow((xm11 - xm22), 2) + 4 * pow(xm12, 2))); // Calculating the eigen values Eq(6) mH^2
    double xmh2 = 0.50 * (xm11 + xm22 + sqrt(pow((xm11 - xm22), 2) + 4 * pow(xm12, 2))); // Calculating the eigen values Eq(6) mh^2

    double amlr = sqrt(xml2); // Eq(6) mh
    double amhr = sqrt(xmh2); // Eq(6) mH


    ama = amar;
    a = atan(xm12 / (pow(amz, 2) * pow(cb, 2) + pow(ama, 2) * pow(sb, 2) - s11 - pow(aml, 2))); // Alpha_effective Eq(7)
 
    double sa = sin(a);
    double ca = cos(a);

    if (ca == 0)
        a = pi / 2;
    else
        a = atan(sa / ca);

    if (ca < 0.0)
    {
        if (sa < 0.0)
            a = a - pi;
        else
            a = a + pi;
    }

    sa = sin(a);
    ca = cos(a);

    // =====================================================================
    //  Now calculate the higgs boson coupling to sfermions and gauginos:
    // =====================================================================

    double sbma = sb * ca - cb * sa;
    double cbma = cb * ca + sb * sa;
    double sbpa = sb * ca + cb * sa;
    double cbpa = cb * ca - sb * sa;
    double mstl2 = pow(amsq, 2) + (0.50 - 2.0 / 3.0 * sw2) * pow(amz, 2) * cos(2.0 * b);
    double mstr2 = pow(amur, 2) + 2.0 / 3.0 * sw2 * pow(amz, 2) * cos(2.0 * b);
    double mlrt = au - amu / tgbet;
    double delt = pow((mstl2 - mstr2), 2) + 4 * pow(mt, 2) * pow(mlrt, 2);
    double mst12 = pow(mt, 2) + 0.50 * (mstl2 + mstr2 - sqrt(delt));
    double mst22 = pow(mt, 2) + 0.50 * (mstl2 + mstr2 + sqrt(delt));

    if (mst12 < 0.0 && imodel == 0)
        return 0;

    mst[1] = sqrt(mst12);
    mst[2] = sqrt(mst22);

    double thet;
    if (mstl2 == mstr2)
        thet = pi / 4;
    else
    {
        thet = 0.50 * atan(2.0 * mt * mlrt / (mstl2 - mstr2));
        if (mstl2 > mstr2)
            thet = thet + pi / 2;
    }

    double cst = cos(thet);
    double sst = sin(thet);

    //===== sbottom masses
    double msbl2 = pow(amsq, 2) + (-0.50 + 1.0 / 3.0 * sw2) * pow(amz, 2) * cos(2.0 * b);
    double msbr2 = pow(amdr, 2) - 1.0 / 3.0 * sw2 * pow(amz, 2) * cos(2.0 * b);
    double mlrb = ad - amu * tgbet;
    double delb = pow((msbl2 - msbr2), 2) + 4 * pow(mb, 2) * pow(mlrb, 2);
    double msb12 = pow(mb, 2) + 0.50 * (msbl2 + msbr2 - sqrt(delb));
    double msb22 = pow(mb, 2) + 0.50 * (msbl2 + msbr2 + sqrt(delb));

    if (msb12 < 0.0 && imodel == 0)
        return 0;

    msb[1] = sqrt(msb12);
    msb[2] = sqrt(msb22);
    double theb;
    if (msbl2 == msbr2)
        theb = pi / 4;
    else
    {
        theb = 0.50 * atan(2.0 * mb * mlrb / (msbl2 - msbr2));
        if (msbl2 > msbr2)
            theb = theb + pi / 2;
    }

    double csb = cos(theb);
    double ssb = sin(theb);

    //===== stau masses

    double amel = msl;
    double amer = mtaur;
    double amtau = mtau;

    double msel2 = pow(amel, 2) + (-0.50 + sw2) * pow(amz, 2) * cos(2.0 * b);
    double mser2 = pow(amer, 2) - sw2 * pow(amz, 2) * cos(2.0 * b);
    double msn2 = pow(amel, 2) + 0.50 * pow(amz, 2) * cos(2.0 * b);
    double mlre = al - amu * tgbet;
    double dele = pow((msel2 - mser2), 2) + 4 * pow(amtau, 2) * pow(mlre, 2);
    double mse12 = pow(amtau, 2) + 0.50 * (msel2 + mser2 - sqrt(dele));
    double mse22 = pow(amtau, 2) + 0.50 * (msel2 + mser2 + sqrt(dele));

    if (mse12 < 0.0 && imodel == 0)
        return 0;

    double msn;
    double ml11loop, ml22loop, ml12loop, mh11loop, mh22loop, mh12loop;
    double l12loop;
    double h12loop;

    msllocal[1] = sqrt(mse12);
    msllocal[2] = sqrt(mse22);
    msn = sqrt(msn2);
    double thel;
    if (msel2 == mser2)
        thel = pi / 4;
    else
    {
        thel = 0.50 * atan(2.0 * amtau * mlre / (msel2 - mser2));
        if (msel2 > mser2)
            thel = thel + pi / 2;
    }

    double csl = cos(thel);
    double ssl = sin(thel);

    //===== light higgs couplings to sfermions
    glt = ca / sb;
    glb = -sa / cb;

    gltt[1][1] = -sbpa * (0.50 * pow(cst, 2) - 2.0 / 3.0 * sw2 * cos(2 * thet)) + pow(mt, 2) / pow(amz, 2) * glt + mt * sst * cst / pow(amz, 2) * (au * glt + amu * ght);
    gltt[2][2] = -sbpa * (0.50 * pow(sst, 2) + 2.0 / 3.0 * sw2 * cos(2 * thet)) + pow(mt, 2) / pow(amz, 2) * glt - mt * sst * cst / pow(amz, 2) * (au * glt + amu * ght);
    gltt[1][2] = -2 * sbpa * sst * cst * (2.0 / 3.0 * sw2 - 0.25) + mt * cos(2 * thet) / 2.0 / pow(amz, 2) * (au * glt + amu * ght);
    gltt[2][1] = -2 * sbpa * sst * cst * (2.0 / 3.0 * sw2 - 0.25) + mt * cos(2 * thet) / 2.0 / pow(amz, 2) * (au * glt + amu * ght);

    glbb[1][1] = -sbpa * (-0.50 * pow(csb, 2) + 1.0 / 3.0 * sw2 * cos(2 * theb)) + pow(mb, 2) / pow(amz, 2) * glb + mb * ssb * csb / pow(amz, 2) * (ad * glb - amu * ghb);
    glbb[2][2] = -sbpa * (-0.50 * pow(ssb, 2) - 1.0 / 3.0 * sw2 * cos(2 * theb)) + pow(mb, 2) / pow(amz, 2) * glb - mb * ssb * csb / pow(amz, 2) * (ad * glb - amu * ghb);
    glbb[1][2] = -2 * sbpa * ssb * csb * (-1.0 / 3.0 * sw2 + 0.250) + mb * cos(2 * theb) / 2.0 / pow(amz, 2) * (ad * glb - amu * ghb);
    glbb[2][1] = -2 * sbpa * ssb * csb * (-1.0 / 3.0 * sw2 + 0.250) + mb * cos(2 * theb) / 2.0 / pow(amz, 2) * (ad * glb - amu * ghb);

    glee[1][1] = -sbpa * (-0.50 * pow(csl, 2) + sw2 * cos(2 * thel)) + pow(amtau, 2) / pow(amz, 2) * glb + amtau * ssl * csl / pow(amz, 2) * (al * glb - amu * ghb);
    glee[2][2] = -sbpa * (-0.50 * pow(ssl, 2) - sw2 * cos(2 * thel)) + pow(amtau, 2) / pow(amz, 2) * glb - amtau * ssl * csl / pow(amz, 2) * (al * glb - amu * ghb);
    glee[1][2] = -2 * sbpa * ssl * csl * (-sw2 + 0.250) + amtau * cos(2 * thel) / 2.0 / pow(amz, 2) * (al * glb - amu * ghb);
    glee[2][1] = -2 * sbpa * ssl * csl * (-sw2 + 0.250) + amtau * cos(2 * thel) / 2.0 / pow(amz, 2) * (al * glb - amu * ghb);

    //===== heavy higgs couplings to sfermions
    ght = sa / sb;
    ghb = ca / cb;

    ghtt[1][1] = cbpa * (0.50 * pow(cst, 2) - 2.0 / 3.0 * sw2 * cos(2 * thet)) + pow(mt, 2) / pow(amz, 2) * ght + mt * sst * cst / pow(amz, 2) * (au * ght - amu * glt);
    ghtt[2][2] = cbpa * (0.50 * pow(sst, 2) + 2.0 / 3.0 * sw2 * cos(2 * thet)) + pow(mt, 2) / pow(amz, 2) * ght - mt * sst * cst / pow(amz, 2) * (au * ght - amu * glt);
    ghtt[1][2] = 2 * cbpa * sst * cst * (2.0 / 3.0 * sw2 - 0.250) + mt * cos(2 * thet) / 2.0 / pow(amz, 2) * (au * ght - amu * glt);
    ghtt[2][1] = 2 * cbpa * sst * cst * (2.0 / 3.0 * sw2 - 0.250) + mt * cos(2 * thet) / 2.0 / pow(amz, 2) * (au * ght - amu * glt);

    ghbb[1][1] = cbpa * (-0.50 * pow(csb, 2) + 1.0 / 3.0 * sw2 * cos(2 * theb)) + pow(mb, 2) / pow(amz, 2) * ghb + mb * ssb * csb / pow(amz, 2) * (ad * ghb + amu * glb);
    ghbb[2][2] = cbpa * (-0.50 * pow(ssb, 2) - 1.0 / 3.0 * sw2 * cos(2 * theb)) + pow(mb, 2) / pow(amz, 2) * ghb - mb * ssb * csb / pow(amz, 2) * (ad * ghb + amu * glb);
    ghbb[1][2] = 2 * cbpa * ssb * csb * (-1.0 / 3.0 * sw2 + 0.250) + mb * cos(2 * theb) / 2.0 / pow(amz, 2) * (ad * ghb + amu * glb);
    ghbb[2][1] = 2 * cbpa * ssb * csb * (-1.0 / 3.0 * sw2 + 0.250) + mb * cos(2 * theb) / 2.0 / pow(amz, 2) * (ad * ghb + amu * glb);

    ghee[1][1] = cbpa * (-0.50 * pow(csl, 2) + sw2 * cos(2 * thel)) + pow(amtau, 2) / pow(amz, 2) * ghb + amtau * ssl * csl / pow(amz, 2) * (al * ghb + amu * glb);
    ghee[2][2] = cbpa * (-0.50 * pow(ssb, 2) - sw2 * cos(2 * thel)) + pow(amtau, 2) / pow(amz, 2) * ghb - amtau * ssl * csl / pow(amz, 2) * (al * ghb + amu * glb);
    ghee[1][2] = 2 * cbpa * ssl * csl * (-sw2 + 0.250) + amtau * cos(2 * thel) / 2.0 / pow(amz, 2) * (al * ghb + amu * glb);
    ghee[2][1] = 2 * cbpa * ssl * csl * (-sw2 + 0.250) + amtau * cos(2 * thel) / 2.0 / pow(amz, 2) * (al * ghb + amu * glb);

    //===== pseudoscalar higgs couplings to sfermions
    gat = 1.0 / tgbet;
    gab = tgbet;
    gatt = -mt / 2.0 / pow(amz, 2) * (amu + au * gat);
    gabb = -mb / 2.0 / pow(amz, 2) * (amu + ad * gab);
    gaee = -amtau / 2.0 / pow(amz, 2) * (amu + al * gab);

    double amw = mw;
    double mlpole, mhpole, mchpole;

    //===== charged higgs couplings sfermions
    double cll3 = (pow(amw, 2) * sin(2 * b) - pow(mt, 2) * gat - pow(mb, 2) * gab) / sqrt(2.0) / pow(amw, 2);
    double crr3 = -mt * mb * (gat + gab) / sqrt(2.0) / pow(amw, 2);
    double clr3 = -mb * (amu + ad * gab) / sqrt(2.0) / pow(amw, 2);
    double crl3 = -mt * (amu + au * gat) / sqrt(2.0) / pow(amw, 2);
    gctb[1][1] = +cst * csb * cll3 + sst * ssb * crr3 + cst * ssb * clr3 + sst * csb * crl3;
    gctb[1][2] = -cst * ssb * cll3 + sst * csb * crr3 + cst * csb * clr3 - sst * ssb * crl3;
    gctb[2][1] = -sst * csb * cll3 + cst * ssb * crr3 - sst * ssb * clr3 + cst * csb * crl3;
    gctb[2][2] = +sst * ssb * cll3 + cst * csb * crr3 - sst * csb * clr3 - cst * ssb * crl3;
    double cll1 = (pow(amw, 2) * sin(2 * b) - pow(amtau, 2) * gab) / sqrt(2.0) / pow(amw, 2);
    double clr1 = -amtau * (amu + al * gab) / sqrt(2.0) / pow(amw, 2);
    gcen[1][1] = csl * cll1 + ssl * clr1;
    gcen[1][2] = -ssl * cll1 + csl * clr1;
    gcen[2][1] = 0.0;
    gcen[2][2] = 0.0;

    double p2l;
    double mhd2 = mhd2;
    double mhu2 = mhu2;


    //=====  neutral higgs couplings to neutralinos
    double tanw = sqrt(sw2) / sqrt(1.0 - sw2);

    for (int i = 1; i <= 4; i++)
        for (int j = 1; j <= 4; j++)
        {
            qqn[i][j] = 1.0 / 2.0 * (z[i][3] * (z[j][2] - tanw * z[j][1]) + z[j][3] * (z[i][2] - tanw * z[i][1]));
            ssn[i][j] = 1.0 / 2.0 * (z[i][4] * (z[j][2] - tanw * z[j][1]) + z[j][4] * (z[i][2] - tanw * z[i][1]));
        }

    for (int i = 1; i <= 4; i++)
        for (int j = 1; j <= 4; j++)
        {
            an1[i][j] = qqn[i][j] * cos(a) - ssn[i][j] * sin(a);
            an2[i][j] = -qqn[i][j] * sin(a) - ssn[i][j] * cos(a);
            an3[i][j] = qqn[i][j] * sin(bet) - ssn[i][j] * cos(bet);
        }

    //=====  neutral higgs couplings to charginos
    for (int i = 1; i <= 2; i++)
        for (int j = 1; j <= 2; j++)
        {
            qqc[i][j] = sqrt(1.0 / 2.0) * u[j][2] * v[i][1];
            ssc[i][j] = sqrt(1.0 / 2.0) * u[j][1] * v[i][2];
        }


    for (int i = 1; i <= 2; i++)
        for (int j = 1; j <= 2; j++)
        {
            ac1[i][j] = qqc[i][j] * cos(a) + ssc[i][j] * sin(a);
            ac2[i][j] = -qqc[i][j] * sin(a) + ssc[i][j] * cos(a);
            ac3[i][j] = qqc[i][j] * sin(bet) + ssc[i][j] * cos(bet);
        }


    for (int i = 1; i <= 2; i++)
        for (int j = 1; j <= 2; j++)
        {
            //=====  charged higgs couplings to charginos-neutralinos
            acnl[i][j] = cos(bet) * (z[j][4] * v[i][1] + (z[j][2] + z[j][1] * tanw) * v[i][2] / sqrt(2.0));
            acnr[i][j] = sin(bet) * (z[j][3] * u[i][1] - (z[j][2] + z[j][1] * tanw) * u[i][2] / sqrt(2.0));
        }

    double amchi2, amchr, tanba;
    double s11l, s12l, s22l;
    double tad1bl, tad2bl;
    double tad1l, tad2l;
    double tad1loop, tad2loop;

    amchr = 0.0;

    if (imodel >= 1)
    {
        printf("\n\nSorry for two loop or Accurate 1-loop for Higgs are not added in this calculatons");
        exit(1);
    }
    /*
   {
        // ============ gluino and heaviest chargino mass needed for subh ======
    amchi2 = pow(am2,2)+pow(amu,2)+2.0*pow(amw,2)+sqrt(pow((pow(am2,2)-pow(amu,2)),2) +4.0*pow(amw,4)*pow(cos(2.0*bet),2)+4.0*pow(amw,2)* (pow(am2,2)+pow(amu,2)+2.0*amu*am2*sin(2.0*bet) ) );
      double amchi = sqrt(0.50*amchi2);
      amglu = am3;
      printf("\n Hi : %e\n",amglu);

        //--use carena et al. for some things not included in the higgs routine:
      subh_hdec(ama,tgbet,amsq,amur,amdr,amt,au,ad,amu,amchi,amlr,amhr,amchr,sa,ca,tanba,amglu);

        //- now call the routine for the full one-loop or 2-loop calculation:
        //=======================================================================
      double q2 = pow(scale,2);

      tbeta = vuewsb/vdewsb;

      if(su_isnan(pizzp) || pow(amz,2)+pizzp <= 0.0)
      {
            // !!! protections added
            // non-pert or nan pb, uses tree-level values temporarily:
        pizzp = 0.0;
        if(irge == irgmax)
            inonpert=-1;
        }

      double rmz=sqrt(pow(amz,2)+pizzp);
      double rmw= rmz*sqrt(1.0-sw2);
      ama = marsave;

      if(kmaflag == 0 && irge >= 2)
      {
        if(madr2save > 0.0)
            ama =sqrt(madr2save);
        amar=ama;
        }

        printf("\n\nATan : %e %e %e %e %e %e %e %e",acnl[1][1],acnl[1][2],acnl[2][1],acnl[2][2],acnr[1][1],acnr[1][2],acnr[2][1],acnr[2][2]);exit(1);

      double amdelta=pow((pow(ama,2)+pow(rmz,2)),2)-4.0*pow(ama,2)*pow(rmz,2)*pow(cos(2.0*bet),2);
      double aml0=sqrt(0.50*(pow(ama,2)+pow(rmz,2)-sqrt(amdelta)));
      double amh0=sqrt(0.50*(pow(ama,2)+pow(rmz,2)+sqrt(amdelta)));
      double amch0=sqrt(pow(ama,2)+pow(rmw,2));

        // defining running higgs masses here for common (added):
      double marunp = ama;
      double mlrunp = aml0;
      double mhrunp = amh0;
      double mchrunp = amch0;
      double alfarunp =0.5*atan(tan(2.0*bet)*(pow(ama,2)+pow(rmz,2))/(pow(ama,2)-pow(rmz,2)));

      if(cos(2.0*bet)*(pow(ama,2)-pow(rmz,2)) > 0)
            alfarunp = alfarunp - pi/2;


      if(aml == 0.0)
      {
        aml=aml0;
        mlpole = aml0;
        }

        if(amh == 0.0)
        {
        amh=amh0;
        mhpole = amh0;
        }
      double amch;
      if(amch == 0.0)
      {
        amch=amch0;
        mchpole = amch0;
        }

        double amlight=aml;
      double amheavy=amh;

      aml = aml0;
      amh = amh0;
      amch = amch0;
      double ama0 = ama;

        //cccccccccccccccccccccc

        double pis1s1h,pis1s2h,pis2s2h,pis1s1l,pis1s2l,pis2s2l;
        double tad1,tad2;

        double picc,pizz,piww;
      su_hloop(q2,amlight,amu,au,ad,al,pis1s1l,pis1s2l,pis2s2l,piaa,picc,pizz,piww,tad1,tad2);
      su_hloop(q2,amheavy,amu,au,ad,al,pis1s1h,pis1s2h,pis2s2h,piaa,picc,pizz,piww,tad1,tad2);

      double vd2 = 2*(pow(amz,2)+pizz)/(pow(g1ewsb,2)+pow(g2ewsb,2))/(1.0+pow(tbeta,2));
      double vu2 = vd2*pow(tbeta,2);
      double vev2=2.0*(vu2+vd2);
        double rmldr = ytau*sqrt(vd2);
      double rmbdr = yb*sqrt(vd2);
      double rmtdr = yt*sqrt(vu2);
      double gstrong=sqrt(4.0*pi*alsewsb);
      double sxt=sin(thett);
      double cxt=cos(thett);
      double sxb=sin(thetb);
      double cxb=cos(thetb);
      double cxl=cos(thetl);
      double sxl=sin(thetl);

      pizzp = pizz;
      double ihdr=0.0;

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%   two--loop alphas corrections (p. slavich)
        double ctx,s11s,s22s,s12s,p2s;
        double s22b,s11b,s12b,p2b;
        double s11w,s12w,s22w,p2w;
        double tad1st,tad2st,tad2sb,tad1sb,tad1w,tad2w;
        double s11bl,s12bl,s22bl;
        double p2bl;

l1:   if(imodel >= 2)
        {
            // full one-loop higgs calculation but neglecting any two-loops
            printf("Sorry two loop Higgs corrections are not added in this code ");
            s11s = 0.0;
        s12s= 0.0;
        s22s=0.0;
        s11b = 0.0;
        s12b= 0.0;
        s22b=0.0;
        p2s=0.0;
        p2b=0.0;
        s11w = 0.0;
        s12w= 0.0;
        s22w=0.0;
        p2w =0.0;
        s11bl = 0.0;
        s12bl= 0.0;
        s22bl=0.0;
        p2bl=0.0;
            s11l = 0.0;
        s12l= 0.0;
        s22l=0.0;
        p2l=0.0;

        //su_dszhiggs(pow(rmtdr,2),am3,pow(mst1,2),pow(mst2,2),sxt,ctx,pow(scale,2), -amu,tbeta,vev2,gstrong,0,s11s,s22s,s12s);
        //su_dszhiggs(pow(rmbdr,2),am3,pow(msb1,2),pow(msb2,2),sxb,cxb,pow(scale,2), -amu,1.0/tbeta,vev2,gstrong,0,s22b,s11b,s12b);
        //su_dszodd(pow(rmtdr,2),am3,pow(mst1,2),pow(mst2,2),sxt,cxt,pow(scale,2),  -amu,tbeta,vev2,gstrong,p2s);
        //su_dszodd(pow(rmbdr,2),am3,pow(msb1,2),pow(msb2,2),sxb,cxb,pow(scale,2), -amu,1.0/tbeta,vev2,gstrong,p2b);

            //		two-loop electroweak corrections (p. slavich routines)
        //su_ddshiggs(pow(rmtdr,2),pow(rmbdr,2),pow(amar,2),pow(mst1,2),pow(mst2,2),pow(msb1,2),pow(msb2,2),sxt,cxt,sxb,cxpow(b,scale,2),-amu,tbeta,vev2, s11w,s12w,s22w);
        //su_ddsodd(pow(rmtdr,2),pow(rmbdr,2),pow(amar,2),pow(mst1,2),pow(mst2,2),pow(msb1,2),pow(msb2,2),sxt,cxt,sxb,cxpow(b,scale,2),-amu,tbeta,vev2,p2w);

            //		now add the tau-lepton contributions.

            //su_taubot(pow(rmldr,2),pow(rmbdr,2),pow(msta1,2),pow(msta2,2),pow(msb1,2),pow(msb2,2),sxl,cxl,sxb,cxpow(b,scale,2),-amu,tbeta,vev2, s11bl,s12bl,s22bl);
        //su_taubotodd(pow(rmldr,2),pow(rmbdr,2),pow(msta1,2),pow(msta2,2),pow(msb1,2),pow(msb2,2),sxl,cxl,sxb,cxpow(b,scale,2),-amu,tbeta,vev2, p2bl);
        //su_tausqhiggs(pow(rmldr,2),pow(amar,2),pow(msntau,2),pow(msta1,2),pow(msta2,2),sxl,cxl,pow(scale,2),-amu,tbeta,vev2,0,s11l,s22l,s12l);
        //su_tausqodd(pow(rmldr,2),pow(amar,2),pow(msntau,2),pow(msta1,2),pow(msta2,2), sxl,cxl,pow(scale,2),-amu,tbeta,vev2,p2l);

            //		2-loop tadpole corrections (p. slavich routines)
        //su_ewsb2loop(pow(rmtdr,2),apow(m3,mst1,2),pow(mst2,2),sxt,cxt, pow(scale,2),-amu,tbeta,vev2,gstrong,tad1st,tad2st);
        //su_ewsb2loop(pow(rmbdr,2),apow(m3,msb1,2),pow(msb2,2),sxb,cxb, pow(scale,2),-amu,1.0/tbeta,vev2,gstrong,tad2sb,tad1sb);
        //su_ddstad(pow(rmtdr,2),pow(rmbdr,2),pow(amar,2),pow(mst1,2),pow(mst2,2), pow(msb1,2),pow(msb2,2),sxt,cxt,sxb,cxpow(b,scale,2),-amu,tbeta,vev2, tad1w,tad2w);
        //su_taubottad(pow(rmldr,2),pow(rmbdr,2),pow(msta1,2),pow(msta2,2),pow(msb1,2), pow(msb2,2),sxl,cxl,sxb,cxpow(b,scale,2),-amu,tbeta,vev2,tad1bl,tad2bl);
        //su_tausqtad(pow(rmldr,2),pow(amar,2),pow(msntau,2),pow(msta1,2),pow(msta2,2),sxl,cxl,pow(scale,2),-amu,tbeta,vev2, tad1l,tad2l);
            }
      else
      {
            // full one-loop higgs calculation but neglecting any two-loops
        s11s = 0.0;
        s12s= 0.0;
        s22s=0.0;
        s11b = 0.0;
        s12b= 0.0;
        s22b=0.0;
        p2s=0.0;
        p2b=0.0;
        s11w = 0.0;
        s12w= 0.0;
        s22w=0.0;
        p2w =0.0;
        s11bl = 0.0;
        s12bl= 0.0;
        s22bl=0.0;
        p2bl=0.0;
            s11l = 0.0;
        s12l= 0.0;
        s22l=0.0;
        p2l=0.0;
        }

        //     add two-loop tadpoles in running ma:

      if(imodel >= 2)
        tad1loop= tad1st+tad1sb+tad1w+tad1l+tad1bl;

      double dvdvd2=-tad1;

      if(imodel >= 2)
        dvdvd2=dvdvd2+tad1loop;

      if(imodel >= 2)
        tad2loop=tad2st+tad2sb+tad2w+tad2l+tad2bl;

      double dvdvu2=-tad2;

      if(imodel >= 2)
        dvdvu2=dvdvu2+tad2loop;

      double mz2dr = pow(amz,2)+pizz;
      double madr2=(mhu2+dvdvu2 -mhd2-dvdvd2)/cos(2*bet)-mz2dr;

      dma = p2s+p2w+p2b+p2l+p2bl;

      if(kmaflag == 0)
      {   								 //!! then ama is really ma pole input
        ama=marsave;
            madr2 = pow(ama,2) +piaa -pow(sb,2)*tad1-pow(cb,2)*tad2 -dma;
        madr2save=madr2;
        }

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ml11loop = mz2dr*pow(cb,2)+madr2*pow(sb,2)-pis1s1l+tad1+ s11s+s11w+s11b+s11l+s11bl+dma*pow(sb,2);
      ml22loop = mz2dr*pow(sb,2)+madr2*pow(cb,2)-pis2s2l+tad2+ s22s+s22w+s22b+s22l+s22bl+dma*pow(cb,2);
      ml12loop = -(mz2dr+madr2)*sb*cb-pis1s2l+ s12s+s12w+s12b+s12l+s12bl-dma*sb*cb;

      mh11loop = mz2dr*pow(cb,2)+madr2*pow(sb,2)-pis1s1h+tad1+ s11s+s11w+s11b+s11l+s11bl+dma*pow(sb,2);
      mh22loop = mz2dr*pow(sb,2)+madr2*pow(cb,2)-pis2s2h+tad2+ s22s+s22w+s22b+s22l+s22bl+dma*pow(cb,2);
      mh12loop = -(mz2dr+madr2)*sb*cb-pis1s2h+s12s+s12w+s12b+s12l+s12bl-dma*sb*cb;

      double mlcr2=0.50*(ml11loop+ml22loop-sqrt(pow((ml11loop-ml22loop),2) +4*mpow(l12loop,2)));
      double mhcr2=0.50*(mh11loop+mh22loop+sqrt(pow((mh11loop-mh22loop),2) +4*mpow(h12loop,2)));

      if(mlcr2 >= 0.0)
         mlpole=sqrt(mlcr2);
      else
      {
         mlpole=aml0;
        if(irge == irgmax && ifix >= 3)
            mlpole=sqrt(mlcr2);
        }

      if(mhcr2 >= 0.0)
         mhpole=sqrt(mhcr2);
      else
      {
         mhpole=amh0;
        if(irge == irgmax && ifix >= 3)
            mhpole=sqrt(mhcr2);
        }

      double mh2dum=0.50*(ml11loop+ml22loop+sqrt(pow((ml11loop-ml22loop),2) +4*mpow(l12loop,2)));

      double s2alfa=2.0*ml12loop/(mh2dum-pow(mlpole,2));
      double c2alfa= (ml11loop-ml22loop)/(mh2dum-pow(mlpole,2));
      double t2alfa=s2alfa/c2alfa;

        // Next is to have correct alpha angle convention:

      if(c2alfa > 0.0)
            a=0.50*atan(t2alfa);


      if(c2alfa < 0.0)
      {
        if(s2alfa < 0.0)
            a=0.50*atan(t2alfa)-pi/2;
        else
            a=0.50*atan(t2alfa)+pi/2;
        }

      tadba=pow(sb,2)*tad1+pow(cb,2)*tad2;
      double macr2 =madr2-piaa+tadba+dma;
      double mchcr2=macr2+pow(amw,2)+piaa-picc+piww;

      if(macr2 >= 0.0)
          mapole = sqrt(macr2);
      else
      {
        mapole = ama;
          if(irge == irgmax && ifix >= 3)
            mapole=sqrt(macr2);
        }

      if(mchcr2 >= 0.0)
        mchpole = sqrt(mchcr2);
      else
      {
        mchpole = amch0;
        if(irge == irgmax && ifix >= 3)
            mchpole=sqrt(mchcr2);
        }
        //=========  end of the full calculation
        }
    */

    double la3t = la3 + la4 + la5;
    double ama2 = pow(amar, 2);
    double aml2 = pow(amlr, 2);
    double amh2 = pow(amhr, 2);
    double amp2 = pow(amchr, 2);

    //  Higgs Couplings
    sbma = sb * ca - cb * sa;
    cbma = cb * ca + sb * sa;
    sbpa = sb * ca + cb * sa;
    cbpa = cb * ca - sb * sa;

    double s2a = 2 * sa * ca;
    double c2a = pow(ca, 2) - pow(sa, 2);
    double s2b = 2 * sb * cb;
    double c2b = pow(cb, 2) - pow(sb, 2);

    double glzz = 1.0 / vlocal / 2 * aml2 * sbma;
    double ghzz = 1.0 / vlocal / 2 * amh2 * cbma;
    double glww = 2 * glzz;
    double ghww = 2 * ghzz;
    double glaz = 1.0 / vlocal * (aml2 - ama2) * cbma;
    double ghaz = -1.0 / vlocal * (amh2 - ama2) * sbma;
    double glpw = -1.0 / vlocal * (amp2 - aml2) * cbma;
    double glmw = glpw;
    double ghpw = 1.0 / vlocal * (amp2 - amh2) * sbma;
    double ghmw = ghpw;
    double gapw = 1.0 / vlocal * (amp2 - ama2);
    double gamw = -gapw;
    double ghhh =     vlocal / 2 * (la1 * pow(ca, 3) * cb + la2 * pow(sa, 3) * sb + la3t * sa * ca * sbpa + la6 * pow(ca, 2) * (3 * sa * cb + ca * sb) + la7 * pow(sa, 2) * (3 * ca * sb + sa * cb));
    double glll =   - vlocal / 2 * (la1 * pow(sa, 3) * cb - la2 * pow(ca, 3) * sb + la3t * sa * ca * cbpa - la6 * pow(sa, 2) * (3 * ca * cb - sa * sb) + la7 * pow(ca, 2) * (3 * sa * sb - ca * cb));
    double glhh = -3 *vlocal / 2 * (la1 * pow(ca, 2) * cb * sa - la2 * pow(sa, 2) * sb * ca + la3t * (pow(sa, 3) * cb - pow(ca, 3) * sb + 2 * sbma / 3) - la6 * ca * (cb * c2a - sa * sbpa) - la7 * sa * (c2a * sb + ca * sbpa));
    double ghll = 3 * vlocal / 2 * (la1 * pow(sa, 2) * cb * ca + la2 * pow(ca, 2) * sb * sa + la3t * (pow(sa, 3) * sb + pow(ca, 3) * cb - 2 * cbma / 3) - la6 * sa * (cb * c2a + ca * cbpa) + la7 * ca * (c2a * sb + sa * cbpa));
    double glaa =   - vlocal / 2 * (la1 * pow(sb, 2) * cb * sa - la2 * pow(cb, 2) * sb * ca - la3t * (pow(sb, 3) * ca - pow(cb, 3) * sa) + 2 * la5 * sbma - la6 * sb * (cb * sbpa + sa * c2b) - la7 * cb * (c2b * ca - sb * sbpa));
    double ghaa =     vlocal / 2 * (la1 * pow(sb, 2) * cb * ca + la2 * pow(cb, 2) * sb * sa + la3t * (pow(sb, 3) * sa + pow(cb, 3) * ca) - 2 * la5 * cbma - la6 * sb * (cb * cbpa + ca * c2b) + la7 * cb * (sb * cbpa + sa * c2b));
    double glpm = 2 * glaa + vlocal * (la5 - la4) * sbma;
    double ghpm = 2 * ghaa + vlocal * (la5 - la4) * cbma;

    glzz = 2 * glzz;
    ghzz = 2 * ghzz;
    glll = 6 * glll;
    ghhh = 6 * ghhh;
    glhh = 2 * glhh;
    ghll = 2 * ghll;
    glaa = 2 * glaa;
    ghaa = 2 * ghaa;
    double xnorm = pow(amz, 2) / vlocal;
    glll = glll / xnorm;
    ghll = ghll / xnorm;
    glhh = glhh / xnorm;
    ghhh = ghhh / xnorm;
    ghaa = ghaa / xnorm;
    glaa = glaa / xnorm;
    glpm = glpm / xnorm;
    ghpm = ghpm / xnorm;
    gat = 1.0 / tgbet;
    gab = tgbet;
    glt = ca / sb;
    glb = -sa / cb;
    ght = sa / sb;
    ghb = ca / cb;
    double gzal = -cbma;
    double gzah = sbma;
    glvv = sbma;
    ghvv = cbma;
    b = bet;

    //  higgs couplings needed in suspect:
    alfa = a;
    double xgat = gat;
    double xgab = gab;
    double xglt = glt;
    double xglb = glb;
    double xght = ght;
    double xghb = ghb;
    double xghvv = ghvv;
    double xglvv = glvv;

    // ===============================================================
    // ========== pole higgs masses (but at one-loop)
    //  if(imodel == 1 || imodel == 0)then
    //  affects only pure 1-loop higgs mass choice:

    double xdlb, xdht, xdhb, xdat, xdab, xdlsb, xdhst, xdhsb, xdlst;
    double xdlt;

    if (imodel == 0)
    {
        xdlt = gf / (2.0 * sqrt(2.0) * pow(pi, 2)) * pow(glt, 2) * (-2.0 * pow(mt, 2) + 0.50 * aml2) * creal(f0_hdec(mt, mt, aml2)) * 3 * pow(mt, 2);
        xdlb = gf / (2.0 * sqrt(2.0) * pow(pi, 2)) * pow(glb, 2) * (-2.0 * pow(mb, 2) + 0.50 * aml2) * creal(f0_hdec(mb, mb, aml2)) * 3 * pow(mb, 2) + gf / (2.0 * sqrt(2.0) * pow(pi, 2)) * pow(glb, 2) * (0.50 * aml2) * log(pow(mb, 2) / pow(mt, 2)) * 3 * pow(mb, 2);
        xdht = gf / (2.0 * sqrt(2.0) * pow(pi, 2)) * pow(ght, 2) * (-2.0 * pow(mt, 2) + 0.50 * amh2) * creal(f0_hdec(mt, mt, amh2)) * 3 * pow(mt, 2);
        xdhb = gf / (2.0 * sqrt(2.0) * pow(pi, 2)) * pow(ghb, 2) * (-2.0 * pow(mb, 2) + 0.50 * amh2) * creal(f0_hdec(mb, mb, amh2)) * 3 * pow(mb, 2) + gf / (2.0 * sqrt(2.0) * pow(pi, 2)) * pow(ghb, 2) * (0.50 * amh2) * log(pow(mb, 2) / pow(mt, 2)) * 3 * pow(mb, 2);
        xdat = gf / (2.0 * sqrt(2.0) * pow(pi, 2)) * pow(gat, 2) * (-0.50 * ama2) * creal(f0_hdec(mt, mt, ama2)) * 3 * pow(mt, 2);
        xdab = gf / (2.0 * sqrt(2.0) * pow(pi, 2)) * pow(gab, 2) * (-0.50 * ama2) * creal(f0_hdec(mb, mb, ama2)) * 3 * pow(mb, 2) + gf / (2.0 * sqrt(2.0) * pow(pi, 2)) * pow(gab, 2) * (-0.50 * ama2) * log(pow(mb, 2) / pow(mt, 2)) * 3 * pow(mb, 2);
        xdlst = 0.0;
        xdlsb = 0.0;
        xdhst = 0.0;
        xdhsb = 0.0;

        for (int i = 1; i <= 2; i++)
            for (int j = 1; j <= 2; j++)
            {
                xdlst = xdlst + gf / (2.0 * sqrt(2.0) * pow(pi, 2)) * pow(gltt[i][j], 2) * creal(f0_hdec(mst[i], mst[j], aml2)) * 3 * pow(amz, 4);
                xdlsb = xdlsb + gf / (2.0 * sqrt(2.0) * pow(pi, 2)) * pow(glbb[i][j], 2) * creal(f0_hdec(msb[i], msb[j], aml2)) * 3 * pow(amz, 4);
                xdhst = xdhst + gf / (2.0 * sqrt(2.0) * pow(pi, 2)) * pow(ghtt[i][j], 2) * creal(f0_hdec(mst[i], mst[j], amh2)) * 3 * pow(amz, 4);
                xdhsb = xdhsb + gf / (2.0 * sqrt(2.0) * pow(pi, 2)) * pow(ghbb[i][j], 2) * creal(f0_hdec(msb[i], msb[j], amh2)) * 3 * pow(amz, 4);
            }

        double xdast = gf / (1.0 * sqrt(2.0) * pow(pi, 2)) * pow(gatt, 2) * creal(f0_hdec(mst[1], mst[2], ama2)) * 3 * pow(amz, 4);
        double xdasb = gf / (1.0 * sqrt(2.0) * pow(pi, 2)) * pow(gabb, 2) * creal(f0_hdec(msb[1], msb[2], ama2)) * 3 * pow(amz, 4);

        aml = sqrt(aml2 + xdlt + xdlb + xdlst + xdlsb);
        amh = sqrt(amh2 + xdht + xdhb + xdhst + xdhsb);
        ama = sqrt(ama2 + xdat + xdab + xdast + xdasb);
        amch = sqrt(pow(ama, 2) + pow(amw, 2));
    }

    else
    {
        aml = mlpole;
        amh = mhpole;
        ama = mapole;
        amch = mchpole;
        alfa = a;
    }

    ma = ama;
    marun = amar;
    aama = ama;
    aaml = aml;
    aamh = amh;
    aamch = amch;

    return 0;
}
