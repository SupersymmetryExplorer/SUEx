#include <stdio.h>
#include <math.h>
#include "variables.h"
#include "complex.h"
#include "electroweak.h"
#include "numericx.h"
#include "functions.h"


static double zfunz;

// ============================================================================================================
// 	All the calculations of Lambda
// ============================================================================================================

// Function xitla()
// ------------------------------------------------------------------------------------------------------------

//	On the running coupling constant in QCD
// M. Prosperi, M. Raciti and C. Simolo
// arXiv:hep-ph/0607209v2
// ------------------------------------------------------------------------------------------------------------

// nf : Number of active quarks, which changes by +/- 1 anytime \mu crosses a quark thresold mf   ... Page No. -4

#define b0(nf) (33.0 - 2.0 * nf)                                                                                                                                   // ... Eq(16)
#define b1(nf) (6.0 * (153.0 - 19.0 * nf) / pow(b0(nf), 2))                                                                                                        // ... Eq(16)
#define als2(nf, x, xlb) (12.0 * pi / (b0(nf) * log(pow(x, 2) / pow(xlb, 2))) * (1.0 - b1(nf) * log(log(pow(x, 2) / pow(xlb, 2))) / log(pow(x, 2) / pow(xlb, 2)))) // ... Eq(29)
#define aa(nf) (12.0 * pi / b0(nf))
#define bb(nf) (b1(nf) / aa(nf))
#define xit(a, b, x) (a / 2.0 * (1.0 + sqrt(1.0 - 4 * b * log(x)))) //... eq(29)  little modified

//	On the running coupling constant in QCD
// M. Prosperi, M. Raciti and C. Simolo
// arXiv:hep-ph/0607209v2
// ------------------------------------------------------------------------------------------------------------

double xitla(int no, double alp, double acc) // acc : is the tolerance value
{                                            // no : Number of loops.
                                             // Iteration function to determine improved lambda's
    // This function can determine 1-loop or 2-loop lambda's

    double q, xlb, xx, y1, y2, dy, x;
    double a, b;

    double nf = 5;
    q = mz;
    xlb = q * exp(-aa(nf) / alp / 2.0); // Scale parameter Lambda. Page 4, eq (4)
                                        // For one loop
    if (no == 1)
        return (xlb);

    do // For two loop. Page 13, eq(23), eq(25)
    {
        x = log(pow(q, 2) / pow(xlb, 2)); // We choose mu to be mz. Determine Lambda such that alpha remains constant.
        a = aa(nf) / alp;                 //
        b = bb(nf) * alp;
        xx = xit(a, b, x);
        xlb = q * exp(-xx / 2.0);
        y1 = alp;
        y2 = als2(nf, q, xlb);
        dy = fabs(y2 - y1) / y1;
    } while (dy >= acc);

    return (xlb);
}

//	Function : xiter()
// ------------------------------------------------------------------------------------------------------------

//	On the running coupling constant in QCD
// M. Prosperi, M. Raciti and C. Simolo
// arXiv:hep-ph/0607209v2
// ------------------------------------------------------------------------------------------------------------

double xiter(double q, double xlb1, int nf1, double xlb, int nf2, double acc)
{
    double y1, y2, dy, xlb2, xx, x, alp;
    double a, b;

    xlb2 = xlb;

    do
    {
        x = log(pow(q, 2) / pow(xlb2, 2));
        alp = als2(nf1, q, xlb1);
        a = aa(nf2) / alp;
        b = bb(nf2) * alp;
        xx = xit(a, b, x);
        xlb2 = q * exp(-xx / 2);
        y1 = als2(nf1, q, xlb1);
        y2 = als2(nf2, q, xlb2);
        dy = fabs(y2 - y1) / y1;
    } while (dy >= acc);

    return (xlb2);
}

#undef b0
#undef b1
#undef als2
#undef aa
#undef bb
#undef xit

//	Function : alsini()
// ------------------------------------------------------------------------------------------------------------

//	On the running coupling constant in QCD
// M. Prosperi, M. Raciti and C. Simolo
// arXiv:hep-ph/0607209v2
// ------------------------------------------------------------------------------------------------------------

void alsini(double acc)
{
    double xlb[7];
    double amc = mc0;
    double amb = mb0;
    double amt = mt0;

    xlb1[1] = 0.0;
    xlb1[2] = 0.0;
    xlb2[1] = 0.0;
    xlb2[2] = 0.0;

    // Equation (35), Page No. 17
    // For 1- loop calculations.

    if (n0 == 3)
    {
        xlb[3] = xlambda;
        xlb[4] = xlb[3] * pow((xlb[3] / amc), (2.0 / 25.0));
        xlb[5] = xlb[4] * pow((xlb[4] / amb), (2.0 / 23.0));
        xlb[6] = xlb[5] * pow((xlb[5] / amt), (2.0 / 21.0));
    }
    else if (n0 == 4)
    {
        xlb[4] = xlambda;
        xlb[5] = xlb[4] * pow((xlb[4] / amb), (2.0 / 23.0));
        xlb[3] = xlb[4] * pow((xlb[4] / amc), (-2.0 / 27.0));
        xlb[6] = xlb[5] * pow((xlb[5] / amt), (2.0 / 21.0));
    }
    else if (n0 == 5)
    {
        xlb[5] = xlambda;
        xlb[4] = xlb[5] * pow((xlb[5] / amb), (-2.0 / 25.0));
        xlb[3] = xlb[4] * pow((xlb[4] / amc), (-2.0 / 27.0));
        xlb[6] = xlb[5] * pow((xlb[5] / amt), (2.0 / 21.0));
    }
    else if (n0 == 6)
    {
        xlb[6] = xlambda;
        xlb[5] = xlb[6] * pow((xlb[6] / amt), (-2.0 / 23.0));
        xlb[4] = xlb[5] * pow((xlb[5] / amb), (-2.0 / 25.0));
        xlb[3] = xlb[4] * pow((xlb[4] / amc), (-2.0 / 27.0));
    }

    for (int i = 1; i <= 6; i++)
        xlb1[i] = xlb[i];

    // Same thing as before for 2-loop calculations.
    // For Details go through Page 17

    if (n0 == 3)
    {
        xlb[3] = xlambda;
        xlb[4] = xlb[3] * pow((xlb[3] / amc), (2.0 / 25.0)) * pow((2.0 * log(amc / xlb[3])), (-107.0 / 1875.0));
        xlb[4] = xiter(amc, xlb[3], 3, xlb[4], 4, acc);
        xlb[5] = xlb[4] * pow((xlb[4] / amb), (2.0 / 23.0)) * pow((2.0 * log(amb / xlb[4])), (-963.0 / 13225.0));
        xlb[5] = xiter(amb, xlb[4], 4, xlb[5], 5, acc);
        xlb[6] = xlb[5] * pow((xlb[5] / amt), (2.0 / 21.0)) * pow((2.0 * log(amt / xlb[5])), (-321.0 / 3381.0));
        xlb[6] = xiter(amt, xlb[5], 5, xlb[6], 6, acc);
    }
    else if (n0 == 4)
    {
        xlb[4] = xlambda;
        xlb[5] = xlb[4] * pow((xlb[4] / amb), (2.0 / 23.0)) * pow((2.0 * log(amb / xlb[4])), (-963.0 / 13225.0));
        xlb[5] = xiter(amb, xlb[4], 4, xlb[5], 5, acc);
        xlb[3] = xlb[4] * pow((xlb[4] / amc), (-2.0 / 27.0)) * pow((2.0 * log(amc / xlb[4])), (107.0 / 2025.0));
        xlb[3] = xiter(amc, xlb[4], 4, xlb[3], 3, acc);
        xlb[6] = xlb[5] * pow((xlb[5] / amt), (2.0 / 21.0)) * pow((2.0 * log(amt / xlb[5])), (-321.0 / 3381.0));
        xlb[6] = xiter(amt, xlb[5], 5, xlb[6], 6, acc);
    }
    else if (n0 == 5)
    {
        xlb[5] = xlambda;
        xlb[4] = xlb[5] * pow((xlb[5] / amb), (-2.0 / 25.0)) * pow((2.0 * log(amb / xlb[5])), (963.0 / 14375.0));
        xlb[4] = xiter(amb, xlb[5], 5, xlb[4], 4, acc);
        xlb[3] = xlb[4] * pow((xlb[4] / amc), (-2.0 / 27.0)) * pow((2.0 * log(amc / xlb[4])), (107.0 / 2025.0));
        xlb[3] = xiter(amc, xlb[4], 4, xlb[3], 3, acc);
        xlb[6] = xlb[5] * pow((xlb[5] / amt), (2.0 / 21.0)) * pow((2.0 * log(amt / xlb[5])), (-321.0 / 3381.0));
        xlb[6] = xiter(amt, xlb[5], 5, xlb[6], 6, acc);
    }
    else if (n0 == 6)
    {
        xlb[6] = xlambda;
        xlb[5] = xlb[6] * pow((xlb[6] / amt), (-2.0 / 23.0)) * pow((2.0 * log(amt / xlb[6])), (321.0 / 3703.0));
        xlb[5] = xiter(amt, xlb[6], 6, xlb[5], 5, acc);
        xlb[4] = xlb[5] * pow((xlb[5] / amb), (-2.0 / 25.0)) * pow((2.0 * log(amb / xlb[5])), (963.0 / 14375.0));
        xlb[4] = xiter(amb, xlb[5], 5, xlb[4], 4, acc);
        xlb[3] = xlb[4] * pow((xlb[4] / amc), (-2.0 / 27.0)) * pow((2.0 * log(amc / xlb[4])), (107.0 / 2025.0));
        xlb[3] = xiter(amc, xlb[4], 4, xlb[3], 3, acc);
    }

    for (int i = 1; i <= 6; i++)
        xlb2[i] = xlb[i];
}

// ============================================================================================================
// -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . --
// ============================================================================================================

// ============================================================================================================
// 	All the calculations of Running QCD Coupling constants : alpha_s
// ============================================================================================================

//  Running of the QCD coupling at scale Q and perturbative order N:
//  alphas()

//	On the running coupling constant in QCD
// M. Prosperi, M. Raciti and C. Simolo
// arXiv:hep-ph/0607209v2
// ------------------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------------------

#define b0(nf) (33.0 - 2.0 * nf)
#define b1(nf) (6.0 * (153.0 - 19.0 * nf) / pow(b0(nf), 2))
#define als1(nf, x) (12.0 * pi / (b0(nf) * log(pow(x, 2) / pow(xlb[nf], 2))))
#define als2(nf, x) (12.0 * pi / (b0(nf) * log(pow(x, 2) / pow(xlb[nf], 2))) * (1.0 - b1(nf) * log(log(pow(x, 2) / pow(xlb[nf], 2))) / log(pow(x, 2) / pow(xlb[nf], 2))))

double falphas(double q, int n)
{
    double xlb[7];
    int nf;

    if (n == 1)
    {
        for (int i = 1; i <= 6; i++)
            xlb[i] = xlb1[i];
    }
    else
    {
        for (int i = 1; i <= 6; i++)
            xlb[i] = xlb2[i];
    }

    if (q < mc0)
        nf = 3;
    else if (q <= mb0)
        nf = 4;
    else if (q <= mt0)
        nf = 5;
    else
        nf = 6;

    if (n == 1)
        return (als1(nf, q)); // Page 4, eq(6)
    else
        return (als2(nf, q)); // page 14, eq(29)
}

#undef b0
#undef b1
#undef als1
#undef als2

// ============================================================================================================
// -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . --
// ============================================================================================================

// ============================================================================================================
// 	Calculating the pole mass of the quarks. Running masses of the quarks
// ============================================================================================================

// Function : runm()
// The running of the quark masses at scale Q and with NF quark flavors:
// ------------------------------------------------------------------------------------------------------------

// QCD Corrections to Hardonic Higgs Decays
// A. Djouadi, M. Spira, P.M.Zerwas
//
// Scheme dependence of the next-to-next-to-leading QCD corrections to Lambda_tot (H0 -> hardons) and the spurious QCD infrared fixed point
// S.G.Gorishny, A.L.Kataev and S.A.Larin
// Physical Review D, Vol 43, No 5

// A fresh look into the heavy quark-mass values
// S. Narison
// Physics LeRe~B 341 (1994) 73-83
// ------------------------------------------------------------------------------------------------------------

#define b0(nf) ((33.0 - 2.0 * nf) / 12.0)                                                                 // Ref 2.. Eq(2.7)
#define b1(nf) ((102.0 - 38.0 / 3.0 * nf) / 16.0)                                                         // Ref 2.. Eq(2.7)
#define b2(nf) ((2857.0 / 2.0 - 5033.0 / 18.0 * nf + 325.0 / 54.0 * pow(nf, 2)) / 64.0)                   // Ref 2.. Eq(2.7)
#define g0(nf) (1.0)                                                                                      // Ref 2.. Eq(2.8)
#define g1(nf) ((202.0 / 3.0 - 20.0 / 9.0 * nf) / 16.0)                                                   // Ref 2.. Eq(2.8)
#define g2(nf) ((1249.0 - (2216.0 / 27.0 + 160.0 / 3.0 * zeta3) * nf - 140.0 / 81.0 * pow(nf, 2)) / 64.0) // Ref 2.. Eq(2.8)
#define c1(nf) (g1(nf) / b0(nf) - b1(nf) * g0(nf) / pow(b0(nf), 2))                                       // Ref 2.. Eq(3.6)
#define c2(nf) ((pow((g1(nf) / b0(nf) - b1(nf) * g0(nf) / pow(b0(nf), 2)), 2) + g2(nf) / b0(nf) + pow(b1(nf), 2) * g0(nf) / pow(b0(nf), 3) - b1(nf) * g1(nf) / pow(b0(nf), 2) - b2(nf) * g0(nf) / pow(b0(nf), 2)) / 2.0)
#define tran(x, xk) (1.0 + coeff1 * falphas(x, 2) / pi + xk * pow((falphas(x, 2) / pi), 2) + coeff3 * pow((falphas(x, 2) / pi), 3))
#define cq(x, nf) (pow((2.0 * b0(nf) * x), (g0(nf) / b0(nf))) * (1.0 + c1(nf) * x + c2(nf) * pow(x, 2))) // Ref 2.. Eq(3.6)

double runm(double q, int nf)
{
    static int nn = 6;
    static double zeta3 = 1.202056903159594;
    int istrange = 0;

    // save istrange
    // 3-loop coeff3 in m_pole/m_running relation added 10/12/03

    double am[nn + 1], ymsb[nn + 1];
    double amt = mtpole;
    double coeff1, coeff3;

    double xmhat; // Check what this variable is doing here. I guess its not needed...
    double xmsb;
    double xkfac;
    double xk;
    double q0;

    nnlo = 1;
    // (always use nnlo now)
    // define the light quark masses
    double amc = 1.42;
    double amsb = 0.19;

    double ams = amsb;

    double acc = 1.0e-8; // Accuracy parameter

    if (idrflag != 1)
        coeff1 = 4.0 / 3.0;
    else
        coeff1 = 5.0 / 3.0;

    if (nnlo == 0)
        coeff3 = 0.0;
    else
        coeff3 = 101.45424;

    am[1] = 0.0;
    am[2] = 0.0;

    double amsd, amsu;
    double dd;

    int imsbar = 0; // Some mistake is there.. CHeck..

    if (imsbar == 1)
    {
        if (istrange == 0)
        {
            // --strange pole mass from msbar-mass at 1 GeV
            amsd = xlambda;
            amsu = 1.e8;
        L123:
            ams = (amsu + amsd) / 2;
            am[3] = ams;
            xmsb = ams / cq(falphas(ams, 2) / pi, 3) * cq(falphas(1.0, 2) / pi, 3) / tran(ams, 0.0);
            dd = (xmsb - amsb) / amsb;
            if (fabs(dd) >= acc)
            {
                if (dd <= 0.0)
                    amsd = am[3];
                else
                    amsu = am[3];
                goto L123;
            }
            istrange = 1;
        }
        am[3] = amsb;
    }
    else
    {
        ams = amsb;
        am[3] = ams;
    }

    //------------------------------------------------------------------------------------------------------------------------------------
    // Modifs jlk: to determine (perturbatively, at an order consistent with the pert level used in runm) mb(pole) from mb(mb)_msbar input:
    // mbmb= mb(mb)_msbar ; mbpole determined iteratively to acc. d-8
    //------------------------------------------------------------------------------------------------------------------------------------
    double xkb;
    double mbmbpole;
    double mbsave;
    if (imbmb == 0)
    {
        // imbmb is just a flag because this calculation is only needed once
        for (int i = 1; i <= 20; i++)
        {
            if (i == 1)
            {
                mbsave = 0.0;
                mbpole = mbmb;
            }
            if (nnlo == 0)
                xkb = 0.0;
            else
                xkb = 16.11 - 1.04 * (4.0 - (amsb + amc) / mbpole);

            if (i >= 3)
            {
                mb0 = mbpole; // amba=mbpole;
                alsini(1.e-8);
            }

            mbmbpole = mbmb * cq(falphas(mbpole, 2) / pi, 4) / cq(falphas(mbmb, 2) / pi, 4); // mbmbpole is mb(mbpole)
            mbpole = mbmbpole * tran(mbpole, xkb);                                           // tran(q,xk) is the usual pert.
                                                                                             // relation between mpole and mrun(mpole), see its def. above

            if (fabs(1.0 - mbsave / mbpole) < 1.e-8)
                break;
            mbsave = mbpole;
        }
        imbmb = 1;
    }

    // rest of calculation follows as before:
    //------------------------------------------------------------------------------------------------------

    am[3] = amsb;   // Strange quark mass
    am[4] = amc;    // Charm quark mass
    am[5] = mbpole; // amb									// Bottom quark mass
    am[6] = amt;    // Top quark mass
    xk = 16.110;

    for (int i = 1; i <= nf - 1; i++)
        xk = xk - 1.040 * (1.0 - am[i] / am[nf]);

    if (nf >= 4)
    {
        xmsb = am[nf] / tran(am[nf], 0.0);
        xmhat = xmsb / cq(falphas(am[nf], 2) / pi, nf);
    }
    else
    {
        xmsb = 0.0;
        xmhat = 0.0;
    }

    // Mass evolution of the quarks
    // Ref 1. Page - 3 Eq(4)
    // -------------------------------------------------------------------------------------------------------

    ymsb[3] = amsb;

    if (nf == 3)
    {
        ymsb[4] = ymsb[3] * cq(falphas(am[4], 2) / pi, 3) / cq(falphas(1.0, 2) / pi, 3);
        ymsb[5] = ymsb[4] * cq(falphas(am[5], 2) / pi, 4) / cq(falphas(am[4], 2) / pi, 4);
        ymsb[6] = ymsb[5] * cq(falphas(am[6], 2) / pi, 5) / cq(falphas(am[5], 2) / pi, 5);
    }
    else if (nf == 4)
    {
        ymsb[4] = xmsb;
        ymsb[5] = ymsb[4] * cq(falphas(am[5], 2) / pi, 4) / cq(falphas(am[4], 2) / pi, 4);
        ymsb[6] = ymsb[5] * cq(falphas(am[6], 2) / pi, 5) / cq(falphas(am[5], 2) / pi, 5);
    }
    else if (nf == 5)
    {
        ymsb[5] = xmsb;
        ymsb[4] = ymsb[5] * cq(falphas(am[4], 2) / pi, 4) / cq(falphas(am[5], 2) / pi, 4);
        ymsb[6] = ymsb[5] * cq(falphas(am[6], 2) / pi, 5) / cq(falphas(am[5], 2) / pi, 5);
    }
    else if (nf == 6)
    {
        ymsb[6] = xmsb;
        ymsb[5] = ymsb[6] * cq(falphas(am[5], 2) / pi, 5) / cq(falphas(am[6], 2) / pi, 5);
        ymsb[4] = ymsb[5] * cq(falphas(am[4], 2) / pi, 4) / cq(falphas(am[5], 2) / pi, 4);
    }

    // Number of quaeks which will contribute in the calculations : n0
    // The normalization scale : q0
    // -------------------------------------------------------------------------------------------------------

    if (q < amc)
    {
        n0 = 3;
        q0 = 1.0;
    }
    else if (q <= mbpole)
    {
        n0 = 4;
        q0 = amc;
    }
    else if (q <= amt)
    {
        n0 = 5;
        q0 = mbpole; // amb;
    }
    else
    {
        n0 = 6;
        q0 = amt;
    }

    if (nnlo == 1 && nf > 3)
        xkfac = tran(am[nf], 0.0) / tran(am[nf], xk);
    else
        xkfac = 1.0;

    return (ymsb[n0] * cq(falphas(q, 2) / pi, n0) / cq(falphas(q0, 2) / pi, n0) * xkfac);
}

#undef b0
#undef b1
#undef b2
#undef g0
#undef g1
#undef g2
#undef c1
#undef c2
#undef tran
#undef cq

//
// The following routines are for the evaluation of the chargino/neutralino and the squark, slepton masses including the radiative corrections.
// a la pierce, bagger, matchev, zhang, hep-ph/9606211.
//	For these corrections, one needs the one- and two-loop passarino-veltman functions discussed above.
// -----------------------------------------------------------------------------------------------------------------------------------------

//
// Calculates the chargino and neutralino masses and mixing angles (with analytical expressions) including radiative corrections in the higgsino
// and gaugino limits). the input parameters at ewsb scale are:
// 	mu,m1.m2,m3: higgs mass parameter and gaugino mass parameters,
// 	b,a: atan(tan(beta)) and the mixing angle alpha in the higgs sector.
//
//	The output parameters are:
// 	mc: the two chargino masses,
// 	mn: the four neutralino masses (absolute values),
// 	mx: the four neutralino masses (including signs).
// The  mass values are ordered with increasing value. the diagonalizing (ordered) mass matrices u,v for charginos and z for neutralinos are
// given in the common block su_matino/u,v,z
// -----------------------------------------------------------------------------------------------------------------------------------------

void su_gaugino(double mu, double m1, double m2, double m3, double b, double a, double mc[], double mn[], double xmn[])
{

    int iord[5], irem[3];
    double xmc[3];
    double x[3][3];
    double bu[3], bv[3];
    double yz[5][5], zx[5][5];
    double ymn[5];
    double mzsave = mz;
    double mwsave = mw;
    double cw = 1.0 / sqrt(1.0 + pow((g1ewsb / g2ewsb), 2));
    double sw = g1ewsb / g2ewsb * cw;
    sw2 = pow(sw, 2);

    double pizz = pizzp;

    if (isnan(pizz) || mz * mz + pizz <= 0.0)
    {
        // protections added non-pert or nan pb, uses tree-level values temporarily:
        pizz = 0.0;
        if (irge == irgmax)
            inonpert = -1;
    }

    double rmz = sqrt(mz * mz + pizz);
    double rmw = rmz * cw;
    mz = rmz;
    mw = rmw;

    double sb = sin(b);
    double cb = cos(b);

    double m1save, m2save, musave;
    double rcm1, rcm2, rcmu;

    double mq3 = msq;
    double mu3 = mtr;
    double md3 = mbr;

    if (inorc == 1) // only at very end of calculation
    {
        // Adding r.c to m1, m2, mu:
        m1save = m1;
        m2save = m2;
        musave = mu;

        su_radcino(mel, muq, mq3, mu3, md3, ma, ytewsb, ybewsb, m1, m2, -mu, tan(b), &rcm1, &rcm2, &rcmu);

        m1 = m1 + rcm1;
        m2 = m2 + rcm2; // r.c to m1,m2,mu at very last iter
        mu = mu - rcmu;
    }

    //
    // Neutralino masses and matrix elements
    // Adding protection for the special(unrealistic) problem m1=m2:

    double m1eqm2 = 0.0;
    double m1sav2;

    if (fabs(m1 - m2) < 1.0e-4)
    {
        m1eqm2 = 1.0;
        m1sav2 = m1;
        m1 = m1 + 1.0e-3;
    }

    //
    // Analytic formulas for the neutralino masses and the neutralino mixing matrix
    // Ahmed A Maarouf, Amr Aboshousha
    // Physical Review D - July 1992
    // (Equations for Neutralinos are taken from this. Check Page 3, 4, 5 )
    // -------------------------------------------------------------------------------------------------------------------------------

    //
    // SUPERSYMMETRY,PARTI (THEORY)
    //  Revised October 2007 by Howard E. Haber (UC Santa Cruz)
    // -------------------------------------------------------------------------------------------------------------------------------

    double eps = -1.0e-10;
    double xc2 = (m1 * m2 - pow(mz, 2) - pow(mu, 2)) - 3.0 / 8.0 * pow((m1 + m2), 2);
    double xc3 = -1.0 / 8.0 * pow((m1 + m2), 3) + 1.0 / 2.0 * (m1 + m2) * (m1 * m2 - pow(mz, 2) - pow(mu, 2)) + (m1 + m2) * pow(mu, 2);
    xc3 += +(m1 * pow(cw, 2) + m2 * pow(sw, 2)) * pow(mz, 2) - mu * pow(mz, 2) * sin(2.0 * b);
    double xc4 = +(m1 * pow(cw, 2) + m2 * pow(sw, 2)) * mu * pow(mz, 2) * sin(2.0 * b) - m1 * m2 * pow(mu, 2);
    xc4 += +1.0 / 4.0 * (m1 + m2) * ((m1 + m2) * pow(mu, 2) + (m1 * pow(cw, 2) + m2 * pow(sw, 2)) * pow(mz, 2) - mu * pow(mz, 2) * sin(2.0 * b));
    xc4 += +1.0 / 16.0 * pow((m1 + m2), 2) * (m1 * m2 - pow(mz, 2) - pow(mu, 2)) - 3.0 / 256.0 * pow((m1 + m2), 4);
    double xs = -pow(xc3, 2) - 2.0 / 27.0 * pow(xc2, 3) + 8.0 / 3.0 * xc2 * xc4;
    double xu = -1.0 / 3.0 * pow(xc2, 2) - 4.0 * xc4;

    complex cxd, cxc, cxa, cxb, cx1, cx2, cx3;
    cxd = mult((-4 * pow(xu, 3) - 27 * pow(xs, 2)), initiallize(1.0, eps));
    cxc = mult(1.0 / 2.0, adddouble((cmult(initiallize(0.0, 1.0), csqrt(mult(1.0 / 27.0, cxd)))), -xs));

    cxa = mult(creal(cpow(cxc, (1.0 / 3.0))), initiallize(1.0, -eps));
    cxb = add(mult(8.0, cxa), mult(-8.0 / 3.0 * xc2, initiallize(1.0, -eps)));

    //   Masses and couplings:
    double x0 = (m1 + m2) / 4.0;
    cx1 = add(mult(1.0 / 2.0, cxa), mult(-xc2 / 6.0, initiallize(1.0, -eps)));
    cx2 = add(mult(-1.0 / 2.0, cxa), mult(-xc2 / 3.0, initiallize(1.0, -eps)));
    cx3 = cmult(initiallize(1.0, -eps), cdiv(xc3, csqrt(cxb)));

    xmn[1] = x0 - cabs(csqrt(cx1)) + cabs(csqrt(add(cx2, cx3)));
    xmn[2] = x0 + cabs(csqrt(cx1)) - cabs(csqrt(add(cx2, mult(-1.0, cx3))));
    xmn[3] = x0 - cabs(csqrt(cx1)) - cabs(csqrt(add(cx2, cx3)));
    xmn[4] = x0 + cabs(csqrt(cx1)) + cabs(csqrt(add(cx2, mult(-1.0, cx3))));

    for (int i = 1; i <= 4; i++)
    {
        mn[i] = fabs(xmn[i]);
        ymn[i] = xmn[i];
        zx[i][2] = -cw / sw * (m1 - xmn[i]) / (m2 - xmn[i]);
        zx[i][3] = (mu * (m2 - xmn[i]) * (m1 - xmn[i]) - pow(mz, 2) * sb * cb * ((m1 - m2) * pow(cw, 2) + m2 - xmn[i])) / mz / (m2 - xmn[i]) / sw / (mu * cb + xmn[i] * sb);
        zx[i][4] = (-xmn[i] * (m2 - xmn[i]) * (m1 - xmn[i]) - pow(mz, 2) * cb * cb * ((m1 - m2) * pow(cw, 2) + m2 - xmn[i])) / mz / (m2 - xmn[i]) / sw / (mu * cb + xmn[i] * sb);
        zx[i][1] = 1.0 / sqrt(1.0 + pow(zx[i][2], 2) + pow(zx[i][3], 2) + pow(zx[i][4], 2));
        yz[i][1] = zx[i][1];
        yz[i][2] = zx[i][2] * zx[i][1];
        yz[i][3] = zx[i][3] * zx[i][1];
        yz[i][4] = zx[i][4] * zx[i][1];
    }

    //     ============ ordering the disorder

    if (mn[3] == mn[4])
        mn[4] = mn[4] * (1.0 + 1.e-8); // !protection

    // (such a degeneracy (within d.p. accuracy) may happen for very large mu)
    double xx0 = min1(mn[1], mn[2], mn[3], mn[4]);
    double xx1 = max1(mn[1], mn[2], mn[3], mn[4]);
    int idummy = 1;

    for (int i = 1; i <= 4; i++)
    {
        if (mn[i] == xx0)
            iord[1] = i;
        else if (mn[i] == xx1)
            iord[4] = i;
        else
        {
            irem[idummy] = i;
            idummy = idummy + 1;
        }
    }

    if (mn[irem[1]] <= mn[irem[2]])
    {
        iord[2] = irem[1];
        iord[3] = irem[2];
    }
    else
    {
        iord[2] = irem[2];
        iord[3] = irem[1];
    }

    int i;

    for (int j = 1; j <= 4; j++)
    {
        i = iord[j];
        xmn[j] = ymn[i];

        dxmn[j] = xmn[j];
        mn[j] = fabs(ymn[i]);

        for (int i1 = 1; i1 <= 4; i1++)
            z[j][i1] = yz[i][i1];
    }

    //
    // Chargino masses and matrix elements
    // -----------------------------------------------------------------------------------------------------------------------

    //
    // The Search for Supersymmetry: Probing Physics Beyond the Standard Model
    // Howard E. Haber, Gordon L. Kane
    // Page No: 137, Appendix C
    // ------------------------------------------------------------------------------------------------------------------------

    double delta = fabs(b - .25 * pi);
    double ddd = mu * cos(b) + m2 * sin(b);
    double ccc = mu * sin(b) + m2 * cos(b);
    double phim, phip;

    if (delta < 0.01)
    {
        phim = pi / 4.0 - .50 * atan((m2 - mu) / (2.0 * mw));
        phip = phim;
    }
    else if (fabs(ccc) < 1.e-5)
    {
        phim = 0.0;
        phip = atan(sqrt(2.0) * mw * sin(b) / (m2 + 1.e-5));
    }
    else if (fabs(ddd) < 1.e-5)
    {
        phip = 0.0;
        phim = atan(sqrt(2.0) * mw * cos(b) / (m2 + 1.e-5));
    }
    else
    {
        double rad = sqrt(pow((pow(m2, 2) - pow(mu, 2)), 2) + 4.0 * pow(mw, 4) * pow(cos(2.0 * b), 2) + 4.0 * pow(mw, 2) * (pow(m2, 2) + pow(mu, 2) + 2.0 * m2 * mu * sin(2.0 * b)));
        phip = atan((rad - (pow(m2, 2) - pow(mu, 2) + 2.0 * pow(mw, 2) * cos(2.0 * b))) / (2.0 * sqrt(2.0) * mw * (mu * cos(b) + m2 * sin(b))));
        phim = atan((rad - (pow(m2, 2) - pow(mu, 2) - 2.0 * pow(mw, 2) * cos(2.0 * b))) / (2.0 * sqrt(2.0) * mw * (mu * sin(b) + m2 * cos(b))));
    }

    double cp = cos(phip);
    double sp = sin(phip);
    double cm = cos(phim);
    double sm = sin(phim);

    // My convention

    u[2][2] = cm;
    u[2][1] = -sm;
    u[1][2] = sm;
    u[1][1] = cm;
    v[1][1] = cp;
    v[1][2] = sp;
    v[2][1] = -sp;
    v[2][2] = cp;

    x[1][1] = m2;
    x[1][2] = sqrt(2.0) * mw * sin(b);
    x[2][1] = sqrt(2.0) * mw * cos(b);
    x[2][2] = mu;

    while (1)
    {
        xmc[1] = (u[1][1] * x[1][1] + u[1][2] * x[2][1]) * v[1][1] + (u[1][1] * x[1][2] + u[1][2] * x[2][2]) * v[1][2];
        xmc[2] = (u[2][1] * x[1][1] + u[2][2] * x[2][1]) * v[2][1] + (u[2][1] * x[1][2] + u[2][2] * x[2][2]) * v[2][2];

        if (xmc[1] < 0.0)
        {
            // some corrections to deal with case where both m1,m2 <0:
            v[1][1] = -v[1][1];
            v[1][2] = -v[1][2];
            continue;
        }

        if (xmc[2] < 0.0)
        {
            v[2][1] = -v[2][1];
            v[2][2] = -v[2][2];
            continue;
        }
        break;
    }

    double mtemp;
    if (xmc[1] > xmc[2])
    {
        mtemp = xmc[1];
        xmc[1] = xmc[2];
        xmc[2] = mtemp;

        for (int j = 1; j <= 2; j++)
        {
            bu[j] = u[1][j];
            u[1][j] = u[2][j];
            u[2][j] = bu[j];
            bv[j] = v[1][j];
            v[1][j] = v[2][j];
            v[2][j] = bv[j];
        }
    }

    mc[1] = fabs(xmc[1]);
    mc[2] = fabs(xmc[2]);

    // some saving
    if (m1eqm2 == 1.0)
        m1 = m1sav2; // added for m1=m2 pbs (see above)

    if (inorc == 1)
    {
        m1 = m1save;
        m2 = m2save;
        mu = musave;
    }

    mz = mzsave;
    mw = mwsave;
}

//
// Calculates the sfermion masses including corrections, and the mixing angles for the 3d generation sfermions.
//
// The input parameters at EWSB scale are:
// 	mql,mur,mdr,mel,mer,mql1,mur1,mdr1,mel1,mer1: sfermion mass terms,
// 	al,at,ab,mu: 3d generation trilinear couplings and the parameter mu
//
// The outputs are the sfermions masses: mst,msb,msl,msu,msd,mse,msn.
//
//	The masses are ordered such that the lightest is 1 and the heaviest is 2. The mixing angles of 3 generation sfermion are in the common block
// common/su_outmix/thet,theb,thel (to be treated with care because of the ordering of the sfermion masses, when compared to other calculations).
// nb this routine also calculates sfermion masses and mixing in a different (bpmz) convention used in several other subroutines
// (the latter are passed via common/su_bpew/..)
// ----------------------------------------------------------------------------------------------------------------------------------------------------

int su_sfermion(double mql, double murpass, double mdrpass, double melpass, double merpass, double al, double at, double ab, double mu, double mst[], double msb[], double msl[], double msu[], double msd[], double mse[], double msn[])
{
    double mst1, mst2;
    mst1 = msttr1;
    mst2 = msttr2;

    double ct, st, mst1true, mst2true, mst22, mst12, mlrt, delt;
    double crll, crrr, crlr, rmb, rmt;
    double thettree, pscale;
    double b = beta;

    double mql1 = muq;
    double mur1 = mur;
    double mdr1 = mdr;
    double mel1 = mel;
    double mer1 = mer;

    double mur = murpass;
    double mdr = mdrpass;
    double mel = melpass;
    double mer = merpass;

    double tb = tan(b);

    // redefining s^2_w at EWSB scale:
    double cw = 1.0 / sqrt(1.0 + pow((g1ewsb / g2ewsb), 2));
    double sw = g1ewsb / g2ewsb * cw;
    double sw2ew = pow(sw, 2);
    sw2 = sw2ew;

    double pizz = 0.0;

    if (isnan(pizz) || pow(mz, 2) + pizz <= 0.0)
    {
        // Protections added non-pert or nan pb, uses tree-level values temporarily:
        pizz = 0.0;
        if (irge == irgmax)
            inonpert = -1;
    }

    double rmz = sqrt(pow(mz, 2) + pizz);
    double rmw = cw * rmz;
    double mzsave = mz;
    double mwsave = mw;
    mz = rmz;
    mw = rmw;

    // First two generations:  no mixing included

    // Up squarks:

    double mstl2 = pow(mql1, 2) + (0.50 - 2.0 / 3.0 * sw2) * pow(mz, 2) * cos(2.0 * b);
    double mstr2 = pow(mur1, 2) + 2.0 / 3.0 * sw2 * pow(mz, 2) * cos(2.0 * b);
    msu[1] = sqrt(mstl2);
    msu[2] = sqrt(mstr2);
    msu1 = msu[1]; // This variable is passing some value to oher functions
    msu2 = msu[2];

    //
    // Down squarks:

    double msbl2 = pow(mql1, 2) + (-0.50 + 1.0 / 3.0 * sw2) * pow(mz, 2) * cos(2.0 * b);
    double msbr2 = pow(mdr1, 2) - 1.0 / 3.0 * sw2 * pow(mz, 2) * cos(2.0 * b);
    msd[1] = sqrt(msbl2);
    msd[2] = sqrt(msbr2);
    msd1 = msd[1]; // This variable is passing some value to oher functions
    msd2 = msd[2];

    // Sleptons:

    double msel2 = pow(mel1, 2) + (-0.50 + sw2) * pow(mz, 2) * cos(2.0 * b);
    double mser2 = pow(mer1, 2) - sw2 * pow(mz, 2) * cos(2.0 * b);
    double msnl2 = pow(mel1, 2) + 0.50 * pow(mz, 2) * cos(2.0 * b);
    mse[1] = sqrt(msel2);
    mse[2] = sqrt(mser2);
    mse1 = mse[1]; // This variable is passing some value to oher functions
    mse2 = mse[2];

    if (msnl2 < 0.0)
    {
        msn[1] = 1.0;
        if (irge == irgmax)
            stnuerr = -1.0;
    }
    else
        msn[1] = sqrt(msnl2);

    msn1bp = msn[1];
    msn[2] = 1.0e+15;

    // Add radiative corrections to first gen. squarks

    if (isfrc == 1 && irge >= 2)
    {
        su_sqcr(alsewsb, m3, msu[1], &dmsu1);
        su_sqcr(alsewsb, m3, msu[2], &dmsu2);
        su_sqcr(alsewsb, m3, msd[1], &dmsd1);
        su_sqcr(alsewsb, m3, msd[2], &dmsd2);
        msu[1] = msu[1] + dmsu1;
        msu[2] = msu[2] + dmsu2;
        msd[1] = msd[1] + dmsd1;
        msd[2] = msd[2] + dmsd2;
    }

    //
    // Now the third generation sfermions:
    // stop masses/mixing
    // first some reinitializations:

    ifirst = 0;
    crll = 0.0;
    crrr = 0.0;
    crlr = 0.0;
    istflip = 0;

    // mb, mt, ml used in sfermion matrix elements should be running masses
    // at ewsb scale, including susy radiative corrections

    vd = sqrt(2 * pow(rmz, 2) / (pow(g1ewsb, 2) + pow(g2ewsb, 2)) / (1.0 + pow(tb, 2)));
    vu = vd * tb;
    rmb = ybewsb * vd;
    rmt = ytewsb * vu;
    rml = ytauewsb * vd;

l1:
    mstl2 = pow(mql, 2) + (0.50 - 2.0 / 3 * sw2ew) * pow(mz, 2) * cos(2 * b) + crll;
    mstr2 = pow(mur, 2) + 2.0 / 3 * sw2ew * pow(mz, 2) * cos(2 * b) + crrr;
    mlrt = at - mu / tb + crlr / rmt;

    delt = pow((mstl2 - mstr2), 2) + 4 * pow(rmt, 2) * pow(mlrt, 2);
    mst12 = pow(rmt, 2) + 0.50 * (mstl2 + mstr2 - sqrt(delt));
    mst22 = pow(rmt, 2) + 0.50 * (mstl2 + mstr2 + sqrt(delt));

    if (mst12 < 0.0)
    {
        //  tachyonic sfermion 1 mass
        mst[1] = 1.0;

        if (irge == irgmax)
            sterr = -1.0;
    }
    else
        mst[1] = sqrt(mst12);

    mst[2] = sqrt(mst22);
    thetout = atan(2 * rmt * mlrt / (mstl2 - mstr2)) / 2;

    if (ifirst == 1)
        mst1true = mst[1];
    if (ifirst == 2)
        mst2true = mst[2];

    if (ifirst == 3)
    {
        mst[1] = mst1true;
        mst[2] = mst2true;
    }

    ct = cos(thetout);
    st = sin(thetout);

    // Defining stop parameters at ewsb scale in bpmz conventions
    if (ifirst == 0)
    {
        mst1bp = sqrt(pow(ct, 2) * (pow(rmt, 2) + mstl2) + pow(st, 2) * (pow(rmt, 2) + mstr2) + 2 * ct * st * rmt * mlrt);
        mst2bp = sqrt(pow(st, 2) * (pow(rmt, 2) + mstl2) + pow(ct, 2) * (pow(rmt, 2) + mstr2) - 2 * ct * st * rmt * mlrt);

        if (isnan(mst1bp))
            mst1bp = 1.0; // Added protection
        if (isnan(mst2bp))
            mst2bp = 1.0;
        thetbp = thetout;
    }

    if (mstl2 > mstr2)
    {
        thetout = thetout + pi / 2;
        istflip = 1;
    }

    if (ifirst == 0)
    {
        // Save tree-level values for other uses:
        mst1 = mst[1];
        mst2 = mst[2];
        msttr1 = mst1;
        msttr2 = mst2;
        thettree = thetout;
    }

    // Adding rad. corr.
    if (isfrc == 1 && irge >= 2 && ifirst < 3)
    {
        ifirst = ifirst + 1;

        // Calculating stop rad. corr with 3 different momenta scales:
        if (ifirst == 1)
            pscale = mst1;
        if (ifirst == 2)
            pscale = mst2;
        if (ifirst == 3)
            pscale = sqrt(mst1 * mst2);

        su_stopcr(pscale, mu, at, ab, m3, &crll, &crlr, &crrr);
        goto l1;
    }

    ifirst = 0;

    // Sbottom masses/mixing
    isbflip = 0;
    msbl2 = pow(mql, 2) + (-0.50 + 1.0 / 3 * sw2ew) * pow(mz, 2) * cos(2 * b);
    msbr2 = pow(mdr, 2) - 1.0 / 3 * sw2ew * pow(mz, 2) * cos(2 * b);
    double mlrb = ab - mu * tb;
    double delb = pow((msbl2 - msbr2), 2) + 4 * pow(rmb, 2) * pow(mlrb, 2);
    double msb12 = pow(rmb, 2) + 0.50 * (msbl2 + msbr2 - sqrt(delb));
    double msb22 = pow(rmb, 2) + 0.50 * (msbl2 + msbr2 + sqrt(delb));

    if (msb12 < 0.0)
    {
        // Tachyonic sfermion mass
        msb[1] = 1.0;
        if (irge == irgmax)
            sberr = -1.0;
    }
    else
        msb[1] = sqrt(msb12);

    msb[2] = sqrt(msb22);

    thebout = atan(2 * rmb * mlrb / (msbl2 - msbr2)) / 2;
    double cb = cos(thebout);
    double sb = sin(thebout);
    double msb1bp, msb2bp;

    // Defining sbottom parameters at EWSB scale in bpmz conventions
    if (ifirst == 0)
    {
        msb1bp = sqrt(pow(cb, 2) * (pow(rmb, 2) + msbl2) + pow(sb, 2) * (pow(rmb, 2) + msbr2) + 2 * cb * sb * rmb * mlrb);
        msb2bp = sqrt(pow(sb, 2) * (pow(rmb, 2) + msbl2) + pow(cb, 2) * (pow(rmb, 2) + msbr2) - 2 * cb * sb * rmb * mlrb);
        if (isnan(msb1bp))
            msb1bp = 1.0; // Added protection
        if (isnan(msb2bp))
            msb2bp = 1.0;

        thebbp = thebout;
    }

    if (msbl2 > msbr2)
    {
        thebout = thebout + pi / 2;
        isbflip = 1;
    }
    double thebtree, msntau2, mlre, dele, mse12, mse22, sl;
    if (ifirst == 0)
    {
        // Save tree-level values for other uses:
        msb1 = msb[1];
        msb2 = msb[2];
        msbtr1 = msb1;
        msbtr2 = msb2;
        thebtree = thebout;
    }

    //
    // Add radiative corrections to sbottom  quarks
    double dmsb1, dmsb2;
    if (isfrc == 1 && irge >= 2)
    {
        su_sqcr(alsewsb, m3, msb[1], &dmsb1);
        su_sqcr(alsewsb, m3, msb[2], &dmsb2);

        msb[1] = msb[1] + dmsb1;
        msb[2] = msb[2] + dmsb2;
    }

    //
    // Stau masses/mixing

    msel2 = pow(mel, 2) + (-0.50 + sw2ew) * pow(mz, 2) * cos(2 * b);
    mser2 = pow(mer, 2) - sw2ew * pow(mz, 2) * cos(2 * b);
    msntau2 = pow(mel, 2) + 0.50 * pow(mz, 2) * cos(2 * b);
    mlre = al - mu * tb;
    dele = pow((msel2 - mser2), 2) + 4 * pow(rml, 2) * pow(mlre, 2);
    mse12 = pow(rml, 2) + 0.50 * (msel2 + mser2 - sqrt(dele));
    mse22 = pow(rml, 2) + 0.50 * (msel2 + mser2 + sqrt(dele));

    if (mse12 < 0.0)
    {
        // Tachyonic sfermion mass
        msl[1] = 1.0;
        if (irge == irgmax)
            stauerr = -1.0;
    }
    else
    {
        msl[1] = sqrt(mse12);
    }

    msl[2] = sqrt(mse22);
    thelout = atan(2 * rml * mlre / (msel2 - mser2)) / 2;
    cl = cos(thelout);
    sl = sin(thelout);

    if (msntau2 < 0.0)
    {
        stnuerr = -1.0;
        msn[3] = 1.0;
        if (irge == irgmax)
            return 0;
    }
    else
    {
        // Tau sneutrino:
        msn[3] = sqrt(msntau2);
    }

    msntau = msn[3];
    msn[4] = 1.0e+15;
    // Defining stau parameters at ewsb scale in bpmz conventions

    if (ifirst == 0)
    {
        msta1bp = sqrt(pow(cl, 2) * (pow(rml, 2) + msel2) + pow(sl, 2) * (pow(rml, 2) + mser2) + 2 * cl * sl * rml * mlre);
        msta2bp = sqrt(pow(sl, 2) * (pow(rml, 2) + msel2) + pow(cl, 2) * (pow(rml, 2) + mser2) - 2 * cl * sl * rml * mlre);
        if (isnan(msta1bp))
            msta1bp = 1.0; // Added protection
        if (isnan(msta2bp))
            msta2bp = 1.0;
        thelbp = thelout;
    }

    if (msel2 > mser2)
        thelout = thelout + pi / 2;

    // nb: for convenience msn(1--4) contains:
    // msn_{e,mu}[1],msn_{e,mu}[2], msn_{tau}[1],msn_{tau}[2]

    mz = mzsave;
    mw = mwsave;

    msbtr1 = msb1;
    msbtr2 = msb2;
    msttr1 = mst1;
    msttr2 = mst2;

    return 0;
}

//
// Calculates leading (chargino and neutralino loops) susy contributions to g_mu -2
//
// Input: mel,mer, al: relevant soft terms (i.e. 2d generation muon sector)
//  		mu, tb
//  		u,v,z, mn, mc1,mc2: chargino and neutralino masses and mixing matrices
// Ouptut:  gmuon, is a_mu = g_mu -2 in standard units
//-------------------------------------------------------------------------------------------------------------------------------------

#define fgm2a(x) (-(1.0 - 6 * x + 3 * pow(x, 2) + 2 * pow(x, 3) - 6 * pow(x, 2) * log(x)) / pow((1.0 - x), 4) / 6)
#define fgm2b(x) ((1.0 - pow(x, 2) + 2 * x * log(x)) / pow((1.0 - x), 3))
#define fgm2c(x) ((1.0 + 1.5 * x - 3 * pow(x, 2) + 0.50 * pow(x, 3) + 3 * x * log(x)) / pow((1.0 - x), 4) / 3)
#define fgm2d(x) (-3 * (1.0 - 4.0 / 3 * x + pow(x, 2) / 3 + 2.0 / 3 * log(x)) / pow((1.0 - x), 3))

void su_gminus2(double mel, double mer, double amu, double mu, double tb, double mc1, double mc2, double mn[]) //,double u[3][],double v[3][],double z[5][],double gmuon)
{
    double mc[3], msl[3], acl[3];
    ml = 0.105658357;
    mc[1] = mc1;
    mc[2] = mc2;
    double b = atan(tb);
    double cw = mw / mz;
    double sw = sqrt(1.0 - pow(cw, 2));
    sw2 = pow(sw, 2);

    //  calculation of the slepton masses and mixing
    double dt = cos(2.0 * b) * mz * mz;
    double msel2 = pow(mel, 2) + (-0.50 + sw2) * dt;
    double mser2 = pow(mer, 2) - sw2 * dt;
    double msnl2 = pow(mel, 2) + 0.50 * dt;
    double mlre = amu - mu * tb;
    double dele = pow((msel2 - mser2), 2) + 4.0 * pow((ml * mlre), 2);
    double mse12 = pow(ml, 2) + 0.50 * (msel2 + mser2 + sqrt(dele));
    double mse22 = pow(ml, 2) + 0.50 * (msel2 + mser2 - sqrt(dele));
    double msn = sqrt(msnl2);
    double thel = 0.50 * atan(2.0 * ml * mlre / (msel2 - mser2));
    double ccl = cos(thel);
    double ssl = sin(thel);

    // def. of mass eingenvalues in terms of angle:
    msl[1] = sqrt(pow(ccl, 2) * (pow(ml, 2) + msel2) + pow(ssl, 2) * (pow(ml, 2) + mser2) + 2 * ccl * ssl * ml * mlre);
    msl[2] = sqrt(pow(ssl, 2) * (pow(ml, 2) + msel2) + pow(ccl, 2) * (pow(ml, 2) + mser2) - 2 * ccl * ssl * ml * mlre);

    //  calculation of the chargino and neutralino couplings

    double rt2 = sqrt(2.0);
    double yuk = ml / (rt2 * mw * sw * cos(b));
    double xcl[3], xcr[3];
    xcl[1] = yuk * u[1][2];
    xcl[2] = yuk * u[2][2];
    xcr[1] = -v[1][1] / sw;
    xcr[2] = -v[2][1] / sw;

    //
    // i.e. for small mixing angle, slepton_1 should be mostly slepton_l

    double xnr[5][3], xnl[5][3], gau;
    for (int ii = 1; ii <= 4; ii++)
    {
        xnr[ii][1] = yuk * z[ii][3] * ccl + rt2 / cw * z[ii][1] * ssl;
        xnr[ii][2] = -yuk * z[ii][3] * ssl + rt2 / cw * z[ii][1] * ccl;
        gau = (z[ii][1] / cw + z[ii][2] / sw) / rt2;
        xnl[ii][1] = -yuk * z[ii][3] * ssl + gau * ccl;
        xnl[ii][2] = -yuk * z[ii][3] * ccl - gau * ssl;
    }
    double anl[5][3];
    for (int i = 1; i <= 4; i++)
        for (int j = 1; j <= 2; j++)
        {
            anl[i][j] = ml / pow(msl[j], 2) * (xnl[i][j] * xnl[i][j] + xnr[i][j] * xnr[i][j]) * fgm2a(pow(mn[i], 2) / pow(msl[j], 2));
            anl[i][j] += +mn[i] / pow(msl[j], 2) * xnl[i][j] * xnr[i][j] * fgm2b(pow(mn[i], 2) / pow(msl[j], 2));
        }

    double anltot = anl[1][1] + anl[2][1] + anl[3][1] + anl[4][1] + anl[1][2] + anl[2][2] + anl[3][2] + anl[4][2];

    for (int i = 1; i <= 2; i++)
        acl[i] = ml / pow(msn, 2) * (pow(xcl[i], 2) + pow(xcr[i], 2)) * fgm2c(pow(mc[i], 2) / pow(msn, 2)) + mc[i] / pow(msn, 2) * xcl[i] * xcr[i] * fgm2d(pow(mc[i], 2) / pow(msn, 2));

    double acltot = acl[1] + acl[2];
    gmuon = ml / 4.0 / pi / 137.0 * (anltot + acltot);

    // Include leading-log 2-loop qed correction

    double msusy = 0.2 * (msl[1] + msl[2] + msn + mc[1] + mc[2]);
    gmuon = gmuon * (1.0 - 4.0 / (pi * 137.0) * log(msusy / ml));
}

#undef fgm2a
#undef fgm2b
#undef fgm2c
#undef fgm2d

//
//  Input (all masses in GeV):
//    TGB       vev's ratio
//    CH_M      lightest chargino mass
//    AMU       MU parameter
//  Output:
//    AMG       gaugino mass
//    CHM(2)    chargino masses ( CHM(1) > CHM(2) )
//    U,V (2,2) chargino diagonalization matrices
//    IERR      1 (no solution for MU) 2 (divergent solution for MU)
//========================================================================

int chargino(double tgb, double ch_m, double amu, double amg, double chm[], int ierr)
{
    wm = 80.4190;
    ierr = 0;

    double sqrt2 = sqrt(2.0);
    double eps = 1.0e-2;
    //========================================================================
    //       definition of chargino's parameters
    //========================================================================
    double cbe = 1.0 / sqrt(1. + tgb * tgb);
    double sbe = tgb * cbe;
    double s2be = 2.0 * sbe * cbe;
    double c2be = cbe * cbe - sbe * sbe;

    //*********** chargino masses   *******************************************
    double amuq = amu * amu;
    double ch_m2 = pow(ch_m, 2);
    double wm2 = wm * wm;
    double xx = amuq - ch_m2;
    double sqrt_mu = xx * (xx + 2. * wm2) + pow((wm2 * s2be), 2);
    double amg1, amg2;

    if (sqrt_mu < 0.0)
    {
        ierr = 1;
        return 0;
    }
    if (fabs(xx) < pow(eps, 6))
    {
        ierr = 2;
        return 0;
    }
    if (fabs(xx) < eps)
    {
        amg1 = (wm2 * s2be - 2. * amuq / s2be) / (2. * amu);
        amg2 = 2. * amu * wm2 * s2be / xx;
    }
    else
    {
        amg1 = (wm2 * amu * s2be + ch_m * sqrt(sqrt_mu)) / xx;
        amg2 = (wm2 * amu * s2be - ch_m * sqrt(sqrt_mu)) / xx;
    }
    // nb here we can choose either the lighter or the heavier gaugino mass
    // according to whether we prefer a light or heavy heavier chargino
    //       however, choosing the max it is more likely that it will be positive
    // as m2 should be (see notes)

    amg = max(amg1, amg2);

    if (amg < 0.0)
    {
        ierr = 1;
        return 0;
    }

    double amgq = amg * amg;
    chm[2] = sqrt((amgq + amuq + 2. * wm2 - sqrt(pow((amgq - amuq), 2) + 4. * wm2 * (wm2 * c2be * c2be + amgq + amuq + 2. * amg * amu * s2be))) / 2.);
    double diff_ch = fabs(2.0 * (chm[2] - ch_m) / (chm[2] + ch_m));

    if (diff_ch > pow(eps, 2))
    {
        amg = (amg1 < amg2) ? amg1 : amg2;
        if (amg < 0.0)
        {
            ierr = 1;
            return 0;
        }
        amgq = pow(amg, 2);
        chm[2] = sqrt((amgq + amuq + 2. * wm2 - sqrt(pow((amgq - amuq), 2) + 4. * wm2 * (wm2 * c2be * c2be + amgq + amuq + 2. * amg * amu * s2be))) / 2.);
        diff_ch = fabs(2. * (chm[2] - ch_m) / (chm[2] + ch_m));
        if (diff_ch > eps * eps)
        {
            ierr = 1;
            return 0;
        }
    }

    chm[1] = sqrt((amgq + amuq + 2. * wm2 + sqrt(pow((amgq - amuq), 2) + 4. * wm2 * (wm2 * c2be * c2be + amgq + amuq + 2. * amg * amu * s2be))) / 2.);

    //************* u,v matrices *****************************************
    double hhh = (amg * amu - wm2 * s2be);
    double eep;
    if (hhh > 0.0)
        eep = 1.0;
    else
        eep = -1.0;

    double d_p = pow((amg + amu), 2) + 2. * wm2 * (1. - s2be);
    double d_m = pow((amg - amu), 2) + 2. * wm2 * (1. + s2be);
    double aa_p = (chm[1] + eep * chm[2]) * (amg + amu) / d_p;
    double aa_m = (chm[1] - eep * chm[2]) * (amg - amu) / d_m;
    double bb_p = wm * sqrt2 * (chm[1] + eep * chm[2]) * (sbe - cbe) / d_p;
    double bb_m = wm * sqrt2 * (chm[1] - eep * chm[2]) * (sbe + cbe) / d_m;

    double sqr = sqrt(pow((aa_p + aa_m), 2) + pow((bb_p + bb_m), 2));
    double cteta = sqr / 2.0;
    double steta = cteta * (aa_p - aa_m) / (bb_p + bb_m);
    double cphi = (aa_p + aa_m) / sqr;
    double sphi = (bb_p + bb_m) / sqr;

    v[1][1] = cphi;
    v[1][2] = sphi;
    v[2][1] = -sphi * eep;
    v[2][2] = cphi * eep;
    u[1][1] = cteta;
    u[1][2] = steta;
    u[2][1] = -steta;

    // Check
    // ------------------------------------------------
    double cck1 = u[1][1] * amg + u[1][2] * sqrt2 * wm * cbe;
    double cck2 = u[2][1] * amg + u[2][2] * sqrt2 * wm * cbe;
    double cbk1 = u[1][2] * amu + u[1][1] * sqrt2 * wm * sbe;
    double cbk2 = u[2][2] * amu + u[2][1] * sqrt2 * wm * sbe;
    double c1 = cck1 * v[1][1] + cbk1 * v[1][2] - chm[1];
    double c2 = cck1 * v[2][1] + cbk1 * v[2][2];
    double c3 = cck2 * v[1][1] + cbk2 * v[1][2];
    double c4 = cck2 * v[2][1] + cbk2 * v[2][2] - chm[2];
    u[2][2] = cteta;

    if (fabs(c1) > eps)
        printf("\ncheck failed c1= %e", c1);
    if (fabs(c2) > eps)
        printf("\ncheck failed c2= %e", c2);
    if (fabs(c3) > eps)
        printf("\ncheck failed c3= %e", c3);
    if (fabs(c4) > eps)
        printf("\ncheck failed c4= %e", c4);
    return 0;
}

//
//
// ------------------------------------------------------------------------------

double su_fr(double x,double y)
{
	return x+y-2*x*y/(x-y)*log(x/y);
	}

void su_delrho(double mt,double gmst[],double gmsb[],double gmstau[],double msn,double thetat,double thetab,double thel,double *drho) 
{
	//
	//   calculates leading one-loop SUSY delta_rho contributions of 3rd gen
	// sfermions (plus leading two-loop QCD contributions) 
	//  INPUT: MT, gmst(2), gmsb(2),gmstau(2),msn: top,stop,sbottom,
	//  stau, stau neutrino masses and stop, sbottom, stau mixing angles
	//  OUTPUT: drho = rho-1 
	//----------------------------------------------------------------------------
	double alph = alpha;
	double ct = cos(thetat);
	double st = sin(thetat);
	double cb = cos(thetab);
	double sb = sin(thetab);
	double ctau = cos(thel);
	double stau = sin(thel);
	double cta2 = pow(ctau,2);
	double sta2 = pow(stau,2);
	double ct2 = pow(ct,2);
	double st2 = pow(st,2);
	double cb2 = pow(cb,2);
	double sb2 = pow(sb,2);
	double mt1 = pow(gmst[1],2);
	double mt2 = pow(gmst[2],2);
	double mb1 = pow(gmsb[1],2);
	double mb2 = pow(gmsb[2],2);
	double mta1 = pow(gmstau[1],2);
	double mta2 = pow(gmstau[2],2);
	
	double drhotb  =   ct2*(cb2*su_fr(mt1,mb1)+sb2*su_fr(mt1,mb2)); 
	drhotb += + st2*(cb2*su_fr(mt2,mb1)+sb2*su_fr(mt2,mb2)); 
	drhotb += - ct2*st2*su_fr(mt1,mt2)-cb2*sb2*su_fr(mb1,mb2);

	double drhotau= -cta2*sta2*su_fr(mta1,mta2)+cta2*su_fr(mta1,msn*msn) + sta2*su_fr(mta2,msn*msn);
	*drho = 3*drhotb*(1.0 +2*0.12/3/pi*(1.0+pow(pi,2)/3))+drhotau;
	*drho = gf/(8*pi*pi* sqrt(2.0))*(*drho);
	}


//
// ---------------------------------------------------------------------------------------------------

//
//  B -> sp in supersymmetry: large contributions beyond the leading order
//  arXiv:hep-ph/0009337v2  (Eq. 5)
// --------------------------------------------------------------------------------------
double ff1(double x)
{

    if (fabs(x - 1.0) > 1.0e-3)
    {
        double d = 1. / pow((x - 1), 3);
        double dg = log(x) / pow((x - 1), 4);
        return x * (7 - 5 * x - 8 * x * x) * d / 24.0 + x * x * (3 * x - 2) * dg / 4.0;
    }
    else
        return -5.0 / 48.0;
}

//
//  arXiv:hep-ph/0009337v2  (Eq. 6)
// --------------------------------------------------------------------------------------
double fg1(double x)
{
    if (fabs(x - 1.0) > 1.0e-2)
    {
        double d = 1. / pow((x - 1), 3);
        double dg = log(x) / pow((x - 1), 4);
        return x * (2 + 5 * x - x * x) * d / 8.0 - 3 * x * x * dg / 4.0;
    }
    else
        return -1.0 / 16.0;
}

//
//  arXiv:hep-ph/0009337v2  (Eq. 7)
// --------------------------------------------------------------------------------------
double ff2(double x)
{
    if (fabs(x - 1.0) > 1.0e-2)
    {
        double d = 1.0 / pow((x - 1), 2);
        double dg = log(x) / pow((x - 1), 3);
        return x * (3.0 - 5.0 * x) * d / 12.0 + x * (3.0 * x - 2.0) * dg / 6.0;
    }
    else
    {
        return -7.0 / 36.0;
    }
}

//
//  arXiv:hep-ph/0009337v2  (Eq. 8)
// --------------------------------------------------------------------------------------
double fg2(double x)
{
    if (fabs(x - 1.) > 1.0e-2)
    {
        double d = 1.0 / pow((x - 1), 2);
        double dg = log(x) / pow((x - 1), 3);
        return x * (3.0 - x) * d / 4.0 - x * dg / 2.0;
    }
    else
        return -1.0 / 6.0;
}

//
//  arXiv:hep-ph/0009337v2  (Eq. 21)
// --------------------------------------------------------------------------------------
double ff3(double x)
{
    return 2.0 / 3.0 * (1 - 1.0 / x) * ff1(x) + ff2(x) + 23.0 / 36.0;
}

//
//  arXiv:hep-ph/0009337v2  (Eq. 22)
// --------------------------------------------------------------------------------------
double fg3(double x)
{
    return 2.0 / 3.0 * (1 - 1.0 / x) * fg1(x) + fg2(x) + 1.0 / 3.0;
}

double ff(double x)
{
    double b[12];
    double z = -log(1.0 - x);
    b[1] = -.50;
    b[2] = 1.0 / 6.0;
    b[3] = 0.0;
    b[4] = -1.0 / 30.0;
    b[5] = 0.0;
    b[6] = 1.0 / 42.0;
    b[7] = 0.0;
    b[8] = -1.0 / 30.0;
    b[9] = 0.0;
    b[10] = 5.0 / 66.0;
    b[11] = 0.0;
    b[12] = -691.0 / 2730.0;
    double cc = z;
    double sum = z;
    for (int i = 1; i <= 12; i++)
    {
        cc = cc * z / (i + 1.0);
        sum = sum + b[i] * cc;
    }
    return sum;
}

//
//  Calculating the Spence Function
//  https://arxiv.org/pdf/hep-ph/9710335  (Used in maany equations Li_2(x))
// --------------------------------------------------------------------------------------
double sp(double x)
{
    double sp1;
    if (x > -1.0 && x < 0.50)
        sp1 = ff(x);
    else if (x > 0.50 && x < 1.0)
        sp1 = -ff(1.0 - x) + pi * pi / 6.0 - log(x) * log(1.0 - x);
    else if (x < -1.0)
        sp1 = -ff(1.0 / x) - pi * pi / 6 - 0.50 * pow(log(-x), 2);
    else
        printf("\nError in dilog : %e", x);

    return sp1;
}

//
// https://arxiv.org/pdf/hep-ph/9710335 (Eq. 38)
// -------------------------------------------------------------------------------------
double fgsm(double x)
{
    if (fabs(x - 1.0) > 1.0e-3)
    {
        double x2 = x * x;
        double x3 = x2 * x;
        double x4 = x3 * x;
        double x5 = x4 * x;
        double d4 = pow((x - 1), 4);
        double d5 = d4 * (x - 1);
        double spx = sp(1 - 1.0 / x);
        double xl = log(x);
        double xl2 = xl * xl;
        double fgsm1 = (-4 * x4 + 40 * x3 + 41 * x2 + x) * spx / (6 * d4) + (-17 * x3 - 31 * x2) * xl2 / (2 * d5) + (-210 * x5 + 1086 * x4 + 4893 * x3 + 2857 * x2 - 1994 * x + 280) * xl / (216 * d5);
        return fgsm1 + (737 * x4 - 14102 * x3 - 28209 * x2 + 610 * x - 508) / (1296 * d4);
    }
    else
        return -9821.0 / 12960.0;
}

//
// https://arxiv.org/pdf/hep-ph/9710335 (Eq. 37)
// -------------------------------------------------------------------------------------
double ffsm(double x)
{
    if (fabs(x - 1.0) > 1.0e-3)
    {
        double x2 = x * x;
        double x3 = x2 * x;
        double x4 = x3 * x;
        double x5 = x4 * x;
        double d4 = pow((x - 1), 4);
        double d5 = d4 * (x - 1);
        double spx = sp(1 - 1. / x);
        double xl = log(x);
        double xl2 = xl * xl;
        double ffsm1 = (-16 * x4 - 122 * x3 + 80 * x2 - 8 * x) * spx / (9 * d4) + (6 * x4 + 46 * x3 - 28 * x2) * xl2 / (3 * d5) + (-102 * x5 - 588 * x4 - 2262 * x3 + 3244 * x2 - 1364 * x + 208) * xl / (81 * d5);
        return ffsm1 + (1646 * x4 + 12205 * x3 - 10740 * x2 + 2509 * x - 436) / (486 * d4);
    }
    else
        return -3451.0 / 9720.0;
}

//
// https://arxiv.org/pdf/hep-ph/9710335 (Eq. 60).
// So terms are slightly changed based on SUSPECT. Need some references for those terms
// -------------------------------------------------------------------------------------
double ffh(double wth)
{
    double th = wth;
    if (fabs(wth - 1) < 1.0e-3)
        th = 1.0 + sign(1.0, (wth - 1.0)) * 1.0e-3;
    double temp = 4.0 * (7.0 * th - 13.0 * th * th + 2.0 * pow(th, 3)) / (3.0 * pow((-1.0 + th), 3));
    temp += 16.0 * (-3.0 * th + 7.0 * th * th - 2.0 * pow(th, 3)) * sp(1.0 - 1.0 / th) / (9.0 * pow((-1.0 + th), 3));
    temp += 8.0 * (-3.0 * th - th * th + 12.0 * pow(th, 3) - 2.0 * pow(th, 4)) * log(th) / (9.0 * pow((-1.0 + th), 4));
    temp += 4.0 * (8.0 * th - 14.0 * th * th - 3.0 * pow(th, 3)) * pow(log(th), 2) / (9.0 * pow((-1.0 + th), 4));
    double ffh1 = ad * au * (temp);
    temp = (893.0 * th - 5706.0 * th * th + 7785.0 * pow(th, 3) - 1244 * pow(th, 4)) / (486.0 * pow((-1.0 + th), 4));
    temp += 2.0 * (18.0 * th * th - 37.0 * pow(th, 3) + 8.0 * pow(th, 4)) * sp(1.0 - 1.0 / th) / (9.0 * pow((-1.0 + th), 4));
    temp += 2.0 * (-56.0 * th + 266.0 * pow(th, 2) - 183.0 * pow(th, 3) - 192.0 * pow(th, 4) + 21.0 * pow(th, 5)) * log(th) / (81.0 * pow((-1.0 + th), 5));
    temp += 2.0 * (-14.0 * th * th + 23.0 * pow(th, 3) + 3.0 * pow(th, 4)) * pow(log(th), 2) / (9.0 * pow((-1.0 + th), 5));
    ffh1 += au * au * (temp);
    return ffh1;
}

//
// https://arxiv.org/pdf/hep-ph/9710335 (Eq. 62).
// So terms are slightly changed based on SUSPECT. Need some references for those terms
// -------------------------------------------------------------------------------------
double fgh(double wth)
{
    double th = wth;
    if (fabs(wth - 1.0) < 1.0e-3)
        th = 1.0 + sign(1.0, (wth - 1.0)) * 1.0e-3;

    double temp = (143.0 * th - 44.0 * th * th + 29.0 * pow(th, 3)) / (8.0 * pow((-1.0 + th), 3));
    temp += (-36.0 * th + 25.0 * th * th - 17.0 * pow(th, 3)) * sp(1.0 - 1.0 / th) / (6.0 * pow((-1.0 + th), 3));
    temp += (-3.0 * th - 187.0 * th * th + 12.0 * pow(th, 3) - 14.0 * pow(th, 4)) * log(th) / (12.0 * pow((-1.0 + th), 4));
    temp += (19.0 * th + 17.0 * th * th) * pow(log(th), 2) / (3.0 * pow((-1.0 + th), 4));
    double fgh1 = ad * au * (temp);

    temp = (1226.0 * th - 18423.0 * th * th + 7866.0 * pow(th, 3) - 4493.0 * pow(th, 4)) / (1296.0 * pow((-1.0 + th), 4));
    temp += (30.0 * pow(th, 2) - 17.0 * pow(th, 3) + 13.0 * pow(th, 4)) * sp(1.0 - 1.0 / th) / (6.0 * pow((-1.0 + th), 4));
    temp += (-238.0 * th + 847.0 * pow(th, 2) + 1335.0 * pow(th, 3) + 318.0 * pow(th, 4) + 42.0 * pow(th, 5)) * log(th) / (216.0 * pow((-1.0 + th), 5));
    temp += -(31.0 * th * th + 17.0 * pow(th, 3)) * pow(log(th), 2.0) / (6.0 * pow((-1.0 + th), 5));

    fgh1 += au * au * (temp);
    return fgh1;
}

//
// https://arxiv.org/pdf/hep-ph/9710335 (Eq. 63).
// -------------------------------------------------------------------------------------

double fghlog(double wth)
{
    double th = wth;

    if (fabs(wth - 1) < 1.0e-3)
        th = 1.0 + sign(1.0, (wth - 1.0)) * 1.0e-3;

    double temp = ((-81.0 * th + 16.0 * pow(th, 2) - 7.0 * pow(th, 3)) / (6.0 * pow((-1.0 + th), 3)) + (19.0 * th + 17.0 * pow(th, 2)) * log(th) / (3.0 * pow((-1.0 + th), 4)));
    double fghlog1 = ad * au * temp;
    temp = (38.0 * th + 261.0 * pow(th, 2) - 18.0 * pow(th, 3) + 7.0 * pow(th, 4)) / (36.0 * pow((-1.0 + th), 4)) - (31.0 * pow(th, 2) + 17.0 * pow(th, 3)) * log(th) / (6.0 * pow((-1.0 + th), 5));
    fghlog1 += pow(au, 2) * temp;
    return -fghlog1;
}

//
// https://arxiv.org/pdf/hep-ph/9710335 (Eq. 61).
// -------------------------------------------------------------------------------------
double ffhlog(double wth)
{
    double th = wth;
    if (fabs(wth - 1) < 1.0e-3)
        th = 1.0 + sign(1.0, (wth - 1.0)) * 1.0e-3;

    double temp = 2.0 * (-21.0 * th + 47.0 * pow(th, 2) - 8.0 * pow(th, 3)) / (9.0 * pow((-1.0 + th), 3));
    temp += 4.0 * (8.0 * th - 14.0 * pow(th, 2) - 3.0 * pow(th, 3)) * log(th) / (9.0 * pow((-1.0 + th), 4));

    double ffhlog1 = ad * au * temp;

    temp = (31.0 * th + 18.0 * pow(th, 2) - 135.0 * pow(th, 3) + 14.0 * pow(th, 4)) / (27.0 * pow((-1.0 + th), 4));
    temp += 2.0 * (-14.0 * pow(th, 2) + 23.0 * pow(th, 3) + 3.0 * pow(th, 4)) * log(th) / (9.0 * pow((-1.0 + th), 5));
    ffhlog1 += pow(au, 2) * temp;

    return -ffhlog1;
}

//
// https://arxiv.org/pdf/hep-ph/9710335 (Eq. 64).
// -------------------------------------------------------------------------------------
double eh(double x)
{
    if (fabs(x - 1.0) > 1.0e-3)
    {
        double d = 1.0 / pow((x - 1.0), 3);
        double dg = log(x) / pow((x - 1.0), 4);
        return x * (16 - 29 * x + 7 * x * x) * d / 36.0 + x * (3 * x - 2) * dg / 6.0;
    }
    else
        return 1.0 / 8.0;
}

//
// https://arxiv.org/pdf/hep-ph/0009337 (Eq. 12).
// -------------------------------------------------------------------------------------
double h2(double x, double y)
{
    if (fabs(x - y) > 1.0e-2)
    {
        if (fabs(x - 1.0) < 1.0e-2)
            return 1.0 / (-1 + y) - y * log(y) / pow((-1 + y), 2);
        else
        {
            if (fabs(y - 1.0) < 1.0e-2)
                return 1.0 / (-1 + x) - x * log(x) / pow((-1 + x), 2);
            else
                return x * log(x) / ((1 - x) * (x - y)) + y * log(y) / ((1 - y) * (-x + y));
        }
    }
    else
    {
        if (fabs(x - 1.0) < 1.0e-2)
            return -1 / 2.0;
        else
            return 1.0 / (1 - x) + log(x) / pow((-1 + x), 2);
    }
    return 0;
}

//
// https://arxiv.org/pdf/hep-ph/9710335 (Eq. 39)
// ---------------------------------------------------------------------------------
double ffsmlog(double x)
{
    if (fabs(x - 1.0) > 1.0e-3)
    {
        double diff = (-x * (18 * x * x + 30 * x - 24) * log(x) + 47 * x * x * x - 63 * x * x + 9 * x + 7) / 24.0 / pow((x - 1), 5);
        return 8 * x * diff + 16.0 / 3.0 * ff1(x) - 16.0 / 9.0 * fg1(x) + 208.0 / 81.0;
    }
    else
        return 1369.0 / 810.0;
}

//
// https://arxiv.org/pdf/hep-ph/9710335 (Eq. 40)
// ---------------------------------------------------------------------------------
double fgsmlog(double x)
{
    if (fabs(x - 1.0) > 1.0e-3)
    {
        double diff = (6 * x * (1 + x) * log(x) - x * x * x - 9 * x * x + 9 * x + 1) / 4.0 / pow((x - 1), 5);
        return 8 * x * diff + 14. / 3.0 * fg1(x) + 35. / 27.0;
    }
    else
        return 869.0 / 1080.0;
}

//
// https://arxiv.org/pdf/hep-ph/9710335 (Eq. 36)
// ---------------------------------------------------------------------------------
double esm(double x)
{
    if (fabs(x - 1.0) > 1.0e-3)
    {
        double d = 1.0 / pow((x - 1), 3);
        double dg = log(x) / pow((x - 1), 4);
        double esm1 = x * (-18 + 11 * x + x * x) * d / 12.0 + x * x * (15 - 16 * x + 4 * x * x) * dg / 6.0;
        return esm1 - 2.0 / 3.0 * log(x);
    }
    else
        return 43.0 / 72.0;
}

//
// https://arxiv.org/pdf/hep-ph/9612313 (Eq. 23, Eq. 24)
// -------------------------------------------------------------------------------------
double asf(double x)
{
    //     CALCOLA ALPHA_S (X)
    //     NB it uses 5 flavors for x<mt and 6 if x>=mt

    double fn = 5.0;
    double b0 = 11.0 - 2 * fn / 3.0;
    double b1 = 102.0 - 38. * fn / 3.0;
    double vvv = 1 - b0 * asc / (2 * pi) * log(zm / x);
    double asf1, ast1;

    if (ioc == 0)
        asf1 = asc / vvv;
    else
        asf1 = asc / vvv * (1 - b1 / b0 * asc / (4 * pi * vvv) * log(vvv));

    if (x > tc)
    {
        vvv = 1 - b0 * asc / (2 * pi) * log(zm / tc);

        if (ioc == 0)
            ast1 = asc / vvv;
        else
            ast1 = asc / vvv * (1 - b1 / b0 * asc / (4 * pi * vvv) * log(vvv));

        double b0t = b0 - 2.0 / 3.0;
        double b1t = b1 - 38.0 / 3.0;
        vvv = 1 - b0t * ast1 / (2 * pi) * log(tc / x);

        if (ioc == 0)
            asf1 = ast1 / vvv;
        else
            asf1 = ast1 / vvv * (1 - b1t / b0t * ast1 / (4 * pi * vvv) * log(vvv));
    }
    return asf1;
}

//
// https://arxiv.org/pdf/hep-ph/9710335 (Eq 32)
// ------------------------------------------------------------------------------
double runt(double x)
{
    // CALCOLA M_top^{running} (X)

    double fn;
    if (x <= tc)
        fn = 5.0;
    else
        fn = 6.0;

    b0 = 11.0 - 2 * fn / 3.0;
    b1 = 102.0 - 38.0 * fn / 3.0;
    double g0 = 8.0;
    double g1 = 404.0 / 3.0 - 40 * fn / 9.0;
    double asx = asf(x);
    double ast = asf(tc);
    double rrr = tc * pow((asx / ast), (g0 / (2 * b0)));

    if (ioc == 0)
        return rrr;

    double corr = 1 + ast * g0 / (4 * pi * 2 * b0) * (-b1 / b0 + g1 / g0) * (asx / ast - 1);
    // this is the relation between mpole/mt(mpole)
    // this is the relation between mpole/mtrun(mtrun)
    double pol = 1 + 4 * ast / (3 * pi) + 8.2430 * pow((ast / pi), 2);
    // this is the 1loop order result
    return rrr * corr / pol;
}

double h3(double x)
{
    if (fabs(x - 1.0) > 1.0e-2)
        return 2.0 * x * log(x) / (1 - x);
    return -2.0;
}

double h4(double x)
{
    if (fabs(x - 1.0) > 1.0e-2)
        return 2.0 * log(x) / (1 - x);
    return -2.0;
}

double gg7(double x)
{
    if (fabs(x - 1.0) > 1.0e-2)
    {
        double d = 1.0 / pow((x - 1), 3);
        double dg = log(x) / pow((x - 1), 4);
        return -(23 - 67 * x + 50 * x * x) * d / 36.0 + x * (2 - 7 * x + 6 * x * x) * dg / 6.0;
    }
    else
        return 3.0 / 8.0;
}

double gg8(double x)
{
    if (fabs(x - 1.0) > 1.0e-3)
    {
        double d = 1.0 / pow((x - 1), 3);
        double dg = log(x) / pow((x - 1), 4);
        return -(4 - 5 * x - 5 * x * x) * d / 12.0 + x * (1 - 2 * x) * dg / 2.0;
    }
    else
        return 1.0 / 8.0;
}

double ffs1(double wth)
{
    double x = wth;
    if (fabs(wth - 1.0) < 1.0e-3)
        x = 1.0 + sign(1.0, (wth - 1.0)) * 1.0e-3;
    double ffs11 = (-4.0 * (-5.0 + 3.0 * x)) / (9.0 * pow((-1.0 + x), 2));
    ffs11 += (16.0 * (3.0 - 7.0 * x) * x * sp(1.0 - 1.0 / x)) / (9.0 * pow((-1.0 + x), 3));
    ffs11 += (4.0 * (4.0 - 30.0 * x + 40.0 * x * x) * log(x)) / (9.0 * pow((-1.0 + x), 3));
    ffs11 += (16.0 * (1.0 - 3.0 * x) * x * pow(log(x), 2)) / (9.0 * pow((-1.0 + x), 3));
    return ffs11;
}

double ffs5(double wth)
{
    double x = wth;
    if (fabs(wth - 1.0) < 1.0e-3)
        x = 1.0 + sign(1.0, (wth - 1.0)) * 1.0e-3;
    double ffs51 = (-85.0 + 347.0 * x - 526.0 * x * x) / (243.0 * pow((-1.0 + x), 3));
    ffs51 += (4.0 * x * (-8.0 + 13.0 * x + 6.0 * x * x) * sp(1.0 - 1.0 / x)) / (9.0 * pow((-1.0 + x), 4));
    ffs51 += (4.0 * (-20.0 + 126.0 * x - 144.0 * x * x - 39.0 * x * x * x) * log(x)) / (81.0 * pow((-1.0 + x), 4));
    ffs51 += (2.0 * x * (-10.0 + 21. * x) * pow(log(x), 2)) / (9.0 * pow((-1.0 + x), 4));
    return ffs51;
}

double fgs1(double wth)
{
    double x = wth;
    if (fabs(wth - 1.0) < 1.0e-3)
        x = 1.0 + sign(1.0, (wth - 1.0)) * 1.0e-3;
    double fgs11 = (61.0 - 39.0 * x) / (12.0 * pow((-1.0 + x), 2));
    fgs11 += (4.0 * x * (3.0 + 4.0 * x) * sp(1.0 - 1.0 / x)) / (3. * pow((-1.0 + x), 3));
    fgs11 += ((7.0 - 60.0 * x - 14.0 * x * x) * log(x)) / (6.0 * pow((-1.0 + x), 3));
    fgs11 += (14.0 * x * pow(log(x), 2)) / (3.0 * pow((-1.0 + x), 3));
    return fgs11;
}

double fgs5(double wth)
{
    double x = wth;
    if (fabs(wth - 1.0) < 1.0e-3)
        x = 1.0 + sign(1.0, (wth - 1.0)) * 1.0e-3;

    double fgs51 = (-1210.0 + 437.0 * x + 1427.0 * x * x) / (648.0 * pow((-1.0 + x), 3));
    fgs51 += -(x * (49.0 + 46.0 * x + 9.0 * x * x) * sp(1.0 - 1.0 / x)) / (12.0 * pow((-1.0 + x), 4));
    fgs51 += +((-85.0 + 603.0 * x + 387.0 * x * x - 78.0 * x * x * x) * log(x)) / (108.0 * pow((-1.0 + x), 4));
    fgs51 += -(13.0 * x * pow(log(x), 2)) / (3.0 * pow((-1.0 + x), 4));
    return fgs51;
}

double h1fun(double x)
{
    if (fabs(x - 1.0) > 1.0e-2)
        return 1.0 / (1 - x) + (2 * x - x * x) * log(x) / pow((1 - x), 2);
    return -1.0 / 2.0;
}

//
// NB there is a misprint in Gabrielli-Giudice a factor2
// moreover this corresponds to E(1/x) in GG notation
//---------------------------------------------------------------

double echi(double x)
{
    if (fabs(x - 1.0) > 1.0e-3)
    {
        double d = 1.0 / pow((x - 1.0), 3);
        double dg = log(x) / pow((x - 1.0), 4);
        return x * (11.0 - 7.0 * x + 2.0 * x * x) * d / 18.0 - x * dg / 3.0;
    }
    else
        return 1.0 / 12.0;
}

//
// in the following function we have used h1(0)=1 and h2(1,1)=-1/2
//--------------------------------------------------------------------

double utd(double sqd1)
{
    double t1 = pow((stophc / glc), 2);
    double d1 = pow((sqd1 / glc), 2);
    double sst = sin(pi * xstopc);
    double cst = cos(pi * xstopc);

    double tlr = sst * cst * (pow(stophc, 2) - pow(stoplc, 2));
    return -h1fun(d1) / 2.0 + (cst * cst * h1fun(t1) + sst * sst) / 2.0 - (tlr / tc) * h4(t1) / glc + (tlr / tc) * (cst * cst * h4(d1) - sst * sst) / glc;
}

//*************************************************************************
// for the checks h1[x_]:= 1/(1-x) + (2x-x^2) Log[x]/(1-x)^2
// h4[x_]:= 2 Log[x]/(1-x)
//*************************************************************************

double udd(double sqd1, double sqd2) // to be checked
{
    double t1 = pow((stophc / glc), 2);
    double d1 = pow((sqd1 / glc), 2);
    double d2 = pow((sqd2 / glc), 2);
    double tb = 1.0 / au;
    double sst = sin(pi * xstopc);
    double cst = cos(pi * xstopc);
    //  here we set Ad= At   ! corrected with true Ad now (jlk):
    double tlr = sst * cst * (pow(stophc, 2) - pow(stoplc, 2));
    //        aad= tlr/tc + amuc/tb
    double aad = ad;

    double udd1 = h1fun(d1) / 2.0 - (cst * cst * h1fun(t1) + sst * sst) / 2.0 - 2.0 * (aad - amuc * tb) / glc * h2(d1, d2);
    udd1 += 2.0 * (aad - amuc * tb) / glc * (cst * cst * h2(t1, d2) + sst * sst * h4(d2) / 2.0);

    return udd1;
}

//
// in the paper this is H_td
//----------------------------------------------------------------------
double ftd(double sqd1) // checked
{
    double t1 = pow((stophc / glc), 2);
    double d1 = pow((sqd1 / glc), 2);
    double sst = sin(pi * xstopc);
    double cst = cos(pi * xstopc);
    double tlr = sin(xstopc * pi) * cos(xstopc * pi) * (pow(stophc, 2) - pow(stoplc, 2));

    return -h1fun(d1) / 2.0 + (pow(cst, 2) * h1fun(t1) + pow(sst, 2)) / 2.0 - (tlr / tc) * h4(t1) / glc + (tlr / tc + amuc * (au + 1.0 / au)) * (pow(cst, 2) * h4(d1) - pow(sst, 2)) / glc;
}

//
// In the paper this is H_d
//--------------------------------------------------------------
double fdd(double sqd1, double sqd2) // to be checked
{
    double t1 = pow((stophc / glc), 2);
    double d1 = pow((sqd1 / glc), 2);
    double d2 = pow((sqd2 / glc), 2);
    double sst = sin(pi * xstopc);
    double cst = cos(pi * xstopc);
    double tb = 1.0 / au;
    // here we set Ad= At  ! corrected with true Ad coupling now (jlk)
    double tlr = sin(xstopc * pi) * cos(xstopc * pi) * (pow(stophc, 2) - pow(stoplc, 2));
    double aad = ad;

    double fdd1 = h1fun(d1) / 2.0 - (pow(cst, 2) * h1fun(t1) + pow(sst, 2)) / 2.0 - 2.0 * (aad - amuc * tb) / glc * h2(d1, d2);
    fdd1 += 2.0 * (aad + amuc / tb) / glc * (pow(cst, 2) * h2(t1, d2) + pow(sst, 2) * h4(d2) / 2.0);

    return fdd1;
}

//
//  Here there are functions which are useful for the derivation of
//  the ckm elements from fit
//-----------------------------------------------------------------------
double aab(double x)
{
    if (fabs(x - 1.0) > 1.0e-3)
        return (-4 + 15 * x - 12 * x * x + x * x * x + 6 * x * x * log(x)) / (4.0 * pow((-1 + x), 3));
    else
        return 3.0 / 4.0;
}

double gboxh(double x)
{
    if (fabs(x - 1.0) > 1.0e-3)
        return (-1 + x * x - 2 * x * log(x)) / pow((-1 + x), 3);
    else
        return 1.0 / 3.0;
}

double fpbox(double x, double y)
{
    if (fabs(x - y) > 1.0e-3)
        return (-1 + x - log(x)) / (pow((-1.0 + x), 2) * (-x + y)) + (x * log(x) / (-1 + x) + y * log(y) / (1 - y)) / pow((x - y), 2);
    else
        return (-1 + x * x - 2 * x * log(x)) / (2 * pow((-1 + x), 3) * x);
}

double gpbox(double x, double y)
{
    if (fabs(x - y) > 1.0e-3)
    {
        if (fabs(y - 1) > 1.0e-3)
            return (3 - 4 * x + x * x + 4 * x * log(x) - 2 * x * x * log(x)) / (2.0 * pow((-1 + x), 2) * (-x + y)) - (3 * (-x + y) / 2.0 + x * x * log(x) / (-1 + x) + y * y * log(y) / (1 - y)) / pow((x - y), 2);
        else
            return (-1 + x * x - 2 * x * log(x)) / pow((-1 + x), 3);
    }
    else
    {
        if (fabs(x - 1) > 1.0e-3)
            return (3 / 2 - 2 * x + x * x / 2 + log(x)) / pow((-1 + x), 3);
        else
            return 1.0 / 3.0;
    }
}

double deltat2(double sqd1) // checked
{
    double t1 = pow((stophc / glc), 2);
    double d1 = pow((sqd1 / glc), 2);

    double tlr = sin(xstopc * pi) * cos(xstopc * pi) * (pow(stophc, 2) - pow(stoplc, 2));
    return 3.0 + h3(d1) + h1fun(t1) / 2.0 - (tlr / tc) * h4(t1) / glc;
}

double delta2(double sqd1, double sqd2) //  checked
{

    double d1 = pow((sqd1 / glc), 2);
    double d2 = pow((sqd2 / glc), 2);
    double tb = 1.0 / au;
    double tlr = sin(xstopc * pi) * cos(xstopc * pi) * (pow(stophc, 2) - pow(stoplc, 2));
    // aad= tlr/tc + amuc/tb  ! corrected with true Ad now (jlk):
    double aad = ad;

    double ctt = cos(xstopc * pi) / sin(xstopc * pi);
    return h3(d2) * (1 - ctt * tc * glc / pow(sqd2, 2)) + 5. / 2.0 + (h1fun(d1) + h1fun(d2)) / 2.0 - 2.0 * (aad - amuc * tb) * h2(d1, d2) / glc;
    // left over from a test?
    // delta2= - 2.d0*(-amuc*tb)*h2(d1,d2)/glc
}

double ww(double sq1, double sq2) // checked
{
    double x = pow((sq1 / glc), 2);
    double y = pow((sq2 / glc), 2);
    if (fabs(x - y) > 1.0e-3)
    {
        if (fabs(x - 1.0) < 1.0e-3)
            return (1 + 5 * y) / (2 * (1 - y)) + y * (2 + y) * log(y) / pow((-1 + y), 2);
        else
        {
            if (fabs(y - 1.0) < 1.0e-3)
                return (1 + 5 * x) / (2 * (1 - x)) + x * (2 + x) * log(x) / pow((-1 + x), 2);
            else
                return (x + y - 2 * x * y) / ((-1 + x) * (-1 + y)) + (x * x * x - 2 * x * y + x * x * y) * log(x) / (pow((-1 + x), 2) * (x - y)) + (2 * x * y - x * y * y - y * y * y) * log(y) / ((x - y) * pow((-1 + y), 2));
        }
    }
    else
        return 0.0;
}

double delta1(double sq) // checked
{

    //	double glc = 0;
    double x = pow((sq / glc), 2);
    double sst = sin(pi * xstopc);
    double cst = cos(pi * xstopc);
    return -3.0 / 4.0 - h1fun(x) / 2.0 + 2.0 * sst * cst * tc / glc;
    // The last term has been neglected according to the paper
}

//
//
//       Input:
//       imod        SM (0), Two Higgs (1), SUSY (2)
//       io          LO (io=0) NLO (io=1)
//       nlosusy     susy coeff at LO (=0) or NLO (=1)
//       ihv         if nlosusy=1: 1= heavy sqk/gluino case (new paper),
//                                 0= light stop_R case (old paper)
//       scw         matching scale
//       as          alpha_s (M_Z)
//       t           pole top mass
//       h           Higgs mass
//       tanb        tan(beta)
//       stopl       lighter stop mass
//       stoph       heavier stop mass
//       xstop       stop mixing angle (in units of PI)
//       sb1,sb2     sbottom masses
//       xsbot       sbottom mixing angle (in units of PI)
//       sqk         common left-squark mass
//       gl         mgluino
//       aat        At soft parameter
//       abo        Ab soft parameter
//       amu         MU parameter
//       ch_m         lightest chargino mass
//       u,v         chargino U,V matrices
//       Output:
//       c70         C7 (LO)
//       c71         C7 (NLO)
//       c80         C8 (LO)
//       c81         C8 (NLO)
//       ee          E(x)
//       ierr        OK (0), manca soluzione per M (1,2)
//       R           R_delta = Delta_susy/Delta_SM per bbbar boxes
// -------------------------------------------------------------------------------------------

void matching(int imod, int io, int nlosusy, int ihv, double scw, double as, double t, double h, double tanb, double stopl, double stoph, double xstop, double sb1, double sb2, double xsbot, double sqk, double gl, double aat, double abo, double amu, double chm[], double *c70l, double *c80l, double *c71l, double *c81l, double *eel, double *bboxl, double *ierrl)
{

    double c70 = *c70l;
    double c71 = *c71l;
    double c80 = *c80l;
    double c81 = *c81l;
    double ee = *eel;
    double bbox = *bboxl;
    double ierr = *ierrl;

    double nem[5], nn[5][5];

    // PARAMETRI DI INPUT
    double wm = 80.4190;
    double zm = 91.18670;

    // ALTRI PARAMETRI
    double pi2 = pi * pi;
    double sq2 = sqrt(2.0);
    double r23 = 2.0 / 3.0;
    double sinb = tanb / sqrt(1 + tanb * tanb);
    double cosb = sinb / tanb;
    double asc = as;
    double tc = t;
    double stoplc = stopl;
    double stophc = stoph;
    double xstopc = xstop;
    double xsbotc = xsbot;
    double amuc = amu;
    double sqkc = sqk;
    double glc = gl;
    int ioc = io;

    // SCELTA MODELLO HIGGS
    au = 1.0 / tanb;
    ad = -tanb; // Modello II
    // ad=1./tanb       // Modello I

    // MASSE RUNNING A SCW
    // topr= runt(scw)
    double topr = runm(scw, 6);

    double xt = pow((topr / wm), 2);
    if (imod >= 1)
        yt = pow((topr / h), 2);

    double str[3], sbr[3], tmix[3][3], bmix[3][3];
    double xsqc[3], xstc[3][3];
    if (imod >= 2)
    {
        str[1] = stoph; // no running mass for stops
        str[2] = stopl;

        tmix[1][1] = cos(xstop * pi);
        tmix[1][2] = sin(xstop * pi);
        tmix[2][1] = -tmix[1][2];
        tmix[2][2] = tmix[1][1];

        // cgg add sbottom
        sbr[1] = sb2; // no running mass for sbottoms
        sbr[2] = sb1;
        bmix[1][1] = cos(xsbot * pi);
        bmix[1][2] = sin(xsbot * pi);
        bmix[2][1] = -bmix[1][2];
        bmix[2][2] = bmix[1][1];

        ierr = 0;

        for (int i = 1; i <= 2; i++)
        {
            xsqc[i] = pow((sqk / chm[i]), 2);
            for (int k = 1; k <= 2; k++)
                xstc[k][i] = pow((str[k] / chm[i]), 2);
        }
    }

    // COEFF WILSON A SCW (LO)
    c70 = ff1(xt);
    c80 = fg1(xt);

    // ==> NNB electroweak corrections only for the calculation of c70b
    // be careful that it is not used both here and in the routine bsg
    // c70= c70 *(0.974d0) !-0.005d0)
    // c80= c80 *(0.993d0) ! -0.005d0)

    if (imod > 1)
    {
        c70 = c70 + au * au / 3. * ff1(yt) - au * ad * ff2(yt);
        c80 = c80 + au * au / 3. * fg1(yt) - au * ad * fg2(yt);
    }

    double accd[3], yy[3], acc[3][3];
    double c70s, c80s;
    double sint, cost, sbmx, cbmx;
    double c7ll, c8ll, c7sg, c8sg, c7sh, c8sh;
    if (imod >= 2)
    {
        sint = tmix[1][2];
        cost = tmix[1][1];

        sbmx = bmix[1][2];
        cbmx = bmix[1][1];

        for (int j = 1; j <= 2; j++)
        {
            accd[j] = u[j][2] * wm / (sq2 * cosb * chm[j]);
            yy[j] = v[j][2] * tmix[2][2] * topr / (sq2 * wm * sinb);
            for (int k = 1; k <= 2; k++)
                acc[j][k] = v[j][1] * tmix[k][1] - v[j][2] * tmix[k][2] * topr / (sq2 * wm * sinb);
        }

        c70s = 0.0;
        c80s = 0.0;
        for (int j = 1; j <= 2; j++)
        {
            // NB for the limit of heavy squarks but stop2 comment next4 lines + do k=2,2
            c70s = c70s + r23 * pow(v[j][1], 2) * pow((wm / sqk), 2) * ff1(xsqc[j]);
            c80s = c80s + r23 * pow(v[j][1], 2) * pow((wm / sqk), 2) * fg1(xsqc[j]);
            c70s = c70s + accd[j] * v[j][1] * ff3(xsqc[j]);
            c80s = c80s + accd[j] * v[j][1] * fg3(xsqc[j]);

            for (int k = 1; k <= 2; k++) // k is the stop index
            {
                c70s = c70s - r23 * pow(acc[j][k], 2) * pow((wm / str[k]), 2) * ff1(xstc[k][j]);
                c80s = c80s - r23 * pow(acc[j][k], 2) * pow((wm / str[k]), 2) * fg1(xstc[k][j]);
                c70s = c70s - accd[j] * acc[j][k] * tmix[k][1] * ff3(xstc[k][j]);
                c80s = c80s - accd[j] * acc[j][k] * tmix[k][1] * fg3(xstc[k][j]);
            }
        }
        c70 = c70 + c70s;
        c80 = c80 + c80s;
    }

    // COEFF WILSON A SCW (NLO)
    double ees;
    double g7chi1, g8chi1, g7chi2, g8chi2, g7chi3, g8chi3, g7chi4, g8chi4, c71s, c81s;
    double cf, c71sw, c81sw, c71sp, c81sp, c71sh, c81sh;
    double amsusy, dmb, xx1, xx2, u1, u2, y11, y12, y21, y22, yukt;
    double epsb, epsbew, epsbp, epst;
    double akk, eta;
    double c7imp, c8imp, aacc[3][3];
    double deltaboxs;
    int isb;

    if (io == 1)
    {
        // gg - Standard Model part
        c71 = ffsm(xt) + ffsmlog(xt) * 2 * log(scw / wm); // Eq 35 (https://arxiv.org/pdf/hep-ph/9710335)
        c81 = fgsm(xt) + fgsmlog(xt) * 2 * log(scw / wm); // Eq 35 (https://arxiv.org/pdf/hep-ph/9710335)
        ee = esm(xt);                                     // Eq 36 (https://arxiv.org/pdf/hep-ph/9710335)

        if (imod >= 1)
        {
            // gg - Charged Higgs two-doublets part
            c71 = c71 + ffh(yt) + ffhlog(yt) * 2 * log(scw / h); // Eq 59 (https://arxiv.org/pdf/hep-ph/9710335)
            c81 = c81 + fgh(yt) + fghlog(yt) * 2 * log(scw / h); // Eq 59 (https://arxiv.org/pdf/hep-ph/9710335)
            ee = ee + au * au * eh(yt);
        }

        if (imod >= 2 && nlosusy == 1)
        {
            // gg - SUSY part: depends on which approximation is valid ...
            c71s = 0.0;
            c81s = 0.0;
            ees = 0.0;
            for (int j = 1; j <= 2; j++)
                ees = ees + pow(acc[j][2], 2) * pow((wm / str[2]), 2) * echi(xstc[2][j]);

            if (ihv == 0) // This is the old light stop_R case
            {
                for (int j = 1; j <= 2; j++)
                {
                    // ees = ees + acc(j,2)**2*(wm/str(2))**2*echi(xstc(2,j))
                    // notice that in the following sqk is used for all the squarks of
                    // light flavors
                    // the functions g7,8chi1,2 are not exactly the same as in the paper
                    // but they are combinations

                    g7chi1 = -8.0 / 9.0 * ff1(xstc[2][j]) * (-1.0 + 6.0 * log(glc / chm[j]) + 2.0 * delta1(sqk)) - 4 / 9.0 * echi(xstc[2][j]) + xstc[2][j] * ffs5(xstc[2][j]) + 32. * (fg1(xstc[2][j]) - 3. * ff1(xstc[2][j])) * 2 * log(scw / chm[j]) / 27.0;
                    g8chi1 = -8.0 / 9.0 * fg1(xstc[2][j]) * (-1.0 + 2. * delta1(sqk) + 6.0 * log(glc / chm[j])) - 1 / 6.0 * echi(xstc[2][j]) + xstc[2][j] * fgs5(xstc[2][j]) - 28. * fg1(xstc[2][j]) * 2 * log(scw / chm[j]) / 9.0;
                    g7chi2 = -4.0 / 3.0 * ff3(xstc[2][j]) * (2. * delta1(sqk) + delta2(sqk, sqk) - 2.0) - ffs1(xstc[2][j]) + 16. * (fg3(xstc[2][j]) - 3. * ff3(xstc[2][j])) * 2 * log(scw / chm[j]) / 9.0;
                    g8chi2 = -4.0 / 3.0 * fg3(xstc[2][j]) * (2. * delta1(sqk) + delta2(sqk, sqk) - 2.0) - fgs1(xstc[2][j]) - 14. * fg3(xstc[2][j]) * 2 * log(scw / chm[j]) / 3.0;
                    g7chi3 = 8.0 / 9.0 * ff1(xstc[2][j]) * (12. * log(scw / glc) + 2.0 * deltat2(sqk) - 2.0);
                    g8chi3 = 8.0 / 9.0 * fg1(xstc[2][j]) * (12. * log(scw / glc) + 2.0 * deltat2(sqk) - 2.0);
                    g7chi4 = 4.0 / 3.0 * ff3(xstc[2][j]) * (6. * log(scw / glc) + deltat2(sqk) - 1.0);
                    g8chi4 = 4.0 / 3.0 * fg3(xstc[2][j]) * (6. * log(scw / glc) + deltat2(sqk) - 1.0);
                    c71s = c71s + pow(acc[j][2], 2) * pow((wm / str[2]), 2) * g7chi1 + accd[j] * tmix[2][1] * acc[j][2] * g7chi2 + acc[j][2] * yy[j] * pow((wm / str[2]), 2) * g7chi3 + accd[j] * yy[j] * tmix[2][1] * g7chi4;
                    c81s = c81s + pow(acc[j][2], 2) * pow((wm / str[2]), 2) * g8chi1 + accd[j] * tmix[2][1] * acc[j][2] * g8chi2 + acc[j][2] * yy[j] * pow((wm / str[2]), 2) * g8chi3 + accd[j] * yy[j] * tmix[2][1] * g8chi4;
                }
                // additional contributions from shifts of SM and 2HDM
                // NB in ..sw we use ww(sqk,0)=1/2 for sqk=mgl
                cf = 4.0 / 3.0;
                c71sw = cf * (pow(tmix[1][1], 2) * ww(sqk, stophc) + pow(tmix[1][2], 2) / 2.0) * gg7(xt);
                c81sw = cf * (pow(tmix[1][1], 2) * ww(sqk, stophc) + pow(tmix[1][2], 2) / 2.0) * gg8(xt);
                c71sp = cf * (2 * utd(sqk) * ff1(xt) / 3.0 - (utd(sqk) + udd(sqk, sqk)) * ff2(xt));
                c81sp = cf * (2 * utd(sqk) * fg1(xt) / 3.0 - (utd(sqk) + udd(sqk, sqk)) * fg2(xt));
                c71sh = cf * (2 * pow(au, 2 / 3.0) * ftd(sqk) * ff1(yt) + (ftd(sqk) + fdd(sqk, sqk)) * ff2(yt));
                c81sh = cf * (2 * pow(au, 2 / 3.0) * ftd(sqk) * fg1(yt) + (ftd(sqk) + fdd(sqk, sqk)) * fg2(yt));

                // if you want to have only LO susy coeff comment the following 3 lines
                c71 = c71 + c71s + c71sw + c71sp + c71sh;
                c81 = c81 + c81s + c81sw + c81sp + c81sh;
            }

            else if (ihv == 1) // This is the new heavy colored particles case
            {
                // LARGE LOGS ONLY (ex large tanb scenario)
                // NNB Msusy is set HERE!!!
                // amsusy= 1200.d0 ! msusy, to be set equal to average?
                amsusy = sqk; // set equal to average squark mass

                //  NB in dmb H2=-1/2 this is true only for degenerate sbottoms!!
                dmb = asf(amsusy) / pi / 3.0 * amuc * tanb / glc;

                //  notation as in the new paper
                xx1 = pow((sb1 / glc), 2);
                xx2 = pow((sb2 / glc), 2);
                u1 = pow((stophc / glc), 2);
                u2 = pow((stoplc / glc), 2);
                y11 = pow((stophc / chm[1]), 2);
                y12 = pow((stophc / chm[2]), 2);
                y21 = pow((stoplc / chm[1]), 2);
                y22 = pow((stoplc / chm[2]), 2);
                yukt = 0.650 * runt(amsusy) / wm / sq2;

                epsb = -asf(amsusy) * 2.0 / 3.0 / pi * amuc / glc * h2(xx1, xx2);
                epsbew = -pow((yukt / 4 / pi), 2) * aat / chm[1] * u[1][2] * v[1][2] * h2(y11, y21) - pow((yukt / 4 / pi), 2) * aat / chm[2] * u[2][2] * v[2][2] * h2(y12, y22);

                epsbew = 0.0; //.............................. Check may be there is some gonjamil

                isb = 1;
                if (isb == 0)
                    epsbp = -asf(amsusy) * 2.0 / 3.0 / pi * amuc / glc * (pow(cost, 2) * h2(u1, xx2) + pow(sint, 2) * h2(u2, xx2));
                else if (isb == 1)
                    epsbp = -asf(amsusy) * 2.0 / 3.0 / pi * amuc / glc * (pow(cost, 2) * (h2(u1, xx2) * pow(cbmx, 2) + h2(u1, xx1) * pow(sbmx, 2)) + pow(sint, 2) * (h2(u2, xx2) * pow(cbmx, 2) + h2(u2, xx1) * pow(sbmx, 2)));

                if (isb == 0)
                    epst = -asf(amsusy) * 2.0 / 3.0 / pi * amuc / glc * (pow(cost, 2) * h2(u2, xx1) + pow(sint, 2) * h2(u1, xx1));
                else if (isb == 1)
                    epst = -asf(amsusy) * 2.0 / 3.0 / pi * amuc / glc * (pow(cost, 2) * (h2(u2, xx1) * pow(cbmx, 2) + h2(u2, xx2) * pow(sbmx, 2)) + pow(sint, 2) * (h2(u1, xx1) * pow(cbmx, 2) + h2(u1, xx2) * pow(sbmx, 2)));

                akk = 1 / (1 + (epsb + epsbew) * tanb);
                eta = asf(amsusy) / asf(scw);

                // first redefine the couplings;

                for (int j = 1; j <= 2; j++)
                    for (int k = 1; k <= 2; k++)
                        aacc[j][k] = v[j][1] * tmix[k][1] - v[j][2] * tmix[k][2] * runt(amsusy) / (sq2 * wm * sinb);

                c7imp = 0.0;
                c8imp = 0.0;

                for (int j = 1; j <= 2; j++)
                {
                    // == first we write the improved LO coeff
                    c7imp = c7imp + r23 * pow(v[j][1], 2) * pow((wm / sqk), 2) * ff1(xsqc[j]);
                    c8imp = c8imp + r23 * pow(v[j][1], 2) * pow((wm / sqk), 2) * fg1(xsqc[j]);
                    c7imp = c7imp + akk * accd[j] * v[j][1] * ff3(xsqc[j]);
                    c8imp = c8imp + akk * accd[j] * v[j][1] * fg3(xsqc[j]);
                    for (int k = 1; k <= 2; k++) // k is the stop index
                    {
                        c7imp = c7imp - r23 * pow(aacc[j][k], 2) * pow((wm / str[k]), 2) * ff1(xstc[k][j]);
                        c8imp = c8imp - r23 * pow(aacc[j][k], 2) * pow((wm / str[k]), 2) * fg1(xstc[k][j]);
                        c7imp = c7imp - akk * accd[j] * aacc[j][k] * tmix[k][1] * ff3(xstc[k][j]);
                        c8imp = c8imp - akk * accd[j] * aacc[j][k] * tmix[k][1] * fg3(xstc[k][j]);
                    }
                }

                // ==   First, the LL from the running of the ops
                //        c7ll= -32.d0/3.d0*(c70s-c80s/3.d0)*log(msy/scw)
                //         c8ll= -28.d0/3.d0*c80s*log(msy/scw)
                // == the same but resummed

                c7ll = pow(eta, (16.0 / 21.0)) * c7imp + 8.0 / 3.0 * (pow(eta, (14.0 / 21.0)) - pow(eta, (16.0 / 21.0))) * c8imp;
                c8ll = pow(eta, (14.0 / 21.0)) * c8imp;
                c7sg = tanb * (epsb - epsbp) * akk * ff2(xt) / (asf(amsusy) / 4.0 / pi);
                c8sg = tanb * (epsb - epsbp) * akk * fg2(xt) / (asf(amsusy) / 4.0 / pi);

                c7sh = 8.0 / 3.0 * amuc * tanb / glc * ff2(yt) * (-0.50 + pow(cost, 2) * h2(pow((stophc / glc), 2), 0.990) + sint * sint * h2(pow((stoplc / glc), 2), 0.990));
                c8sh = 8.0 / 3.0 * amuc * tanb / glc * fg2(yt) * (-0.50 + pow(cost, 2) * h2(pow((stophc / glc), 2), 0.990) + sint * sint * h2(pow((stoplc / glc), 2), 0.990));

                // if you want to have only LO susy coeff
                // if(nlosusy.eq.1)then
                c71s = c7sg; //+ c7ll+c7sh
                c81s = c8sg; //+c8ll+c8sh
                c70 = c70 + c7ll - c70s;
                c80 = c80 + c8ll - c80s;
                c71 = c71 + c71s; //+c71sw+c71sp+c71sh these are only for light stop2
                c81 = c81 + c81s; //+c81sw+c81sp+c81sh
            }
            // switch on ihv
            ee = ee + ees;
        }
    }

    // the ratio R_delta for the B-Bbar mixing and epsilon
    double deltabox = aab(xt); // SM contribution

    if (imod >= 1) // 2HDM contribution
        deltabox = deltabox + pow(au, 4) * yt * gboxh(yt) / 4.0 + 2 * au * au * xt * (fpbox(xt, xt / yt) + gpbox(xt, xt / yt) / 4.0);

    if (imod >= 2) // here we include only the light stop
    {
        deltaboxs = 0.0;
        for (int ibox = 1; ibox <= 2; ibox++)
            for (int jbox = 1; jbox <= 2; jbox++)
                deltaboxs = deltaboxs + pow(acc[jbox][2], 2) * pow(acc[ibox][2], 2) * (wm * wm / pow(chm[jbox], 2)) / xt * gpbox(xstc[2][jbox], (pow(chm[ibox], 2) / pow(chm[jbox], 2)));
        deltabox = deltabox + deltaboxs;
    }

    bbox = deltabox / aab(xt);

    *c70l = c70;
    *c71l = c71;
    *c80l = c80;
    *c81l = c81;
    *eel = ee;
    *bboxl = bbox;
    *ierrl = ierr;
}

double gre(double t)
{
    // PARTE REALE DI G(T)
    if (t < 4.0)
        return -2 * pow((atan(sqrt(t / (4 - t)))), 2);

    double pi2 = pi * pi;
    return -pi2 / 2 + 2 * pow((log((sqrt(t) + sqrt(t - 4)) / 2.0)), 2);
}

double gim(double t)
{
    // PARTE IMMAGINARIA DI G(T)
    if (t < 4.0)
        return 0.0;

    return -2 * pi * (log((sqrt(t) + sqrt(t - 4)) / 2.0));
}

double func1(double t)
{
    double z = zfunz;
    return (1 - z * t) * (pow((gre(t) / t + 0.5), 2) + pow((gim(t) / t), 2));
}

double func2(double t)
{
    double z = zfunz;
    return pow((1 - z * t), 2) * (pow((gre(t) / t + 0.5), 2) + pow((gim(t) / t), 2));
}

double func3(double t)
{
    double z = zfunz;
    return gre(t) + t / 2.0;
}

double func4(double t)
{
    double z = zfunz;
    return (1 - z * t) * (gre(t) + t / 2.0);
}

//
//
// -------------------------------------------------------------------------------

double gauss1(double (*f)(double), double a, double b, double eps)
{
    double cons = 1.0e-25;

    double w[] = {
        0.0000000000000000000000000000,
        0.1012285362903762591525313543,
        0.2223810344533744705443559944,
        0.3137066458778872873379622020,
        0.3626837833783619829651504493,
        0.0271524594117540948517805725,
        0.0622535239386478928628438370,
        0.0951585116824927848099251076,
        0.1246289712555338720524762822,
        0.1495959888165767320815017305,
        0.1691565193950025381893120790,
        0.1826034150449235888667636680,
        0.1894506104550684962853967232};

    double x[] = {
        0.0000000000000000000000000000,
        0.9602898564975362316835608686,
        0.7966664774136267395915539365,
        0.5255324099163289858177390492,
        0.1834346424956498049394761424,
        0.9894009349916499325961541735,
        0.9445750230732325760779884155,
        0.8656312023878317438804678977,
        0.7554044083550030338951011948,
        0.6178762444026437484466717640,
        0.4580167776572273863424194430,
        0.2816035507792589132304605015,
        0.0950125098376374401853193354};

    double delta = cons * fabs(a - b);
    double gauss11 = 0;
    double aa = a;
    double y, bb, c1, c2, s8, s16, u;

    for (;;)
    {
        y = b - aa;
        if (fabs(y) <= delta)
            return gauss11;
    l2:
        bb = aa + y;
        c1 = 0.5 * (aa + bb);
        c2 = c1 - aa;
        s8 = 0.0;
        s16 = 0.0;
        for (int i = 1; i <= 4; i++)
        {
            u = x[i] * c2;
            s8 = s8 + w[i] * (f(c1 + u) + f(c1 - u));
        }

        for (int i = 5; i <= 12; i++)
        {
            u = x[i] * c2;
            s16 = s16 + w[i] * (f(c1 + u) + f(c1 - u));
        }
        s8 = s8 * c2;
        s16 = s16 * c2;
        if (fabs(s16 - s8) > eps * (1.0 + fabs(s16)))
            break;

        gauss11 += s16;
        aa = bb;
    }

    y = 0.5 * y;
    if (fabs(y) > delta)
        goto l2;

    printf("\nTOO HIGH ACCURACY REQUIRED");
    gauss11 = 0.0;
    return gauss11;
}

  
//
// ---------------- Beginning of code --------------------------------------------
//       The matching routine calculates the Wilson coefficients for a given model.
//       The bsg routine calculates the branching ratio starting from the
//       Wilson coefficients. Version 1/6/98, bug fixed.
//       Our electroweak corrections are implemented only for mh = 100. (7/2000)

//       Input:
//       as          alpha_s (M_Z)
//       t           top mass [enters only in the calculation of alpha_s(scw)]
//       bc          pole bottom mass - pole charm mass
//       rcb         pole charm mass / pole bottom mass
//       aeinv       1/alpha(muw)
//       scw         matching scale
//       scb         low-energy scale
//       scl         scale for semileptonic
//       vkm         |V_ts V_tb / V_cb |^2
//       bsl         semileptonic BR
//       deltp       cut on photon energy
//       io          LO (io=0) NLO (io=1)
//       c70,c71,c80,
//       c81,ee      matching conditions
//       bbox        ratio box/box_SM for CKM matrix elemts
//       Output:
//       br          BR(b->s gamma)
// ------------------------------------------------------------------------------------------------------

void su_bsg(double as, double t, double bc, double rcb, double aeinv, double scw, double scb, double scl, double vkm, double bsl, double deltp, double io, double c70, double c71, double c80, double c81, double ee, double bbox, double *brs)
{

    double br = *brs;
    double aa[9], eee[9], ff[9], gg[9], hh[9], clead[9], fnum[9][9], bb[9];
    // PARAMETRI DI INPUT
    wm = 80.419;
    zm = 91.187;
    double rbs = 50.0;
    double delt = deltp;
    double hlam = 0.12; // lambda_2 (HQET parameter)
    int ityp = 1;       // calcola f_ij tramite integrali
    // ityp=2      // valori numerici di f_ij per delt=0.99 o 0.9 e sqrt(z)=0.29

    // EW corrections for mh=100 ??
    int iew = 1;

    // ALTRI PARAMETRI
    asc = as;
    ioc = io;
    tc = t;
    double z = rcb * rcb;
    zfunz = z;
    double ae = 1.0 / aeinv;

    double pi2 = pi * pi;

    if (c70 == 0)
    {
        br = 0;
        return;
    }

    // MAGIC NUMBERS
    aa[1] = 14.0 / 23.0;
    aa[2] = 16.0 / 23.0;
    aa[3] = 6.0 / 23.0;
    aa[4] = -12.0 / 23.0;
    aa[5] = 0.4086;
    aa[6] = -0.4230;
    aa[7] = -0.8994;
    aa[8] = 0.1456;
    eee[1] = 4661194.0 / 816831.0;
    eee[2] = -8516.0 / 2217.0;
    eee[3] = 0.0;
    eee[4] = 0.0;
    eee[5] = -1.9043;
    eee[6] = -0.1008;
    eee[7] = 0.1216;
    eee[8] = 0.0183;
    ff[1] = -17.3023;
    ff[2] = 8.5027;
    ff[3] = 4.5508;
    ff[4] = 0.7519;
    ff[5] = 2.0040;
    ff[6] = 0.7476;
    ff[7] = -0.5385;
    ff[8] = 0.0914;
    gg[1] = 14.8088;
    gg[2] = -10.8090;
    gg[3] = -0.8740;
    gg[4] = 0.4218;
    gg[5] = -2.9347;
    gg[6] = 0.3971;
    gg[7] = 0.1600;
    gg[8] = 0.0225;
    hh[1] = 626126.0 / 272277.0;
    hh[2] = -56281.0 / 51730.0;
    hh[3] = -3.0 / 7.0;
    hh[4] = -1.0 / 14.0;
    hh[5] = -0.6494;
    hh[6] = -0.0380;
    hh[7] = -0.0186;
    hh[8] = -0.0057;
    bb[1] = 0.5784;
    bb[2] = -0.3921;
    bb[3] = -0.1429;
    bb[4] = 0.04762;
    bb[5] = -0.1275;
    bb[6] = 0.03174;
    bb[7] = 0.007811;
    bb[8] = -0.003101;
    
    // CALCOLO ALPHA
    double asscw = asf(scw);
    double asscb = asf(scb);
    double eta = asscw / asscb;

    // COEFF WILSON A SCB
    double c10b = -pow(eta, (-12. / 23.)) + pow(eta, (6. / 23.));
    double c20b = (pow(eta, (-12. / 23.)) + 2 * pow(eta, (6. / 23.))) / 3.;
    double sumh = 0.0;

    for (int i = 1; i <= 8; i++)
        sumh = sumh + hh[i] * pow(eta, aa[i]);

    double c70b = pow(eta, (16.0 / 23.0)) * c70 + 8. / 3. * (pow(eta, (14. / 23.)) - pow(eta, (16. / 23.))) * c80 + sumh;

    int iwrongsign = 0;
    if (c70b < 0.0)
        iwrongsign = 0; // Sign as in SM
    else
        iwrongsign = 1; // Sign has been flipped: trouble w/ b --> sll!

    // addition of sign of C7 = sign of SM constraint (equiv b->s l+l- sign):
    // electroweak corrections only for the calculation of c70b
    double c70t, c80t, c70bew;
    if (iew == 1)
    {
        c70t = c70 * (0.974); //-0.005d0)
        c80t = c80 * (0.993); // -0.005d0)
                              //   c70(mb) including ew corrections; last term from ew from mixings
        c70bew = pow(eta, (16. / 23.)) * c70t + 8. / 3. * (pow(eta, (14.0 / 23.0)) - pow(eta, (16. / 23.))) * c80t + sumh + 0.00054;
    }
    else
        c70bew = c70b;

    double c80b = (c80 + 313063.0 / 363036.0) * pow(eta, (14.0 / 23.0)) - 0.9135 * pow(eta, 0.4086) + 0.0873 * pow(eta, (-0.4230)) - 0.0571 * pow(eta, (-0.8994)) + 0.0209 * pow(eta, 0.1456);
    double sumg, c71b;
    if (io == 1)
    {
        sumg = 0.0;
        for (int i = 1; i <= 8; i++)
        {
            // old     sumg=sumg+(eee(i)*eta*ee+
            //     +   ff(i)+gg(i)*eta)*eta**aa(i)
            sumg = sumg + (eee[i] * eta * (ee + 4.0 / 3.0 * log(scw / wm)) + ff[i] + (gg[i] + 12 * bb[i] * log(scw / wm)) * eta) * pow(eta, aa[i]);
        }
        c71b = pow(eta, (39. / 23)) * c71 + 8. / 3 * (pow(eta, (37. / 23)) - pow(eta, (39. / 23))) * c81;
        c71b = c71b + (297664. / 14283. * pow(eta, (16. / 23.)) - 7164416. / 357075. * pow(eta, (14. / 23.)) + 256868. / 14283. * pow(eta, (37. / 23.)) - 6698884. / 357075. * pow(eta, (39. / 23.))) * c80;
        c71b = c71b + 37208. / 4761. * (pow(eta, (39. / 23.)) - pow(eta, (16. / 23.))) * c70 + sumg;
    }
    clead[1] = c10b;
    clead[2] = c20b;
    clead[3] = -pow(eta, (-12.0 / 23.0)) / 27. + 2 * pow(eta, (6. / 23.)) / 63. - 0.0659 * pow(eta, 0.4086) + 0.0595 * pow(eta, (-0.423)) - 0.0218 * pow(eta, (-0.8994)) + 0.0335 * pow(eta, 0.1456);
    clead[4] = -pow(eta, (-12.0 / 23.0)) / 9. + pow(eta, (6. / 23.)) / 21. + 0.0237 * pow(eta, 0.4086) - 0.0173 * pow(eta, (-0.423)) - 0.1336 * pow(eta, (-0.8994)) - 0.0316 * pow(eta, 0.1456);
    clead[5] = pow(eta, (-12.0 / 23.0)) / 108. - pow(eta, (6. / 23.)) / 126. + 0.0094 * pow(eta, 0.4086) - 0.01 * pow(eta, (-0.423)) + 0.0010 * pow(eta, (-0.8994)) - 0.0017 * pow(eta, 0.1456);
    clead[6] = -pow(eta, (-12.0 / 23.0)) / 36. - pow(eta, (6. / 23.)) / 84. + 0.0108 * pow(eta, 0.4086) + 0.0163 * pow(eta, (-0.423)) + 0.0103 * pow(eta, (-0.8994)) + 0.0023 * pow(eta, 0.1456);
    clead[7] = c70b;
    clead[8] = c80b;

    // CORREZIONI HQET
    double bm2 = pow((bc / (1 - sqrt(z))), 2);
    double bm = sqrt(bm2); // bottom quark mass
    double cm2 = z * bm2;  // charm squared mass
    double fkin = 1 - 8 * z + 8 * pow(z, 3) - pow(z, 4) - 12 * pow(z, 2) * log(z);
    double hqb = -6 * hlam * (1 - pow((1 - z), 4) / fkin);
    double hqc = hlam * (c10b / 54.0 - c20b / 9.0) / c70b; // corrected according to NPB version
    double hqet = 1 + hqb / bm2 + hqc / cm2;

    // CALCOLO LEADING
    // vkm from input is the result of a SM fit here we change to other models
    double box = 1.0 + (1.04884 - vkm) * (1.0 - 1.0 / sqrt(bbox)) / vkm;
    // box=1.d0   !in case we don't want that
    double deltg, cb, csl;
    // Czarnecki-Marciano QED correction
    double alpha0 = 1.0 / 137.036;
    double cqed0 = alpha0 / ae * (1 - 2.0 * ae * log(zm / bm) / pi) * (1 - ae / pi * (104.0 / 243.0 - 8.0 * c70 / 9.0) * log(wm / bm) / fabs(c70b));
    double delt2, delt3, deltg1, fff77, fff88, fff78, eps, est0, est1, fff22, fff27, est2;
    br = bsl * vkm * 6 * ae * pow(c70b, 2) / (pi * fkin) * hqet * box;

    if (io == 0)
    {
        br = br * cqed0;
    }

    // CORREZIONE NLO
    if (io == 1)
    {
        cb = -8. * asf(scb) / (3 * pi); 
        // dalla massa polo del b

        double z2 = z * z;
        double z3 = z2 * z;
        double z4 = z3 * z;
        double zs = sqrt(z);
        double zl = log(z);
        // **** ho portato 17./ nella riga sotto
        double hz = -(1 - z2) * (25. / 4. - 239. / 3. * z + 25. / 4. * z2) + z * zl * (20 + 90 * z - 4. / 3. * z2 + 17. / 3. * z3) + z2 * pow(zl, 2) * (36 + z2) + (1 - z2) * (17 - 64 * z + 17 * z2) / 3. * log(1 - z);
        hz = hz - 4 * (1 + 30 * z2 + z4) * zl * log(1 - z) - (1 + 16 * z2 + z4) * (6 * sp(z) - pi2) - 32 * zs * zs * zs * (1 + z) * (pi2 - 4 * sp(zs) + 4 * sp(-zs) - 2 * zl * log((1 - zs) / (1 + zs)));
        csl = 2 * asf(scl) / (3 * pi) * hz / fkin; // dal semileptonico

        double rie3 = 1.202057;
        double r2 = -833 + 144 * pi2 * pow(zs, 3) + (1728 - 180 * pi2 - 1296 * rie3 + (1296 - 324 * pi2) * zl + 108 * pow(zl, 2) + 36 * pow(zl, 3)) * z + (648 + 72 * pi2 + (432 - 216 * pi2) * zl + 36 * pow(zl, 3)) * z2;
        r2 = 2. / 243. * (r2 + (-54 - 84 * pi2 + 1092 * zl - 756 * pow(zl, 2)) * z3);
        double r1 = -r2 / 6.;
        double r7 = -10. / 3. - 8 * pi2 / 9.;
        double r8 = -4. / 27. * (-33 + 2 * pi2);
        double geff1 = -208.0 / 243.0;
        double geff2 = 416.0 / 81.0;
        double geff7 = 32.0 / 3.0;
        double geff8 = -32.0 / 9.0;
        cd = c10b * (r1 + geff1 * log(bm / scb)) + c20b * (r2 + geff2 * log(bm / scb));
        cd = cd + c70b * (r7 + geff7 * log(bm / scb)) + c80b * (r8 + geff8 * log(bm / scb));
        cd = 2.0 * asscb / (4 * pi) * c70b * (c71b + cd); // dall'ampiezza D

        deltg = log(delt);

        if (ityp == 1)
        {
            delt2 = delt * delt;
            delt3 = delt2 * delt;
            deltg1 = log(1 - delt);
            fff77 = (10 * delt + delt2 - 2. / 3. * delt3 + delt * (delt - 4) * deltg) / 3.0;
            //      old formula
            //      fff88=-2*log(rbs)*(delt2+2*delt+2*deltg1)+2*sp(1-delt)-pi2/3.
            //      fff88=(fff88-delt*(1+2*delt)*deltg+8*deltg1-2./3.*delt3+3*delt2
            //     +   +7*delt)/27.
            //      fff78=8./9.*(sp(1-delt)-pi2/6.-delt*deltg+(11*delt-3*delt2+
            //     +   delt3/3.)/4.)
            //      new and correct formula
            fff88 = -2 * log(rbs) * (delt2 + 2 * delt + 4 * deltg1) + 4 * sp(1 - delt) - 2 * pi2 / 3.0;
            fff88 = (fff88 - delt * (2 + delt) * deltg + 8 * deltg1 - 2. / 3. * delt3 + 3 * delt2 + 7 * delt) / 27.0;
            fff78 = 8. / 9. * (sp(1 - delt) - pi2 / 6. - delt * deltg + (9 * delt - delt2 + delt3 / 3.) / 4.);
            eps = 1.0e-4;
            est0 = 0.0;
            est1 = (1 - delt) / z;
            est2 = 1. / z;
            printf("\n\n %e %e %e", est0, est1, eps); // exit(1);
            fff22 = 16 * z / 27. * (delt * gauss1(func1, est0, est1, eps) + gauss1(func2, est1, est2, eps));
            fff27 = -8 * z2 / 9. * (delt * gauss1(func3, est0, est1, eps) + gauss1(func4, est1, est2, eps));
        }

        else if (ityp == 2)
        {
            if (fabs(delt - 0.99) < 0.0001)
            {
                fff77 = 3.42106099217725;
                // old numbers
                //   fff78=0.394065893946580
                //   fff88=0.668299012321822
                // new numbers (correct) for delta=0.99
                fff78 = 0.389665893956536;
                fff88 = 3.216166804631025;
                fff22 = 3.402419197348582e-002;
                fff27 = 1.817503700368764e-002;
            }
            else if (fabs(delt - 0.9) < 0.0001)
            {
                fff77 = 3.20598527956;
                fff78 = 0.38734061186;
                fff88 = 1.31742266135;
                fff22 = 3.40367951120922e-002;
                fff27 = 1.81788563733362e-002;
            }
        }

        double fff28 = -fff27 / 3.;
        double fff11 = fff22 / 36.;
        double fff12 = -fff22 / 3.;
        double fff17 = -fff27 / 6.;
        double fff18 = -fff28 / 6.;
        fnum[1][1] = fff11;
        fnum[1][2] = fff12;
        fnum[1][3] = -0.003493067823644294;
        fnum[1][4] = 0.0005821779706073824;
        fnum[1][5] = -0.04589256102467917;
        fnum[1][6] = -0.05997538991269699;
        fnum[1][7] = fff17;
        fnum[1][8] = fff18;
        fnum[2][2] = fff22;
        fnum[2][3] = 0.02095840694186576;
        fnum[2][4] = -0.003493067823644294;
        fnum[2][5] = 0.2753553661480751;
        fnum[2][6] = 0.359852339476182;
        fnum[2][7] = fff27;
        fnum[2][8] = fff28;
        fnum[3][3] = 0.01399990588537613;
        fnum[3][4] = -0.004666635295125378;
        fnum[3][5] = 0.3276651145509532;
        fnum[3][6] = 0.06660649932577872;
        fnum[3][7] = 0.04208256917382194;
        fnum[3][8] = -0.01402752305794064;
        fnum[4][4] = 0.008791470008620995;
        fnum[4][5] = -0.05461085242515887;
        fnum[4][6] = 0.1569505914595811;
        fnum[4][7] = -0.007013761528970321;
        fnum[4][8] = 0.002337920509656775;
        fnum[5][5] = 1.93694064002221;
        fnum[5][6] = 0.9505876392301583;
        fnum[5][7] = 0.5925919999999998;
        fnum[5][8] = -0.1975306666666666;
        fnum[6][6] = 0.9305162079553912;
        fnum[6][7] = 0.0001919555740097084;
        fnum[6][8] = -0.00006398519133657409;
        fnum[7][7] = fff77;
        fnum[7][8] = fff78;
        fnum[8][8] = fff88;
    }

    double a1 = (exp(-asscb * deltg * (7 + 2 * deltg) / (3 * pi)) - 1) * pow(c70b, 2);
    double a2 = 0.0;

    for (int j = 1; j <= 8; j++)
        for (int i = 1; i <= j; i++)
            a2 = a2 + fnum[i][j] * clead[i] * clead[j];

    a2 = asscb / pi * a2;
    double ca = a1 + a2; // dal bremsstrahlung

    // QED correction after eq 11 of Kagan-Neubert
    double c7em = (32 / 75 * pow(eta, (-9 / 23)) - 40 / 69 * pow(eta, (-7 / 23)) + 88. / 575. * pow(eta, (16 / 23))) * c70 + (-32. / 575. * pow(eta, (-9 / 23)) + 32. / 1449. * pow(eta, (-7 / 23)) + 640. / 1449. * pow(eta, (14.0 / 23.0)) - 704. / 1725. * pow(eta, (16 / 23))) * c80 + -190. / 8073. * pow(eta, (-35 / 23)) - 359. / 3105. * pow(eta, (-17. / 23.)) + 4276. / 121095. * pow(eta, (-12. / 23.)) + 350531. / 1009125. * pow(eta, (-9. / 23.)) + 2. / 4347. * pow(eta, (-7. / 23.)) - 5956 / 15525. * pow(eta, (6. / 23.)) + 38380. / 169533. * pow(eta, (14. / 23.)) - 748. / 8625. * pow(eta, (16. / 23.));

    // here we can decide where to put the QED correction: using cqed or in the
    // overall QED factor qed1 below
    double cqed = 2.0 * ae / asscb * c7em * c70b;
    //       TRUNCATED VERSION
    //       deff=sqrt(abs((1+cb+csl+(cd+ca+cqed)/c70bew**2)*c70bew**2))
    //        NB here the term |D|**2 in eq.16 was expanded in powers of alphas
    //        this may cause problems when c70b is small therefore here we put
    //        also the option without truncation
    //       NON TRUNCATED VERSION

    double deff = sqrt(fabs(pow((c70bew + cd / 2.0 / c70b), 2) + ca + cqed + (cb + csl) * pow(c70b, 2)));

    // - for the same reason we need to recalculate hqc
    hqc = hlam * (c10b / 54.0 - c20b / 9.0) * c70b / pow(deff, 2); // corrected according to NPB version
    double hqetnlo = 1 + hqb / bm2 + hqc / cm2;

    br = br * pow(deff, 2) / pow(c70b, 2) * hqetnlo / hqet;

    // semileptonic qed corr
    double cqed1 = (1 - 2.0 * ae * log(zm / bm) / pi);

    br = br * alpha0 / ae * cqed1;
    *brs = br;
}

//
// calculates the degree of fine-tuning in a given model 
//  (at the moment with respect to mz and mtop only).
// input: mu,tbeta, mhd^2, mhu^2 (at the ewsb scale)
// output:   czmu,czbmu,ctmu,ctbmu are (dimensionless) measures of
// the degree of fine-tuning   on mu and b*mu with respect to mz and mtop,
//	respectively. the larger those numbers (>>1), the more it is "fine-tuned"
//--------------------------------------------------------------------------

void su_finetune(double mu,double tb,double mhd2,double mhu2,double *czmu,double *czbmu,double *ctmu,double *ctbmu)
{
   *czmu = 2*pow(mu,2)/pow(mz,2)*(1.0 + (pow(tb,2)+1.0)/pow((pow(tb,2)-1.0),2)*4*pow(tb,2)*(mhd2-mhu2)/((mhd2-mhu2)*(pow(tb,2)+1.0)- pow(mz,2)*(pow(tb,2)-1.0)));  
   *czbmu = 4*pow(tb,2)*(pow(tb,2)+1.0)/pow((pow(tb,2)-1.0),3)*(mhd2-mhu2)/pow(mz,2);   
   *ctmu = *czmu/2 +2*pow(mu,2)/(mhd2+ mhu2+2*pow(mu,2))/(pow(tb,2)-1.0);        
   *ctbmu = *czbmu/2 +1.0/(1.0-pow(tb,2));
   }
   
