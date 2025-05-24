#include <stdio.h>
#include <math.h>
#include "variables.h"
#include "initcond.h"

// ============================================================================================================
// -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . --
// ============================================================================================================

// ============================================================================================================
//  AMSB model -
// ============================================================================================================

//
// Calculates the initial conditions at initial scale (where the rge starts) in the general AMSB model (i.e. including a soft-susy
// breaking scalar mass m_0 with a different weight c_i for every higgs and sfermion scalar mass).
//
// The input parameters at the initial scale are:
// 	m32: the gravitino mass,
// 	m0 : the soft-susy breaking scalar mass term,
// 	cq,cu,cd,cl,ce,chu,chd: weights  of m0 for the different soft terms,
//
// (for the original amsb model: c_i=0 and usual minimal amsb model: c_i=1),
// 	g12,g22,g23: gauge couplings squared,
// 	ytau,yb,yt : third generation yukawa gauge couplings.
//
// The ouputs at the initial scale are:
// 	m1,m2,m3: gaugino mass terms,
// 	au,ad,al,au1,ad1,al1: 3d and 1st/2d generation trilinear couplings,
// 	mhu2,mhd2,mtaur2,msl2,mbr2,mtr2,msq2,mer2,mel2,mdr2,mur2,muq2,
//
// higgs and sfermion soft mass terms squared.
// -----------------------------------------------------------------------------------------------------------------------------------------

//
//	Phenomenological Consequences of Supersymmetry with Anomaly-Induced Masses
//		Tony Gherghetta, Gian F. Giudice, and James D. Wells
// -----------------------------------------------------------------------------------------------------------------------------------------

//
// Here come the various subroutines for the models
// -----------------------------------------------------------------------------------------------------------------------------------------

void su_amsbsub(double m0, double m32, double cq, double cu, double cd, double cl, double ce, double chu, double chd, double y[], double pass_m[])
{
    double g12 = y[1];
    double g22 = y[2];
    double g32 = y[3];
    double ytau = y[4];
    double yb = y[5];
    double yt = y[6];
    double al = y[9];
    double ad = y[10];
    double au = y[11];
    double al1 = y[29];
    double au1 = y[30];
    double ad1 = y[31];
    double mhu2 = y[12];
    double mhd2 = y[13];
    double mtaur2 = y[14];
    double msl2 = y[15];
    double mbr2 = y[16];
    double mtr2 = y[17];
    double msq2 = y[18];
    double mer2 = y[24];
    double mel2 = y[25];
    double mdr2 = y[26];
    double mur2 = y[27];
    double muq2 = y[28];

    double m1 = pass_m[0];
    double m2 = pass_m[1];
    double m3 = pass_m[2];

    double nf = 6.0;

    double cpi = 1.0 / (16 * pi * pi);
    double ytau2 = pow(ytau, 2);
    double yb2 = pow(yb, 2);
    double ytop2 = pow(yt, 2);

    //  gauginos

    // adding 2-loop gauge+yukawa in boundary conditions:
    // --------------------------------------------------------------------------------------------------------------------------------------

    double betg1 = 33.0 / 5 * g12 + cpi * g12 * ((19 * nf / 15 + 9.0 / 25) * g12 + (3 * nf / 5 + 9.0 / 5) * g22 + 44 * nf / 15 * g32 - 26 * ytop2 / 5 - 14 * yb2 / 5 - 18 * ytau2 / 5); // Eq (57) // Here nf = 2
    double betg2 = g22 + cpi * g22 * ((nf / 5 + 3.0 / 5) * g12 + (7 * nf - 17.0) * g22 + 4 * nf * g32 - 6 * ytop2 - 6 * yb2 - 2 * ytau2);                                               // Eq (58)
    double betg3 = -3 * g32 + cpi * g32 * (11 * nf / 30 * g12 + 3 * nf / 2 * g22 + (34 * nf / 3 - 54.0) * g32 - 4 * ytop2 - 4 * yb2);                                                   // Eq (59)  Check may be there is some mistake

    m1 = cpi * betg1 * m32; // Eq(54)
    m2 = cpi * betg2 * m32; // Eq(55)
    m3 = cpi * betg3 * m32; // Eq(56)

    //
    // soft scalar masses and couplings
    // --------------------------------------------------------------------------------------------------------------------------------------
    double temp = 0.0;

    temp = -22 * ytop2 * ytop2 - 5 * yb2 * yb2 - 5 * yb2 * ytop2 - yb2 * ytau2 + (6 * g12 / 5 + 6 * g22 + 16 * g32) * ytop2 + 2 * g12 / 5 * yb2;
    temp += +(13 * nf / 15 + 403.0 / 450) * pow(g12, 2) + (3 * nf - 21.0 / 2) * pow(g22, 2) + (16 * nf / 3 - 304.0 / 9) * pow(g32, 2);
    temp += +g12 * g22 + 136.0 / 45 * g12 * g32 + 8 * g22 * g32;
    double byt = yt * (-13 * g12 / 15 - 3 * g22 - 16 * g32 / 3 + 6 * ytop2 + yb2) + cpi * yt * temp;

    temp = 0.0;
    temp += -22 * yb2 * yb2 - 5 * ytop2 * ytop2 - 5 * yb2 * ytop2 - 3 * yb2 * ytau2 - 3 * ytau2 * ytau2 + 4 * g12 / 5 * ytop2 + (2 * g12 / 5 + 6 * g22 + 16 * g32) * yb2;
    temp += +6 * g12 / 5 * ytau2 + (7 * nf / 15 + 7.0 / 18) * pow(g12, 2) + (3 * nf - 21.0 / 2) * pow(g22, 2) + (16 * nf / 3 - 304.0 / 9) * pow(g32, 2) + g12 * g22 + 8 * g12 * g32 / 9 + 8 * g22 * g32;
    double byb = yb * (-7 * g12 / 15 - 3 * g22 - 16 * g32 / 3 + ytop2 + 6 * yb2 + ytau2) + cpi * yb * temp;

    temp = 0.0;
    temp += -10 * ytau2 * ytau2 - 9 * yb2 * yb2 - 9 * yb2 * ytau2 - 3 * yb2 * ytop2 + (6 * g22 + 6 * g12 / 5) * ytau2 + (-2 * g12 / 5 + 16 * g32) * yb2;
    temp += +(9 * nf / 5 + 27.0 / 10) * pow(g12, 2) + (3 * nf - 21.0 / 2) * pow(g22, 2) + 9 * g12 * g22 / 5;
    double byta = ytau * (-9 * g12 / 5 - 3 * g22 + 3 * yb2 + 4 * ytau2) + cpi * ytau * temp;

    //
    // trilinear a terms (3rd generation)
    // --------------------------------------------------------------------------------------------------------------------------------------
    au = -byt / yt * cpi * m32;    // Eq(54)
    ad = -byb / yb * cpi * m32;    // Eq(55)
    al = -byta / ytau * cpi * m32; // Eq(56)

    //
    // trilinear a terms (1st and 2d generations)
    // --------------------------------------------------------------------------------------------------------------------------------------

    temp = 0.0;
    temp += -3 * ytau2 * ytau2 - 9 * yb2 * yb2 - 3 * yb2 * ytop2 + 6 * g12 / 5 * ytau2 + (-2 * g12 / 5 + 16 * g32) * yb2;
    temp += +(9 * nf / 5 + 27.0 / 10) * pow(g12, 2) + (3 * nf - 21.0 / 2) * pow(g22, 2) + 9 * g12 * g22 / 5;
    double dyovery4 = ytau2 + 3 * yb2 - 9 * g12 / 5 - 3 * g22 + cpi * temp;

    temp = 0.0;
    temp += -9 * yb2 * yb2 - 3 * yb2 * ytop2 - 3 * ytau2 * ytau2 + (-2 * g12 / 5 + 16 * g32) * yb2 + 6 * g12 / 5 * ytau2;
    temp += +(7 * nf / 15 + 7.0 / 18) * pow(g12, 2) + (3 * nf - 21.e0 / 2) * pow(g22, 2) + (16 * nf / 3 - 304.0 / 9) * pow(g32, 2) + g12 * g22 + 8 * g12 * g32 / 9 + 8 * g22 * g32;
    double dyovery5 = 3 * yb2 + ytau2 - 7 * g12 / 15 - 3 * g22 - 16 * g32 / 3 + cpi * temp;

    temp = 0.0;
    temp += -9 * ytop2 * ytop2 - 3 * yb2 * ytop2 + (4 * g12 / 5 + 16 * g32) * ytop2 + (13 * nf / 15 + 403.0 / 450) * pow(g12, 2);
    temp += +(3 * nf - 21.0 / 2) * pow(g22, 2) + (16 * nf / 3 - 304.0 / 9) * pow(g32, 2) + g12 * g22 + 136.0 / 45 * g12 * g32 + 8 * g22 * g32;
    double dyovery6 = 3 * ytop2 - 13 * g12 / 15 - 3 * g22 - 16 * g32 / 3 + cpi * temp;

    au1 = -dyovery6 * cpi * m32; // Eq (54)
    ad1 = -dyovery5 * cpi * m32; // Eq (55)
    al1 = -dyovery4 * cpi * m32; // Eq (56)

    //
    // 3rd generation fermion masses
    // --------------------------------------------------------------------------------------------------------------------------------------

    msl2 = (-99 * pow(g12, 2) / 50 - 3 * pow(g22, 2) / 2 + ytau * byta) * pow(cpi, 2) * pow(m32, 2) + cl * pow(m0, 2);                           // Eq (52)
    mtaur2 = (-198 * pow(g12, 2) / 25 + 2 * ytau * byta) * pow(cpi, 2) * pow(m32, 2) + ce * pow(m0, 2);                                          // Eq (53)
    msq2 = (-11 * pow(g12, 2) / 50 - 3 * pow(g22, 2) / 2 + 8 * pow(g32, 2) + yt * byt + yb * byb) * pow(cpi, 2) * pow(m32, 2) + cq * pow(m0, 2); // Eq (49)
    mtr2 = (-88 * pow(g12, 2) / 25 + 8 * pow(g32, 2) + 2 * yt * byt) * pow(cpi, 2) * pow(m32, 2) + cu * pow(m0, 2);                              // Eq (47)
    mbr2 = (-22 * pow(g12, 2) / 25 + 8 * pow(g32, 2) + 2 * yb * byb) * pow(cpi, 2) * pow(m32, 2) + cd * pow(m0, 2);                              // Eq (48)

    //
    // 1rst and 2d generations (degenerate with 3rd):
    // --------------------------------------------------------------------------------------------------------------------------------------

    mel2 = (-99 * pow(g12, 2) / 50 - 3 * pow(g22, 2) / 2) * pow(cpi, 2) * pow(m32, 2) + cl * pow(m0, 2);                   // Eq (52)
    mer2 = (-198 * pow(g12, 2) / 25) * pow(cpi, 2) * pow(m32, 2) + ce * pow(m0, 2);                                        // Eq (53)
    muq2 = (-11 * pow(g12, 2) / 50 - 3 * pow(g22, 2) / 2 + 8 * pow(g32, 2)) * pow(cpi, 2) * pow(m32, 2) + cq * pow(m0, 2); // Eq (49)
    mur2 = (-88 * pow(g12, 2) / 25 + 8 * pow(g32, 2)) * pow(cpi, 2) * pow(m32, 2) + cu * pow(m0, 2);                       // Eq (47)
    mdr2 = (-22 * pow(g12, 2) / 25 + 8 * pow(g32, 2)) * pow(cpi, 2) * pow(m32, 2) + cd * pow(m0, 2);                       // Eq (48)

    //
    // higgs mass term^2:
    // --------------------------------------------------------------------------------------------------------------------------------------

    mhu2 = (-99 * pow(g12, 2) / 50 - 3 * pow(g22, 2) / 2 + 3 * yt * byt) * pow(cpi, 2) * pow(m32, 2) + chu * pow(m0, 2);               // Eq(50)
    mhd2 = (-99 * pow(g12, 2) / 50 - 3 * pow(g22, 2) / 2 + 3 * yb * byb + ytau * byta) * pow(cpi, 2) * pow(m32, 2) + chd * pow(m0, 2); // Eq(51)

    y[1] = g12;
    y[2] = g22;
    y[3] = g32;
    y[4] = ytau;
    y[5] = yb;
    y[6] = yt;
    y[9] = al;
    y[10] = ad;
    y[11] = au;
    y[29] = al1;
    y[30] = au1;
    y[31] = ad1;
    y[12] = mhu2;
    y[13] = mhd2;
    y[14] = mtaur2;
    y[15] = msl2;
    y[16] = mbr2;
    y[17] = mtr2;
    y[18] = msq2;
    y[24] = mer2;
    y[25] = mel2;
    y[26] = mdr2;
    y[27] = mur2;
    y[28] = muq2;

    pass_m[0] = m1;
    pass_m[1] = m2;
    pass_m[2] = m3;
}

// ============================================================================================================
// -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . --
// ============================================================================================================

// ============================================================================================================
//  GMSB model -
// ============================================================================================================

//
// Here is defined the Spence function that is needed:
// Ref : http://en.wikipedia.org/wiki/Spence's_function
// -----------------------------------------------------------------------------------------------------------------------------------------

static double su_pli2(double x)
{

    double rz, sp;
    int ii;
    double z;

    rz = 0.0;
    if (x <= 0.50)
    {
        z = -log(1.0 - x);
        for (int i = 1; i <= 9; i++)
        {
            ii = 2 * i + 1;
            rz = rz + b[i] / f[ii] * pow(z, ii);
        }
        sp = b0 * z + b1 * z * z / 2.0 + rz;
    }
    else
    {
        z = -log(x);
        for (int i = 1; i <= 9; i++)
        {
            ii = 2 * i + 1;
            rz = rz + b[i] / f[ii] * pow(z, ii);
        }
        sp = -(b0 * z + b1 * z * z / 2.0 + rz) + pi * pi / 6.0 - log(x) * log(1.0 - x);
    }
    return (sp);
}

//
//  Calculates the  gmsb model initial conditions at the messenger scale  (where the rge starts).
//
//	 The input at the messenger scale aremgmmess:
//  	mgmmess,mgmsusy: messenger and susy-breaking scales,
//  	nl, nq number of lepton/ quark messengers (in minimal gmsb, nl=nq=1),
//  	g12,g22,g23: gauge couplings squared.
//
//  The output parameters at the messenger scale are:
//  	m1,m2,m3, the gaugino masses,
//  	au,ad,al,au1,ad1,al1: the trilinear sfermion couplings,
//  	mhu2,mhd2,mtaur2,msl2,mbr2,mtr2,msq2,mer2,mel2,mdr2,mur2,muq2:
//  	higgs and sfermion soft mass terms squared.
//
//  The routine needs to evaluate a spence function which is supplied:
// -----------------------------------------------------------------------------------------------------------------------------------------

//
// SuSpect: a Fortran Code for the Supersymmetric and Higgs Particle Spectrum in the MSSM
//		Abdelhak DJOUADI, Jean KNEUR and Gilbert MOULTAKA
//	arXiv:hep-ph/0211331v2
// -----------------------------------------------------------------------------------------------------------------------------------------

void su_gmsbsub(double mgmmess, double mgmsusy, double nl, double nq, double y[], double pass_m[])
{

    double g12 = y[1];
    double g22 = y[2];
    double g32 = y[3];
    double ytau = y[4];
    double yb = y[5];
    double yt = y[6];
    double al = y[9];
    double ad = y[10];
    double au = y[11];
    double al1 = y[29];
    double au1 = y[30];
    double ad1 = y[31];
    double mhu2 = y[12];
    double mhd2 = y[13];
    double mtaur2 = y[14];
    double msl2 = y[15];
    double mbr2 = y[16];
    double mtr2 = y[17];
    double msq2 = y[18];
    double mer2 = y[24];
    double mel2 = y[25];
    double mdr2 = y[26];
    double mur2 = y[27];
    double muq2 = y[28];

    double m1 = pass_m[0];
    double m2 = pass_m[1];
    double m3 = pass_m[2];

    //===== defining one-loop functions and related material:
    //

    //
    // Computation of factorial for spence fction:
    //	Ref : http://en.wikipedia.org/wiki/Spence's_function
    // --------------------------------------------------------------------------------------------------------------------------------------

    f[1] = 1.0;
    for (int k = 2; k <= 19; k++)
        f[k] = f[k - 1] * k;

    //
    // Computation of the first relevant bernouilli numbers:
    // Ref : http://en.wikipedia.org/wiki/Bernoulli_number
    // --------------------------------------------------------------------------------------------------------------------------------------

    b0 = 1.0;
    b1 = -1.0 / 2.0;
    b[1] = 1.0 / 6.0;
    b[2] = -1.0 / 30.0;
    b[3] = 1.0 / 42.0;
    b[4] = -1.0 / 30.0;
    b[5] = 5.0 / 66.0;
    b[6] = -691.0 / 2730.0;
    b[7] = 7.0 / 6.0;
    b[8] = -3617.0 / 510.0;
    b[9] = 43867.0 / 798.0;

    //
    // Def: mgmsusy= f/s; mgmmess=lambdas, x=f/(lambdas^2)=mgmsusy/mgmmess
    // 	where lambda is the coupling of messengers to the goldstino in the superpotential  gaugino masses:
    //	Ref : Eq(22)
    // --------------------------------------------------------------------------------------------------------------------------------------

    double f1_x, f2_x, x;
    x = mgmsusy / mgmmess;
    f1_x = ((1.0 + x) * log(1.0 + x) + (1.0 - x) * log(1.0 - x)) / pow(x, 2);
    f2_x = (1.0 + x) / pow(x, 2) * (log(1.0 + x) - 2 * su_pli2(x / (1.0 + x)) + 1.0 / 2 * su_pli2(2 * x / (1.0 + x))) + (1.0 - x) / pow(x, 2) * (log(1.0 - x) - 2 * su_pli2(-x / (1.0 - x)) + 1.0 / 2 * su_pli2(-2 * x / (1.0 - x)));

    //
    // nq "quark" messengers, nl "lepton" messengers: under su[3]*su[2]*u[1]
    // we choose:
    //   q ~ ( 3,   1,   -1/3)
    //   l ~  (1,   2,    1/2) *)
    //
    // dc=sum of dynkin(messengers) times casimir(scalar partners) convoluted by
    // the corresponding gauge couplings
    // minimal gmsb: dc for the various mssm scalar partners *)
    // --------------------------------------------------------------------------------------------------------------------------------------

    //
    // Storing th Dynking indices
    // Ref : Eq(24), Eq(25), Eq(26)
    // --------------------------------------------------------------------------------------------------------------------------------------

    double yq = 1.0 / 6;
    double yu = -2.0 / 3;
    double yd = 1.0 / 3;
    double yl = -1.0 / 2;
    double yr = 1.0;
    double yh1 = -1.0 / 2;
    double yh2 = 1.0 / 2;

    al1 = g12 / (16 * pi * pi);
    double al2 = g22 / (16 * pi * pi);
    double al3 = g32 / (16 * pi * pi);
    double dcq = 4.0 / 3 * nq * al3 * al3 + 3.0 / 4 * nl * al2 * al2 + 3.0 / 5 * yq * yq / 5 * (2 * nq + 3 * nl) * pow(al1, 2);
    double dcu = 4.0 / 3 * nq * al3 * al3 + 3.0 / 5 * yu * yu / 5 * (2 * nq + 3 * nl) * al1 * al1;
    double dcd = 4.0 / 3 * nq * al3 * al3 + 3.0 / 5 * yd * yd / 5 * (2 * nq + 3 * nl) * al1 * al1;
    double dcl = 3.0 / 4 * nl * al2 * al3 + 3.0 / 5 * yl * yl / 5 * (2 * nq + 3 * nl) * al1 * al1;
    double dce = 3.0 / 5 * yr * yr / 5 * (2 * nq + 3 * nl) * al1 * al1;
    double dchd = 3.0 / 4 * nl * al2 * al2 + 3.0 / 5 * pow(yh1, 2) / 5 * (2 * nq + 3 * nl) * al1 * al1;
    double dchu = 3.0 / 4 * nl * al2 * al2 + 3.0 / 5 * pow(yh2, 2) / 5 * (2 * nq + 3 * nl) * al1 * al1;

    //
    // dc_i are group factor combinations for soft terms mi^2
    // gauginos
    // --------------------------------------------------------------------------------------------------------------------------------------

    double f1x;

    if (x == 1.0)
        f1x = 1.38629;
    else
        f1x = f1_x;

    m1 = g12 / (16 * pi * pi) * mgmsusy / 5 * (2 * nq + 3 * nl) * f1x;
    m2 = g22 / (16 * pi * pi) * mgmsusy * nl * f1x;
    m3 = g32 / (16 * pi * pi) * mgmsusy * nq * f1x;

    //
    // Soft scalar masses and couplings
    // --------------------------------------------------------------------------------------------------------------------------------------

    au = 0.0;
    ad = 0.0;
    al = 0.0;
    au1 = 0.0;
    ad1 = 0.0;
    al1 = 0.0;

    //
    // 3rd generation
    // --------------------------------------------------------------------------------------------------------------------------------------

    double f2x;

    if (x == 1.0)
        f2x = 0.702266;
    else
        f2x = f2_x;

    msl2 = 2 * pow(mgmsusy, 2) * dcl * f2x;
    mtaur2 = 2 * pow(mgmsusy, 2) * dce * f2x;
    msq2 = 2 * pow(mgmsusy, 2) * dcq * f2x;
    mtr2 = 2 * pow(mgmsusy, 2) * dcu * f2x;
    mbr2 = 2 * pow(mgmsusy, 2) * dcd * f2x;

    //
    // 1rst and 2d generations (degenerate with 3rd):
    // --------------------------------------------------------------------------------------------------------------------------------------

    mel2 = 2 * pow(mgmsusy, 2) * dcl * f2x;
    mer2 = 2 * pow(mgmsusy, 2) * dce * f2x;
    muq2 = 2 * pow(mgmsusy, 2) * dcq * f2x;
    mur2 = 2 * pow(mgmsusy, 2) * dcu * f2x;
    mdr2 = 2 * pow(mgmsusy, 2) * dcd * f2x;

    //
    // higgs mass term^2:
    // --------------------------------------------------------------------------------------------------------------------------------------

    mhd2 = 2 * pow(mgmsusy, 2) * dchd * f2x;
    mhu2 = 2 * pow(mgmsusy, 2) * dchu * f2x;

    y[1] = g12;
    y[2] = g22;
    y[3] = g32;
    y[4] = ytau;
    y[5] = yb;
    y[6] = yt;
    y[9] = al;
    y[10] = ad;
    y[11] = au;
    y[29] = al1;
    y[30] = au1;
    y[31] = ad1;
    y[12] = mhu2;
    y[13] = mhd2;
    y[14] = mtaur2;
    y[15] = msl2;
    y[16] = mbr2;
    y[17] = mtr2;
    y[18] = msq2;
    y[24] = mer2;
    y[25] = mel2;
    y[26] = mdr2;
    y[27] = mur2;
    y[28] = muq2;

    pass_m[0] = m1;
    pass_m[1] = m2;
    pass_m[2] = m3;
}
