#include <math.h>
#include "variables.h"
#include "radiative.h"

//
// Following are the main routine for the RGE evolution of parameter between low and high energy scales. It returns a set of n mass and coupling
// parameters "y" at a specified scale exp(x2) when given at an initial scale exp(x2). its uses the beta functions in the subroutine su_derivs
// and solves the coupled differential equations with the su_rkqc() runge-kutta subroutine. thus y(n) is a  vector containing all the n rg
// evolving parameters at  various possible scales depending on evolution stages.
// -------------------------------------------------------------------------------------------------------------------------------------------------

//
//	The parameters are :
//
//  y[1] = g1^2   		u[1] gauge coupling squared
//  y[2] = g2^2   		su[2]_l gauge coupling squared
//  y[3] = g3^2   		su[3] gauge coupling squared
//  y[4] = y_tau  		tau lepton yukawa coupling
//  y[5] = y_b    		bottom  quark yukawa coupling
//  y[6] = y_top  		top quark yukawa coupling
//  y[7] = ln(vu) 		logarithm of the vev vu
//  y[8] = ln(vd) 		logarithm of the vev vd
//  y[9] = a_tau  		trilinear coupling for stau
//  y[10]= a_b    		trilinear coupling for sbottom
//  y[11]= a_top  		trilinear coupling for stop
//  y[12]= (m_phi_u)^2  scalar phi_u mass term squared
//  y[13]= (m_phi_d)^2  scalar phi_d mass term squared
//  y[14]= mtaur^2 		right-handed stau mass term squared
//  y[15]= msl^2   		left-handed stau mass term squared
//  y[16]= mbr^2   		right-handed sbottom mass term squared
//  y[17]= mtr^2   		right-handed stop mass term squared
//  y[18]= msq^2   		left-handed stop mass term squared
//  y[19]= b       		the (dimensionful) bilinear parameter b
//  y[20]= ln(|m1|) 		logarithm of the bino mass term
//  y[21]= ln(|m2|) 		logarithm of the wino mass term
//  y[22]= ln(|m3|) 		logarithm of the gluino mass term
//  y[23]= ln(|mu|) 		logarithm of the |mu| parameter
//  y[24]= mer^2 			right-handed selectron (smuon) mass term squared
//  y[25]= mel^2 			left-handed selectron (smuon) mass term squared
//  y[26]= mdr^2 			right-handed sdown (sstrange) mass term squared
//  y[27]= mur^2 			right-handed sup (scharm) mass term squared
//  y[28]= muq^2 			left-handed sup (scharm) mass term squared
//  y[29]= a_l   			trilinear coupling for selectron (smuon)
//  y[30]= a_d   			trilinear coupling for sdown (sstrange)
//  y[31]= a_u   			trilinear coupling for sup (scharm).
// -------------------------------------------------------------------------------------------------------------------------------------------------

//
// Note that the number of running parameters consist of :
//		the 22 parameters of the phenomenological mssm;
//		the 3 gauge and the  3 yukawa couplings;
//  	the 3 parameters (vu, vd, b) which are in fact linearly dependent of others.
// -------------------------------------------------------------------------------------------------------------------------------------------------

//
// su_deriv1(x,y,dydx)		  											// All the RGEs are at 1-loop order (faster)
// -------------------------------------------------------------------------------------------------------------------------------------------------

// These are the derivatives of the rg running parameters y(xn), i.e the beta functions beta(y)=d(y)/dln(q). the analytic expressions of the functions
// are taken from (up to some sign conventions which have been changed):
// deriv1 : includes only 1-loop rge with full mssm threshold.
// -------------------------------------------------------------------------------------------------------------------------------------------------

//
// 1. Renormalization group study of the standard model and its extensions: The minimal supersymmetric standard model
// 	By : castano, ramond, piard
// 	Phys. Rev. D49 (1994) 4882,

//
// 2. Supersymmetric particles spectrum
// 	V. Barger, M. S. Berger, P. Ohmann
// 	Phys. Rev. D49 (1994) 4908.
// ---------------------------------------------------------------------------------------------------------------

void su_deriv1(double x, double y[], double dydx[])
{

    static double nn = 3.0;
    static double nd = 3.0;
    static double ne = 3.0;

    // pi=4*atan(1.0)
    double q = exp(x);

    // g1, g2 gauge unif.:
    if (kunif != 0 && h1 > 0.0 && q > 1.0e15)
    {
        if (y[1] >= y[2])
        {
            if (ifirst == 0)
            {
                egut = q;
                // c freeze out the gauge +yukawa+vu,vd couplings after that

                for (int ii = 1; ii <= 31; ii++)
                    ygut[ii] = y[ii];
            }
            ifirst = 1;
        }
    }

    // Simple (unique scale) threshold in beta functions:
    double st1 = 1.0; // (full mssm rges)
    double st2 = 1.0;
    double nu = 2.0 + st1;
    double nsq = 3 * st2;
    double nsu = 3 * st2;
    double nsd = 3 * st2;
    double nsl = 3 * st2;
    double nse = 3 * st2;
    double nhino = 2.0 * st2;
    double nwino = st2;
    double ngino = st2;
    double nh = 1.0 + st2;

    //
    // Coefficient of the beta functions for gauge couplings
    // Eq(B6), Eq(B7a), Eq(B7b)
    // ---------------------------------------------------------------------------------------------------------------

    double b10 = 2.0 / 5 * (17.0 * nu / 12 + 5.0 * nd / 12 + 5.0 * ne / 4.0 + nn / 4.0);
    b10 += +nsq / 30.0;
    b10 += +4.0 * nsu / 15.0 + nsd / 15.0 + nsl / 10.0 + nse / 5.0;
    b10 += +nhino / 5.0 + nh / 10.0;

    double b20 = -22.0 / 3.0 + (nu + nd) / 2.0 + (ne + nn) / 6.0;
    b20 += +nsq / 2.0 + nsl / 6.0 + nhino / 3.0;
    b20 += +nh / 6.0 + 4 * nwino / 3.0;

    double b30 = 2.0 * (nu + nd) / 3.0 + nsq / 3.0 + nsu / 6.0 + nsd / 6.0 + 2 * ngino - 11.0;

    //
    // Gauge coupling beta functions (nb the variables are g^2):
    // Variables :
    // 	-y[1]--y[3]: g1^2,g2^2,g3^2.
    // 	cpi = 1/(4*pi)^2
    // Ref :  Eq(B1)     .. (ps. here we are taking g^2 as variable not g)
    // ---------------------------------------------------------------------------------------------------------------

    dydx[1] = 2 * cpi * b10 * y[1] * y[1];
    dydx[2] = 2 * cpi * b20 * y[2] * y[2];
    dydx[3] = 2 * cpi * b30 * y[3] * y[3];

    //
    // Yukawa coupling beta function (only ytau,yb,ytop included):
    // Variables :
    // 	-y[4]--y[6]: ytau,yb,ytop
    // Ref :  Eq(B16a), Eq(B16b), Eq(B16c)   .. (Here the equations are for tau, bottom, and top. which
    //															will be similar to the electron, up and down respectively )
    // ---------------------------------------------------------------------------------------------------------------

    double ytau2 = y[4] * y[4];
    double yb2 = y[5] * y[5];
    double ytop2 = y[6] * y[6];

    double ytaubeta = 3 * y[4] * ytau2;
    ytaubeta += +(3 * yb2 + ytau2) * y[4];
    ytaubeta += +(-3.0 / 5 * y[1] * (15.0 / 4 + 3.0 / 4 - (1.0 / 4 + 1.0 + 1.0 / 4)) - y[2] * (9.0 / 4 + 9.0 / 4 - 3.0 * 2.0 / 4)) * y[4];

    double ybbeta = 3 * y[5] * yb2;
    ybbeta += y[5] * ytop2 + (3 * yb2 + ytau2) * y[5];
    ybbeta += +(-3.0 / 5 * y[1] * (5.0 / 12 + 3.0 / 4 - (1.0 / 36 + 1.0 / 9 + 1.0 / 4)) - y[2] * (9.0 / 4 + 9.0 / 4 - 3 * 2.0 / 4) - y[3] * (8.0 - 4.0 * 2 / 3)) * y[5];

    double ytbeta = 3 * y[6] * ytop2;
    ytbeta += +y[6] * yb2 + 3 * y[6] * ytop2;
    ytbeta += +(-3.0 / 5 * y[1] * (17.0 / 12 + 3.0 / 4 - (1.0 / 36 + 4.0 / 9 + 1.0 / 4)) - y[2] * (9.0 / 4 + 9.0 / 4 - 3.0 * 2 / 4.0) - y[3] * (8.0 - 4.0 * 2.0 / 3.0)) * y[6];

    /////// I guess the upper part is not required because we are only taking the one loop contributions. So may be we can cut this part later..
    /////// Check what happens in the rest part of this function............

    //
    // Yukawa coupling beta function (only ytau,yb,ytop included):
    // Variables :
    // 	-y[4]--y[6]: ytau,yb,ytop
    // Ref :  Eq(A5), Eq(A9), Eq(A10),Eq(A11)   .. (Here the equations are for tau, bottom, and top. which
    //																will be similar to the electron, up and down respectively )
    //
    // p.s. Here we are only taking the one loop contributions.
    // ---------------------------------------------------------------------------------------------------------------

    dydx[4] = cpi * y[4] * (4 * ytau2 + 3 * yb2 - 9 * y[1] / 5.0 - 3 * y[2]);
    dydx[5] = cpi * y[5] * (6 * yb2 + ytau2 + ytop2 - 7 * y[1] / 15. - 3.0 * y[2] - 16.0 * y[3] / 3.0);
    dydx[6] = cpi * y[6] * (6 * ytop2 + yb2 - 13 * y[1] / 15.0 - 3 * y[2] - 16.0 * y[3] / 3.0);

    //
    // Higgs vev beta functions:
    // 	- y[7], y[8] = ln(vu), ln(vd)
    // Ref : Eq(A6), Eq(A12), Eq(A13)
    // ---------------------------------------------------------------------------------------------------------------

    dydx[7] = cpi * (3.0 / 4 * (y[1] / 5.0 + y[2]) - 3.0 * ytop2);
    dydx[8] = cpi * (3.0 / 4 * (y[1] / 5.0 + y[2]) - 3.0 * yb2 - ytau2);

    //
    // Soft susy-breaking terms beta functions:
    //
    // 	- y[9]--y[11] : atau, ab, atop
    // Ref : Eq(A6), Eq(A12), Eq(A13) [Note that the matrices are here only number as we are not considering the first two flavors and only the third flavour is used]
    // ---------------------------------------------------------------------------------------------------------------

    dydx[9] = cpi * (8.0 * ytau2 * y[9] + 6.0 * yb2 * y[10] + 6.0 * (3.0 * y[1] * sgnm1 * exp(y[20]) / 5.0 + sgnm2 * y[2] * exp(y[21])));
    dydx[10] = cpi * (12.0 * y[10] * yb2 + 2.0 * y[9] * ytau2 + 2.0 * y[11] * ytop2 + 14.0 * y[1] / 15.0 * sgnm1 * exp(y[20]) + 6.0 * y[2] * sgnm2 * exp(y[21]) + 32 * y[3] / 3 * sgnm3 * exp(y[22]));
    dydx[11] = cpi * (12.0 * y[11] * ytop2 + 2.0 * y[10] * yb2 + 26.0 * y[1] / 15.0 * sgnm1 * exp(y[20]) + 6.0 * y[2] * sgnm2 * exp(y[21]) + 32 * y[3] / 3 * sgnm3 * exp(y[22]));

    //
    // Mass terms :
    // 	- y[12]--y[13] : m^2(phi_u), m^2(phi_d)
    // Ref : Eq(A17), Eq(A18)
    // ---------------------------------------------------------------------------------------------------------------

    double trym2 = y[18] - 2 * y[17] + y[16] - y[15] + y[14] + y[12] - y[13] + 2 * (y[28] - 2 * y[27] + y[26] - y[25] + y[24]);
    dydx[12] = 2 * cpi * (3 * ytop2 * (y[12] + y[18] + y[17] + y[11] * y[11]) + 3.0 / 10 * y[1] * trym2 - 3 * y[1] / 5 * exp(2 * y[20]) - 3 * y[2] * exp(2 * y[21]));
    dydx[13] = 2 * cpi * (ytau2 * (y[13] + y[15] + y[14] + y[9] * y[9]) + 3 * yb2 * (y[13] + y[18] + y[16] + y[10] * y[10]) - 3.0 / 10 * y[1] * trym2 - 3 * y[1] / 5 * exp(2 * y[20]) - 3 * y[2] * exp(2 * y[21]));

    //
    // - (1-loop) y[14]--y[19] : m^2_tau, m^2_l, m^2_b, m^2_top, m^2_q, b
    // Ref : Eq(A19), Eq(A20), Eq(A21), Eq(A22), Eq(A23), Eq(A24)
    // ---------------------------------------------------------------------------------------------------------------

    dydx[14] = 2 * cpi * (2 * ytau2 * (y[13] + y[14] + y[15] + y[9] * y[9]) + 3 * y[1] / 5 * trym2 - 12 * y[1] / 5 * exp(2 * y[20]));
    dydx[15] = 2 * cpi * (ytau2 * (y[13] + y[15] + y[14] + y[9] * y[9]) - 3 * y[1] / 10 * trym2 - 3 * y[1] / 5 * exp(2 * y[20]) - 3 * y[2] * exp(2 * y[21]));
    dydx[16] = 2 * cpi * (2 * yb2 * (y[13] + y[16] + y[18] + y[10] * y[10]) + y[1] / 5 * trym2 - 4 * y[1] / 15 * exp(2 * y[20]) - 16 * y[3] / 3 * exp(2 * y[22]));
    dydx[17] = 2 * cpi * (2 * ytop2 * (y[12] + y[17] + y[18] + y[11] * y[11]) - 2 * y[1] / 5 * trym2 - 16 * y[1] / 15 * exp(2 * y[20]) - 16 * y[3] / 3 * exp(2 * y[22]));
    dydx[18] = 2 * cpi * (ytop2 * (y[12] + y[17] + y[18] + y[11] * y[11]) + yb2 * (y[13] + y[18] + y[16] + y[10] * y[10]) + y[1] / 10 * trym2 - y[1] / 15 * exp(2 * y[20]) - 3 * y[2] * exp(2 * y[21]) - 16 * y[3] / 3 * exp(2 * y[22]));
    dydx[19] = 2 * cpi * (3 * y[11] * ytop2 + 3 * y[10] * yb2 + y[9] * ytau2 + 3 * y[1] / 5 * sgnm1 * exp(y[20]) + 3 * y[2] * sgnm2 * exp(y[21]));

    //
    // Gauginos masses beta functions:
    // 	- y[20]--y[22] : ln (m1,m2,m3)
    // Ref : Eq(A26), Eq(A7)
    // ---------------------------------------------------------------------------------------------------------------

    dydx[20] = -2 * cpi * (-3.0 / 5 - nf) * y[1];
    dydx[21] = -2 * cpi * (5.0 - nf) * y[2];
    dydx[22] = -2 * cpi * (9.0 - nf) * y[3];

    //
    // The mu parameter:
    // 	- y(23) = ln mu
    // Ref : Eq(A8)
    // ---------------------------------------------------------------------------------------------------------------

    dydx[23] = cpi * (3 * ytop2 + 3 * yb2 + ytau2 - 3 * y[1] / 5 - 3 * y[2]);

    //
    // 1st and 2d gen. sfermion mass^2 terms:
    // 	- y(24)--y(28) : m^2_er, m^2_el, m^2_dr, m^2_ur, m^2_ul
    // Ref : Eq(A19), Eq(A20), Eq(A21), Eq(A22), Eq(A23)
    // ---------------------------------------------------------------------------------------------------------------

    dydx[24] = 2 * cpi * (3 * y[1] / 5 * trym2 - 12 * y[1] / 5 * exp(2 * y[20]));
    dydx[25] = 2 * cpi * (-3 * y[1] / 10 * trym2 - 3 * y[1] / 5 * exp(2 * y[20]) - 3 * y[2] * exp(2 * y[21]));
    dydx[26] = 2 * cpi * (y[1] / 5 * trym2 - 4 * y[1] / 15 * exp(2 * y[20]) - 16 * y[3] / 3 * exp(2 * y[22]));
    dydx[27] = 2 * cpi * (-2 * y[1] / 5 * trym2 - 16 * y[1] / 15 * exp(2 * y[20]) - 16 * y[3] / 3 * exp(2 * y[22]));
    dydx[28] = 2 * cpi * (y[1] / 10 * trym2 - y[1] / 15 * exp(2 * y[20]) - 3 * y[2] * exp(2 * y[21]) - 16 * y[3] / 3 * exp(2 * y[22]));

    //
    // Trilinear coupling for sstrange scharm and smuon
    // 	- y(29)--y(31) : ae (anu), ad (as), au (ac)
    // Ref : Eq(A14), Eq(A15), Eq(A16)
    // ---------------------------------------------------------------------------------------------------------------

    dydx[29] = cpi * (2 * ytau2 * y[9] + 6 * yb2 * y[10] + 6 * (3 * y[1] * sgnm1 * exp(y[20]) / 5 + sgnm2 * y[2] * exp(y[21])));
    dydx[30] = cpi * (2 * ytau2 * y[9] + 6 * yb2 * y[10] + 14 * y[1] / 15 * sgnm1 * exp(y[20]) + 6 * y[2] * sgnm2 * exp(y[21]) + 32 * y[3] / 3 * sgnm3 * exp(y[22]));
    dydx[31] = cpi * (6 * ytop2 * y[11] + 26 * y[1] / 15 * sgnm1 * exp(y[20]) + 6 * y[2] * sgnm2 * exp(y[21]) + 32 * y[3] / 3 * sgnm3 * exp(y[22]));
}

//
// Variables :
//  y[1] = g1^2   		u[1] gauge coupling squared
//  y[2] = g2^2   		su[2]_l gauge coupling squared
//  y[3] = g3^2   		su[3] gauge coupling squared
//  y[4] = y_tau  		tau lepton yukawa coupling
//  y[5] = y_b    		bottom  quark yukawa coupling
//  y[6] = y_top  		top quark yukawa coupling
//  y[7] = ln(vu) 		logarithm of the vev vu
//  y[8] = ln(vd) 		logarithm of the vev vd
//  y[9] = a_tau  		trilinear coupling for stau
//  y[10]= a_b    		trilinear coupling for sbottom
//  y[11]= a_top  		trilinear coupling for stop
//  y[12]= (m_phi_u)^2  scalar phi_u mass term squared
//  y[13]= (m_phi_d)^2  scalar phi_d mass term squared
//  y[14]= mtaur^2 		right-handed stau mass term squared
//  y[15]= msl^2   		left-handed stau mass term squared
//  y[16]= mbr^2   		right-handed sbottom mass term squared
//  y[17]= mtr^2   		right-handed stop mass term squared
//  y[18]= msq^2   		left-handed stop mass term squared
//  y[19]= b       		the (dimensionful) bilinear parameter b
//  y[20]= ln(|m1|) 		logarithm of the bino mass term
//  y[21]= ln(|m2|) 		logarithm of the wino mass term
//  y[22]= ln(|m3|) 		logarithm of the gluino mass term
//  y[23]= ln(|mu|) 		logarithm of the |mu| parameter
//  y[24]= mer^2 			right-handed selectron (smuon) mass term squared
//  y[25]= mel^2 			left-handed selectron (smuon) mass term squared
//  y[26]= mdr^2 			right-handed sdown (sstrange) mass term squared
//  y[27]= mur^2 			right-handed sup (scharm) mass term squared
//  y[28]= muq^2 			left-handed sup (scharm) mass term squared
//  y[29]= a_l   			trilinear coupling for selectron (smuon)
//  y[30]= a_d   			trilinear coupling for sdown (sstrange)
//  y[31]= a_u   			trilinear coupling for sup (scharm).
// -------------------------------------------------------------------------------------------------------------------------------------------------

//
// 1. Two-Loop Renormalization Group Equations for Soft Supersymmetry-Breaking Couplings
// 	Stephen P. Martin and Michael T. Vaughn
//    hep-ph/9311340, phys.rev. d50 (1994) 2282

//
// 2. Renormalization group study of the standard model and its extensions: The minimal supersymmetric standard model
// 	By : castano, ramond, piard
// 	Phys. Rev. D49 (1994) 4882,
// ---------------------------------------------------------------------------------------------------------------

//
// su_deriv2(x,y,dydx)													// 2-loop RGEs for gauge+Yukawas+gauginos
// deriv2 : includes 2-loop rge for gauge, yukawa cpls, gaugino masses, m_hu,d
// ---------------------------------------------------------------------------------------------------------------

void su_deriv2(double x, double y[], double dydx[])
{

    static double nn = 3.0;
    static double nd = 3.0;
    static double ne = 3.0;

    double q = exp(x);
    // g1, g2 gauge unif.:

    if (kunif != 0 && h1 > 0.0 && q > 1.0e15)
    {
        if (y[1] >= y[2])
        {
            if (ifirst == 0)
            {
                egut = q; // freeze out the gauge +yukawa+vu,vd couplings at gut scale:

                for (int ii = 1; ii <= 31; ii++)
                    ygut[ii] = y[ii];
            }
            ifirst = 1;
        }
    }

    double st1 = 1.0; //	(full mssm rges)
    double st2 = 1.0;
    double nu = 2.0 + st1;
    double nsq = 3 * st2;
    double nsu = 3 * st2;
    double nsd = 3 * st2;
    double nsl = 3 * st2;
    double nse = 3 * st2;
    double nhino = 2.0 * st2;
    double nwino = st2;
    double ngino = st2;
    double nh = 1.0 + st2;

    //
    // 2-loop gauge coupling beta functions (nb the variables are g^2):
    // Variables :
    //		-y[1]--y[3]: g1^2,g2^2,g3^2.
    //		nf = number of flavours = 6
    //		cpi = 1/16pi
    //	Ref : Eq(4.4), Eq(4.5)          .....   Remember that here Y[1] - g1^2 .. whereas in the paper the equations are for g1 etc
    // ------------------------------------------------------------------------------------------------------------------------------------------------

    double b10 = 2.0 / 5 * (17 * nu / 12 + 5 * nd / 12 + 5 * ne / 4 + nn / 4) + nsq / 30 + 4 * nsu / 15 + nsd / 15 + nsl / 10 + nse / 5 + nhino / 5 + nh / 10;
    double b20 = -22.0 / 3 + (nu + nd) / 2 + (ne + nn) / 6 + nsq / 2 + nsl / 6 + nhino / 3 + nh / 6 + 4 * nwino / 3;
    double b30 = 2 * (nu + nd) / 3 + nsq / 3 + nsu / 6 + nsd / 6 + 2 * ngino - 11.0;

    double ytau2 = y[4] * y[4];
    double yb2 = y[5] * y[5];
    double ytop2 = y[6] * y[6];
    double mm1 = sgnm1 * exp(y[20]);
    double mm2 = sgnm2 * exp(y[21]);
    double mm3 = sgnm3 * exp(y[22]);


    dydx[1] = 2.0 * cpi * b10 * y[1] * y[1] + 2.0 * cpi * cpi * y[1] * y[1] * ((19.0 * nf / 15 + 9.0 / 25) * y[1] + (3.0 * nf / 5 + 9.0 / 5) * y[2] + 44.0 * nf / 15 * y[3] - 26.0 * ytop2 / 5.0 - 14.0 * yb2 / 5.0 - 18.0 * ytau2 / 5);
    dydx[2] = 2.0 * cpi * b20 * y[2] * y[2] + 2.0 * cpi * cpi * y[2] * y[2] * ((nf / 5.0 + 3.0 / 5) * y[1] + (7.0 * nf - 17.0) * y[2] + 4.0 * nf * y[3] - 6.0 * ytop2 - 6.0 * yb2 - 2.0 * ytau2);
    dydx[3] = 2.0 * cpi * b30 * y[3] * y[3] + 2.0 * cpi * cpi * y[3] * y[3] * (11.0 * nf / 30.0 * y[1] + 3.0 * nf / 2.0 * y[2] + (34.0 * nf / 3 - 54.0) * y[3] - 4.0 * ytop2 - 4.0 * yb2);

    //
    // 2-loop yukawa coupling beta function (only ytau,yb,ytop included):
    // 	-y[4]--y[6]: ytau,yb,ytop
    // Ref : Eq(4.8), Eq(4.11) - Eq(4.16)  - [ Also Eq(A27) of https://arxiv.org/pdf/hep-ph/9308335]
    // ------------------------------------------------------------------------------------------------------------------------------------------------

    double ytaubeta = 0.0;
    ytaubeta += 3 * y[4] * ytau2;
    ytaubeta += +(3 * yb2 + ytau2) * y[4];
    ytaubeta += +(-3.0 / 5 * y[1] * (15.0 / 4 + 3.0 / 4 - (1.0 / 4 + 1.0 + 1.0 / 4)) - y[2] * (9.0 / 4 + 9.0 / 4 - 3.0 * 2 / 4)) * y[4];

    double ybbeta = 0.0;
    ybbeta += 3 * y[5] * yb2;
    ybbeta += y[5] * ytop2 + (3 * yb2 + ytau2) * y[5];
    ybbeta += +(-3.0 / 5 * y[1] * (5.0 / 12 + 3.0 / 4 - (1.0 / 36 + 1.0 / 9 + 1.0 / 4)) - y[2] * (9.0 / 4 + 9.0 / 4 - 3 * (2.0) / 4) - y[3] * (8.0 - 4.0 * 2 / 3)) * y[5];

    double ytbeta = 0.0;
    ytbeta += 3 * y[6] * ytop2;
    ytbeta += y[6] * yb2 + 3 * y[6] * ytop2;
    ytbeta += +(-3.0 / 5 * y[1] * (17.0 / 12 + 3.0 / 4 - (1.0 / 36 + 4.0 / 9 + 1.0 / 4)) - y[2] * (9.0 / 4 + 9.0 / 4 - 3.0 * 2 / 4) - y[3] * (8.0 - 4.0 * 2 / 3)) * y[6];

    dydx[4] = cpi * (ytaubeta + cpi * y[4] * (-10 * ytau2 * ytau2 - 9 * yb2 * yb2 - 9 * yb2 * ytau2 - 3 * yb2 * ytop2 + (6 * y[2] + 6 * y[1] / 5) * ytau2 + (-2 * y[1] / 5 + 16 * y[3]) * yb2 + (9 * nf / 5 + 27.e0 / 10) * y[1] * y[1] + (3 * nf - 21.0 / 2) * y[2] * y[2] + 9 * y[1] * y[2] / 5));
    dydx[5] = cpi * (ybbeta + cpi * y[5] * (-22 * yb2 * yb2 - 5 * ytop2 * ytop2 - 5 * yb2 * ytop2 - 3 * yb2 * ytau2 - 3 * ytau2 * ytau2 + 4 * y[1] / 5 * ytop2 + (2 * y[1] / 5 + 6 * y[2] + 16 * y[3]) * yb2 + 6 * y[1] / 5 * ytau2 + (7 * nf / 15 + 7.0 / 18) * y[1] * y[1] + (3 * nf - 21.0 / 2) * y[2] * y[2] + (16 * nf / 3 - 304.0 / 9) * y[3] * y[3] + y[1] * y[2] + 8 * y[1] * y[3] / 9 + 8 * y[2] * y[3]));
    dydx[6] = cpi * (ytbeta + cpi * y[6] * (-22 * ytop2 * ytop2 - 5 * yb2 * yb2 - 5 * yb2 * ytop2 - yb2 * ytau2 + (6 * y[1] / 5 + 6 * y[2] + 16 * y[3]) * ytop2 + 2 * y[1] / 5 * yb2 + (13 * nf / 15 + 403.0 / 450) * y[1] * y[1] + (3 * nf - 21.0 / 2) * y[2] * y[2] + (16 * nf / 3 - 304.0 / 9) * y[3] * y[3] + y[1] * y[2] + 136.0 / 45 * y[1] * y[3] + 8 * y[2] * y[3]));

    //
    // (2-loop) higgs vev beta functions:
    // 	- y[7], y[8] = ln vu, ln vd
    // Ref : Check the the 2nd referance Eq(A6), Eq(A12), Eq(A13), Eq(A32), Eq(A33)
    // ------------------------------------------------------------------------------------------------------------------------------------------------

    dydx[7] = cpi * (3.0 / 4 * (y[1] / 5 + y[2]) - 3 * ytop2) + cpi * cpi * (9 * ytop2 * ytop2 / 4 + 9 * ytop2 * yb2 / 4 - (19 * y[1] / 10 + 9 * y[2] / 2 + 20 * y[3]) * ytop2 - (279.0 / 800 + 1803 * nf / 3200) * y[1] * y[1] - (207.0 / 32 + 357 * nf / 128) * y[2] * y[2] - (27.0 / 80 + 9 * nf / 160) * y[1] * y[2]);
    dydx[8] = cpi * (3.0 / 4 * (y[1] / 5 + y[2]) - 3 * yb2 - ytau2) + cpi * cpi * (9 * yb2 * yb2 / 4 + 9 * yb2 * ytop2 / 4 + 3 * ytau2 * ytau2 / 4 - (2 * y[1] / 5 + 9 * y[2] / 2 + 20 * y[3]) * yb2 - (9 * y[1] / 5 + 3 * y[2] / 2) * ytau2 - (279.0 / 800 + 1803 * nf / 3200) * y[1] * y[1] - (207.0 / 32 + 357 * nf / 128) * y[2] * y[2] - (27.0 / 80 + 9 * nf / 160) * y[1] * y[2]);

    //
    // (2-loop) soft susy-breaking terms beta functions:
    // 	- y[9]--y[11] : atau, ab, atop
    // Ref : Eq(4.22), Eq(4.23)
    // Ref : Eq(4.20), Eq(4.21)
    // Ref : Eq(4.18), Eq(4.19)
    // ------------------------------------------------------------------------------------------------------------------------------------------------
    double temp = 0.0;
    double dhtau1loop = y[9] * (3 * yb2 + 6 * ytau2 - 3 * y[2] - 9 * y[1] / 5) + 6 * y[10] * yb2 + 6 * y[9] * ytau2 + 6 * mm2 * y[2] + 18.0 * y[1] / 5 * mm1;

    temp = (-9 * yb2 * yb2 - 3 * yb2 * ytop2 - 14 * ytau2 * ytau2 - 15 * yb2 * ytau2 + (16 * y[3] - 2 * y[1] / 5) * yb2 + 6 * y[1] / 5 * ytau2);
    temp += +(12.0 * y[2] - 6.0 * y[1] / 5.0) * ytau2 + 15.0 * y[2] * y[2] / 2.0 + 9.0 * y[2] * y[1] / 5.0 + 27.0 * y[1] * y[1] / 2.0;
    double dhtau2loop = y[9] * temp;
    dhtau2loop += -6 * (6 * y[10] * yb2 * yb2 + (y[10] + y[11]) * yb2 * ytop2) - 36 * y[9] * ytau2 * ytau2 - yb2 * ytau2 * (12 * y[9] + 18 * y[10]) + (32.0 * y[3] - 4.0 * y[1] / 5.0) * yb2 * y[10];
    dhtau2loop += +12.0 * y[1] / 5.0 * ytau2 * y[9] + (6.0 * y[2] + 6.0 * y[1] / 5.0) * ytau2 * y[9] - (32.0 * y[3] * mm3 - 4.0 * y[1] / 5.0 * mm1) * yb2;
    dhtau2loop += -12.0 * y[1] / 5.0 * mm1 * ytau2 - 12.0 * y[2] * mm2 * ytau2 - 30.0 * y[2] * y[2] * mm2 - 18.0 * y[2] * y[1] / 5.0 * (mm1 + mm2) - 54.0 * y[1] * y[1] * mm1;
    dydx[9] = cpi * dhtau1loop + cpi * cpi * dhtau2loop - y[9] / y[4] * dydx[4];

    double dhb1loop = y[10] * (8 * yb2 + ytau2 + ytop2 - 16 * y[3] / 3 - 3 * y[2] - 7 * y[1] / 15) + 10 * y[10] * yb2;
    dhb1loop += +2 * y[9] * ytau2 + 2 * y[11] * ytop2 + 14 * y[1] / 15 * mm1 + 6 * y[2] * mm2 + 32 * y[3] / 3 * mm3;

    temp = 0.0;
    temp = -30 * yb2 * yb2 - 5 * ytop2 * ytop2 - 7 * yb2 * ytop2 - 5 * yb2 * ytau2 - 3 * ytau2 * ytau2;
    temp += +(16 * y[3] - 2 * y[1] / 5) * yb2 + 6 * y[1] / 5 * ytau2 + 4 * y[1] / 5 * ytop2 + (12 * y[2] + 6 * y[1] / 5) * yb2;
    temp += -16 * y[3] * y[3] / 9 + 8 * y[3] * y[2] / 9 + 15 * y[2] * y[2] / 2 + y[1] * y[2] + 287 * y[1] * y[1] / 90;
    double dhb2loop = y[10] * temp;
    dhb2loop += -80 * y[10] * yb2 * yb2 - 20 * y[11] * ytop2 * ytop2 - (8 * y[10] + 10 * y[11]) * yb2 * ytop2;
    dhb2loop += -12 * y[9] * ytau2 * ytau2 - yb2 * ytau2 * (6 * y[9] + 4 * y[10]) + (32 * y[3] - 4 * y[1] / 5) * yb2 * y[10] + 12 * y[1] / 5 * ytau2 * y[9];
    dhb2loop += +8 * y[1] / 5 * y[11] * ytop2 + (6 * y[2] + 6 * y[1] / 5) * yb2 * y[10] - (32 * y[3] * mm3 - 4 * y[1] / 5 * mm1) * yb2 - 12 * y[1] / 5 * mm1 * ytau2;
    dhb2loop += -(12 * y[2] * mm2 + 8 * y[1] / 5 * mm1) * yb2 - 8 * y[1] / 5 * mm1 * ytop2 + 64 * y[3] * y[3] / 9 * mm3 - 16 * y[3] * y[2] * (mm3 + mm2);
    dhb2loop += -16 * y[3] * y[1] / 9 * (mm3 + mm1) - 30 * y[2] * y[2] * mm2 - 2 * y[2] * y[1] * (mm1 + mm2) - 574 * y[1] * y[1] / 45 * mm1;

    dydx[10] = cpi * dhb1loop + cpi * cpi * dhb2loop - y[10] / y[5] * dydx[5];

    double dht1loop = y[11] * (8 * ytop2 + yb2 - 16 * y[3] / 3 - 3 * y[2] - 13 * y[1] / 15) + 10 * y[11] * ytop2 + 2 * y[10] * yb2 + 26 * y[1] / 15 * mm1 + 6 * y[2] * mm2 + 32 * y[3] / 3 * mm3;

    temp = -5 * yb2 * yb2 - 30 * ytop2 * ytop2 - 7 * yb2 * ytop2 - yb2 * ytau2 + (16 * y[3] + 4 * y[1] / 5) * ytop2 + 12 * y[2] * ytop2 + 2 * y[1] / 5 * yb2;
    temp += -16 * y[3] * y[3] / 9 + 8 * y[3] * y[2] + 15 * y[2] * y[2] / 2 + 136 * y[3] * y[1] / 45 + y[2] * y[1] + 2743 * y[1] * y[1] / 450;
    double dht2loop = y[11] * temp;
    dht2loop += -80 * y[11] * ytop2 * ytop2 - 20 * y[10] * yb2 * yb2 - (8 * y[11] + 10 * y[10]) * yb2 * ytop2 - 2 * yb2 * ytau2 * (y[9] + y[10]);
    dht2loop += +(32 * y[3] + 8 * y[1] / 5) * ytop2 * y[11] + 4 * y[1] / 5 * yb2 * y[10] + (6 * y[2] + 6 * y[1] / 5) * ytop2 * y[11];
    dht2loop += -(32 * y[3] * mm3 + 8 * y[1] / 5 * mm1) * ytop2 - 4 * y[1] / 5 * mm1 * yb2 - (12 * y[2] * mm2 + 4 * y[1] / 5 * mm1) * ytop2;
    dht2loop += +64 * y[3] * y[3] / 9 * mm3 - 16 * y[3] * y[2] * (mm3 + mm2) - 272 * y[3] * y[1] / 45 * (mm3 + mm1);
    dht2loop += -30 * y[2] * y[2] * mm2 - 2 * y[2] * y[1] * (mm1 + mm2) - 5486 * y[1] * y[1] / 225 * mm1;

    dydx[11] = cpi * dht1loop + cpi * cpi * dht2loop - y[11] / y[6] * dydx[6];

    //
    // Scalar phi_u mass term squared :
    // 	y[12]--y[13] : m^2(phi_u), m^2(phi_d)
    //	Eq(4.27), Eq(4.28), Eq(4.29), Eq(4.30), Eq(4.31), Eq(4.32), Eq(4.27), Eq(4.27)
    // ------------------------------------------------------------------------------------------------------------------------------------------------

    //
    // One-loop part
    // ------------------------------------------------------------------------------------------------------------------------------------------------

    double trym2;
    trym2 = y[18] - 2 * y[17] + y[16] - y[15] + y[14] + y[12] - y[13] + 2 * (y[28] - 2 * y[27] + y[26] - y[25] + y[24]);
    dydx[12] = 2 * cpi * (3 * ytop2 * (y[12] + y[18] + y[17] + y[11] * y[11]) + 3.0 / 10 * y[1] * trym2 - 3 * y[1] / 5 * pow(mm1, 2) - 3 * y[2] * pow(mm2, 2));
    dydx[13] = 2 * cpi * (ytau2 * (y[13] + y[15] + y[14] + y[9] * y[9]) + 3 * yb2 * (y[13] + y[18] + y[16] + y[10] * y[10]) - 3.0 / 10 * y[1] * trym2 - 3 * y[1] / 5 * pow(mm1, 2) - 3 * y[2] * pow(mm2, 2));

    //	printf(" (%e %e %e) ",trym2,y[1],y[2]);
    //
    // Two-loop part
    // ------------------------------------------------------------------------------------------------------------------------------------------------

    double sumt = y[12] + y[18] + y[17] + y[11] * y[11];
    double sumb = y[13] + y[18] + y[16] + y[10] * y[10];
    double sumtau = y[13] + y[15] + y[14] + y[9] * y[9];

    double trmq = 2 * y[28] + y[18]; // Yu'Yu
    double trmu = 2 * y[27] + y[17]; // Yul'Yul
    double trmd = 2 * y[26] + y[16]; // Yd'Yd
    double trml = 2 * y[25] + y[15]; // Yel'Yel
    double trme = 2 * y[24] + y[14]; // Ye'Ye

    double curlysp = 0.0;
    curlysp += -ytop2 * (3 * y[12] + y[18] - 4 * y[17]);
    curlysp += +yb2 * (3 * y[13] - y[18] - 2 * y[16]);
    curlysp += +ytau2 * (y[13] + y[15] - 2 * y[14]);
    curlysp += +(1.5 * y[2] + 0.3 * y[1]) * (y[12] - y[13] - trml);
    curlysp += +(8.0 / 3 * y[3] + 1.5 * y[2] + 1.0 / 30 * y[1]) * trmq;
    curlysp += -(16.0 / 3 * y[3] + 16.0 / 15 * y[1]) * trmu;
    curlysp += +(8.0 / 3 * y[3] + 2.0 / 15 * y[1]) * trmd + 1.2 * y[1] * trme;

    double sig1 = y[1] / 5 * (3 * (y[12] + y[13]) + trmq + 3 * trml + 8 * trmu + 2 * trmd + 6 * trme);
    double sig2 = y[2] * (y[12] + y[13] + 3 * trmq + trml);
    double sig3 = y[3] * (2 * trmq + trmu + trmd);

    temp = 0.0;                                                                                                    //	Wq(4.31)
    temp += -6.0 * (6 * pow(ytop2, 2) * (sumt + y[11] * y[11]) + (sumt + sumb + 2 * y[11] * y[10]) * ytop2 * yb2); // The term inside the third bracket
    temp += +32 * y[3] * ytop2 * (sumt + 2.0 * pow(mm3, 2) - 2 * y[11] * mm3);
    temp += +1.6 * y[1] * ytop2 * (sumt + 2.0 * pow(mm1, 2) - 2 * y[11] * mm1);
    temp += +1.2 * y[1] * curlysp + 33 * pow(y[2], 2) * pow(mm2, 2);
    temp += +3.6 * y[1] * y[2] * (pow(mm2, 2) + pow(mm1, 2) + mm2 * mm1);
    temp += +621.0 / 25 * pow(y[1], 2) * pow(mm1, 2) + 3 * y[2] * sig2 + 0.6 * y[1] * sig1;
    dydx[12] = dydx[12] + pow(cpi, 2) * temp;

    temp = 0.0;                                                                                                                                                 // Eq(4.33)
    temp += -6.0 * (6 * pow(yb2, 2) * (sumb + y[10] * y[10]) + (sumt + sumb + 2.0 * y[11] * y[10]) * ytop2 * yb2 + 2 * (sumtau + y[9] * y[9]) * pow(ytau2, 2)); // The term inside the third bracket
    temp += +32 * y[3] * yb2 * (sumb + 2 * pow(mm3, 2) - 2 * y[10] * mm3);
    temp += -0.8 * y[1] * yb2 * (sumb + 2 * pow(mm1, 2) - 2 * y[10] * mm1);
    temp += +2.4 * y[1] * ytau2 * (sumtau + 2.0 * pow(mm1, 2) - 2 * y[9] * mm1);
    temp += -1.2 * y[1] * curlysp + 33 * pow(y[2], 2) * pow(mm2, 2);
    temp += +3.6 * y[1] * y[2] * (pow(mm2, 2) + pow(mm1, 2) + mm2 * mm1);
    temp += +621.0 / 25 * pow(y[1], 2) * pow(mm1, 2) + 3 * y[2] * sig2 + 0.6 * y[1] * sig1;
    dydx[13] = dydx[13] + pow(cpi, 2) * temp;

    //
    // S-particles mass term squared (2-loop)
    // 	y[14]--y[19] : m^2_tau, m^2_l, m^2_b, m^2_top, m^2_q, b
    // Eq(4.42), Eq(4.43)
    // Eq(4.36), Eq(4.37)
    // Eq(4.40), Eq(4.41)
    // Eq(4.38), Eq(4.39)
    // Eq(4.34), Eq(4.35)
    // Eq(4.25), Eq(4.26)
    // ------------------------------------------------------------------------------------------------------------------------------------------------

    dydx[14] = 2 * cpi * (2 * ytau2 * (y[13] + y[14] + y[15] + y[9] * y[9]) + 3 * y[1] / 5 * trym2 - 12 * y[1] / 5 * pow(mm1, 2));
    temp = -16 * pow(ytau2, 2) * (sumtau + y[9] * y[9]);
    temp += -ytau2 * yb2 * (12 * (sumtau + sumb) + 8 * y[9] * y[10]) + 2 * sumtau * ytau2 * (6 * y[2] - 6 * y[1] / 5);
    temp += +12 * y[2] * ytau2 * (2 * pow(mm2, 2) - 2 * y[9] * mm2) - 12 * y[1] / 5 * ytau2 * (2 * pow(mm1, 2) - 2 * y[9] * mm1);
    temp += +12 * y[1] / 5 * curlysp + 2808.0 / 25 * pow(y[1], 2) * pow(mm1, 2) + 12 * y[1] / 5 * sig1;
    dydx[14] = dydx[14] + pow(cpi, 2) * temp;

    dydx[15] = 2 * cpi * (ytau2 * (y[13] + y[15] + y[14] + y[9] * y[9]) - 3 * y[1] / 10 * trym2 - 3 * y[1] / 5 * pow(mm1, 2) - 3 * y[2] * pow(mm2, 2));
    temp = -12 * pow(ytau2, 2) * (sumtau + y[9] * y[9]);
    temp += -6 * ytau2 * yb2 * (sumtau + sumb + 2 * y[9] * y[10]) + 6 * y[1] / 5 * ytau2 * (2 * sumtau - 4 * mm1 * y[9] + 4 * pow(mm1, 2));
    temp += -6 * y[1] / 5 * curlysp + 33 * pow(y[2], 2) * pow(mm2, 2) + 18 * y[1] * y[2] / 5 * (pow(mm2, 2) + pow(mm1, 2) + mm1 * mm2);
    temp += +621.0 / 25 * pow(y[1], 2) * pow(mm1, 2) + 3 * y[1] / 5 * sig1 + 3 * y[2] * sig2;
    dydx[15] = dydx[15] + pow(cpi, 2) * temp;

    dydx[16] = 2 * cpi * (2 * yb2 * (y[13] + y[16] + y[18] + y[10] * y[10]) + y[1] / 5 * trym2 - 4 * y[1] / 15 * pow(mm1, 2) - 16 * y[3] / 3 * pow(mm3, 2));
    temp = -32 * pow(yb2, 2) * (sumb + y[10] * y[10]) - 4 * ytop2 * yb2 * (sumt + sumb + 2 * y[11] * y[10]);
    temp += -4 * ytau2 * yb2 * (sumtau + sumb + 2 * y[9] * y[10]) + 2 * (6 * y[2] + 2 * y[1] / 5) * yb2 * sumb;
    temp += +12 * y[2] * yb2 * (2 * pow(mm2, 2) - 2 * y[10] * mm2) + 4 * y[1] / 5 * yb2 * (2 * pow(mm1, 2) - 2 * y[10] * mm1);
    temp += +4 * y[1] / 5 * curlysp - 128 * pow(y[3], 2) / 3 * pow(mm3, 2) + 128 * y[1] * y[3] / 45 * (pow(mm3, 2) + pow(mm1, 2) + mm1 * mm3);
    temp += +808.0 / 75 * pow(y[1], 2) * pow(mm1, 2) + 4 * y[1] / 15 * sig1 + 16 * y[3] / 3 * sig3;
    dydx[16] = dydx[16] + pow(cpi, 2) * temp;

    dydx[17] = 2 * cpi * (2 * ytop2 * (y[12] + y[17] + y[18] + y[11] * y[11]) - 2 * y[1] / 5 * trym2 - 16 * y[1] / 15 * pow(mm1, 2) - 16 * y[3] / 3 * pow(mm3, 2));
    temp = -32 * pow(ytop2, 2) * (sumt + y[11] * y[11]);
    temp += -4 * ytop2 * yb2 * (sumt + sumb + 2 * y[11] * y[10]) + 2 * (6 * y[2] - 2 * y[1] / 5) * ytop2 * sumt;
    temp += +12 * y[2] * ytop2 * (2 * pow(mm2, 2) - 2 * y[11] * mm2) - 4 * y[1] / 5 * ytop2 * (2 * pow(mm1, 2) - 2 * y[11] * mm1);
    temp += -8 * y[1] / 5 * curlysp - 128 * pow(y[3], 2) / 3 * pow(mm3, 2) + 512 * y[1] * y[3] / 45 * (pow(mm3, 2) + pow(mm1, 2) + mm1 * mm3);
    temp += +3424.0 / 75 * pow(y[1], 2) * pow(mm1, 2) + 16 * y[1] / 15 * sig1 + 16 * y[3] / 3 * sig3;
    dydx[17] = dydx[17] + pow(cpi, 2) * temp;

    dydx[18] = 2 * cpi * (ytop2 * (y[12] + y[17] + y[18] + y[11] * y[11]) + yb2 * (y[13] + y[18] + y[16] + y[10] * y[10]) + y[1] / 10 * trym2 - y[1] / 15 * pow(mm1, 2) - 3 * y[2] * pow(mm2, 2) - 16 * y[3] / 3 * pow(mm3, 2));
    temp = -20 * pow(ytop2, 2) * (sumt + y[11] * y[11]) - 20 * pow(yb2, 2) * (sumb + y[10] * y[10]) - 2 * ytau2 * yb2 * (sumtau + sumb + y[9] * y[10]);
    temp += +2 * y[1] / 5 * (4 * ytop2 * (sumt - 2 * mm1 * y[11] + 2 * pow(mm1, 2)) + 2 * yb2 * (sumb - 2 * mm1 * y[10] + 2 * pow(mm1, 2)));
    temp += +2 * y[1] / 5 * curlysp - 128 * pow(y[3], 2) / 3 * pow(mm3, 2);
    temp += +32 * y[2] * y[3] * (pow(mm3, 2) + pow(mm2, 2) + mm2 * mm3) + 32 * y[1] * y[3] / 45 * (pow(mm3, 2) + pow(mm1, 2) + mm1 * mm3);
    temp += +33 * pow(y[2], 2) * pow(mm2, 2) + 2 * y[1] * y[2] / 5 * (pow(mm2, 2) + pow(mm1, 2) + mm1 * mm2);
    temp += +199 * pow(y[1], 2) / 75 * pow(mm1, 2) + y[1] / 15 * sig1 + 16 * y[3] / 3 * sig3 + 3 * y[2] * sig2;
    dydx[18] = dydx[18] + pow(cpi, 2) * temp;

    double betb2 = 0.0;
    dydx[19] = 2 * cpi * (3 * y[11] * ytop2 + 3 * y[10] * yb2 + y[9] * ytau2 + 3 * y[1] / 5 * mm1 + 3 * y[2] * mm2);
    betb2 += -12 * (3 * y[11] * pow(ytop2, 2) + 3 * y[10] * pow(yb2, 2) + y[11] * yb2 * ytop2 + y[10] * ytop2 * yb2 + y[9] * pow(ytau2, 2));
    betb2 += +(32 * y[3] + 8 * y[1] / 5) * y[11] * ytop2 + (32 * y[3] - 4 * y[1] / 5) * y[10] * yb2 + 12 * y[1] / 5 * y[9] * ytau2;
    betb2 += -(32 * y[3] * mm3 + 8 * y[1] * mm1 / 5) * ytop2 - (32 * y[3] * mm3 - 4 * y[1] * mm1 / 5) * yb2 - 12 * y[1] / 5 * mm1 * ytau2;
    betb2 += -30 * pow(y[2], 2) * mm2 - 18 * y[1] * y[2] / 5 * (mm1 + mm2) - 414 * pow(y[1], 2) / 25 * mm1;
    dydx[19] = dydx[19] + pow(cpi, 2) * betb2;

    //
    // Gauginos masses beta functions (includes two-loop):
    // 	- y[20]--y[22] : ln (m1,m2,m3)
    // Ref : Referance 2 Eq(A26), Eq(A27), Eq(A28)
    // ------------------------------------------------------------------------------------------------------------------------------------------------

    temp = (19.0 * nf / 15 + 9.0 / 25) * y[1] * (1.0 + 1.0) + (3.0 * nf / 5 + 9.0 / 5) * y[2] * (1.0 + mm2 / mm1) + 44.0 * nf / 15 * y[3] * (1.0 + mm3 / mm1);
    temp += -26 * ytop2 * (1.0 - y[11] / mm1) / 5 - 14 * yb2 * (1.0 - y[10] / mm1) / 5 - 18 * ytau2 * (1.0 - y[9] / mm1) / 5;
    dydx[20] = -2 * cpi * (-3.0 / 5 - nf) * y[1] + 2 * cpi * cpi * y[1] * temp;

    temp = (nf / 5.0 + 3.0 / 5) * y[1] * (1.0 + mm1 / mm2) + (7.0 * nf - 17.0) * y[2] * (1.0 + 1.0);
    temp += +4 * nf * y[3] * (1.0 + mm3 / mm2) - 6 * ytop2 * (1.0 - y[11] / mm2);
    temp += -6 * yb2 * (1.0 - y[10] / mm2) - 2 * ytau2 * (1.0 - y[9] / mm2);
    dydx[21] = -2 * cpi * (5.0 - nf) * y[2] + 2 * cpi * cpi * y[2] * temp;

    temp = 11.0 * nf / 30 * y[1] * (1.0 + mm1 / mm3) + 3.0 * nf / 2 * y[2] * (1.0 + mm2 / mm3);
    temp += +(34.0 * nf / 3 - 54.0) * y[3] * (1.0 + 1.0) - 4 * ytop2 * (1.0 - y[11] / mm3) - 4 * yb2 * (1.0 - y[10] / mm3);
    dydx[22] = -2 * cpi * (9.0 - nf) * y[3] + 2 * cpi * cpi * y[3] * temp;

    //
    // The mu parameter:
    // 	- y[23] = ln mu
    // Ref : Eq(4.25), Eq(4.26)
    // ------------------------------------------------------------------------------------------------------------------------------------------------

    dydx[23] = cpi * (3 * ytop2 + 3 * yb2 + ytau2 - 3 * y[1] / 5 - 3 * y[2]);

    //
    //     two-loop part
    // ------------------------------------------------------------------------------------------------------------------------------------------------

    temp = -3 * (3 * pow(ytop2, 2) + 3 * pow(yb2, 2) + 2 * ytop2 * yb2 + pow(ytau2, 2));
    temp += +(16 * y[3] + 4.0 / 5. * y[1]) * ytop2 + (16 * y[3] - 2.0 / 5. * y[1]) * yb2 + 6.0 / 5. * y[1] * ytau2;
    temp += +7.5 * pow(y[2], 2) + 1.8 * y[1] * y[2] + 207.0 / 50. * pow(y[1], 2);
    dydx[23] = dydx[23] + pow(cpi, 2) * temp;

    //
    // (2-loop) y(24)--y(28) : 1st and 2d gen. sfermion mass^2 terms:
    // 	m^2_er, m^2_el, m^2_dr, m^2_ur, m^2_ul
    // Ref :
    //		Eq(4.42), Eq(4.43)
    //		Eq(4.36), Eq(4.37)
    //		Eq(4.40), Eq(4.41)
    //		Eq(4.38), Eq(4.39)
    //		Eq(4.34), Eq(4.35)
    // 																			..  Check why other terms are missing
    // ------------------------------------------------------------------------------------------------------------------------------------------------

    dydx[24] = 2 * cpi * (3 * y[1] / 5 * trym2 - 12 * y[1] / 5 * pow(mm1, 2));
    dydx[24] = dydx[24] + pow(cpi, 2) * (12 * y[1] / 5 * curlysp + 2808.0 / 25 * pow(y[1], 2) * pow(mm1, 2) + 12 * y[1] / 5 * sig1);

    dydx[25] = 2 * cpi * (-3 * y[1] / 10 * trym2 - 3 * y[1] / 5 * pow(mm1, 2) - 3 * y[2] * pow(mm2, 2));
    dydx[25] = dydx[25] + pow(cpi, 2) * (-6.0 * y[1] / 5 * curlysp + 33.0 * pow(y[2], 2) * pow(mm2, 2) + 18.0 * y[1] * y[2] / 5 * (pow(mm2, 2) + pow(mm1, 2) + mm1 * mm2) + 621.0 / 25 * pow(y[1], 2) * pow(mm1, 2) + 3.0 * y[1] / 5 * sig1 + 3.0 * y[2] * sig2);

    dydx[26] = 2 * cpi * (y[1] / 5 * trym2 - 4 * y[1] / 15 * pow(mm1, 2) - 16 * y[3] / 3 * pow(mm3, 2));
    temp = 4 * y[1] / 5 * curlysp - 128 * pow(y[3], 2) / 3 * pow(mm3, 2) + 128 * y[1] * y[3] / 45 * (pow(mm3, 2) + pow(mm1, 2) + mm1 * mm3);
    temp += +808.0 / 75 * pow(y[1], 2) * pow(mm1, 2) + 4 * y[1] / 15 * sig1 + 16 * y[3] / 3 * sig3;
    dydx[26] = dydx[26] + pow(cpi, 2) * temp;

    dydx[27] = 2 * cpi * (-2 * y[1] / 5 * trym2 - 16 * y[1] / 15 * pow(mm1, 2) - 16 * y[3] / 3 * pow(mm3, 2));
    temp = -8 * y[1] / 5 * curlysp - 128 * pow(y[3], 2) / 3 * pow(mm3, 2) + 512 * y[1] * y[3] / 45 * (pow(mm3, 2) + pow(mm1, 2) + mm1 * mm3);
    temp += +3424.0 / 75 * pow(y[1], 2) * pow(mm1, 2) + 16 * y[1] / 15 * sig1 + 16 * y[3] / 3 * sig3;
    dydx[27] = dydx[27] + pow(cpi, 2) * temp;

    dydx[28] = 2 * cpi * (y[1] / 10 * trym2 - y[1] / 15 * pow(mm1, 2) - 3 * y[2] * pow(mm2, 2) - 16 * y[3] / 3 * pow(mm3, 2));
    temp = 2 * y[1] / 5 * curlysp - 128 * pow(y[3], 2) / 3 * pow(mm3, 2);
    temp += +32 * y[2] * y[3] * (pow(mm3, 2) + pow(mm2, 2) + mm2 * mm3);
    temp += +32 * y[1] * y[3] / 45 * (pow(mm3, 2) + pow(mm1, 2) + mm1 * mm3);
    temp += +33 * pow(y[2], 2) * pow(mm2, 2);
    temp += +2 * y[1] * y[2] / 5 * (pow(mm2, 2) + pow(mm1, 2) + mm1 * mm2);
    temp += +199 * pow(y[1], 2) / 75 * pow(mm1, 2);
    temp += +y[1] / 15 * sig1 + 16 * y[3] / 3 * sig3 + 3 * y[2] * sig2;
    dydx[28] = dydx[28] + pow(cpi, 2) * temp;

    //
    // Trilinear coupling for muon, strange and charm
    // 	- (2-loop) y(29)--y(31) : ae (anu), ad (as), au (ac)
    // Ref :
    // 	Eq(4.22), Eq(4.23)
    // 	Eq(4.20), Eq(4.21)
    // 	Eq(4.18), Eq(4.19)
    // ------------------------------------------------------------------------------------------------------------------------------------------------

    double dhe1loop = y[29] * (ytau2 + 3 * yb2 - 3 * y[2] - 9 * y[1] / 5) + 6 * y[10] * yb2 + 2 * y[9] * ytau2 + 6 * mm2 * y[2] + 18 * y[1] / 5 * mm1;
    double dhe2loop = y[29] * (-9 * yb2 * yb2 - 3 * yb2 * ytop2 - 3 * ytau2 * ytau2 + (16 * y[3] - 2 * y[1] / 5) * yb2 + 6 * y[1] / 5 * ytau2 + 15 * y[2] * y[2] / 2 + 9 * y[2] * y[1] / 5 + 27 * y[1] * y[1] / 2);
    dhe2loop += -6 * (6 * y[10] * yb2 * yb2 + (y[10] + y[11]) * yb2 * ytop2);
    dhe2loop += -12 * y[9] * ytau2 * ytau2;
    dhe2loop += +(32 * y[3] - 4 * y[1] / 5) * yb2 * y[10] + 12 * y[1] / 5 * ytau2 * y[9];
    dhe2loop += -(32 * y[3] * mm3 - 4 * y[1] / 5 * mm1) * yb2 - 12 * y[1] / 5 * mm1 * ytau2;
    dhe2loop += -30 * y[2] * y[2] * mm2 - 18 * y[2] * y[1] / 5 * (mm1 + mm2) - 54 * y[1] * y[1] * mm1;
    temp = -3 * ytau2 * ytau2 - 9 * yb2 * yb2 - 3 * yb2 * ytop2 + 6 * y[1] / 5 * ytau2 + (-2 * y[1] / 5 + 16 * y[3]) * yb2;
    temp += +(9 * nf / 5 + 27.0 / 10) * y[1] * y[1] + (3 * nf - 21.0 / 2) * y[2] * y[2] + 9 * y[1] * y[2] / 5;
    double dyovery4 = cpi * (ytau2 + 3 * yb2 - 9 * y[1] / 5 - 3 * y[2]) + pow(cpi, 2) * temp;
    dydx[29] = cpi * dhe1loop + cpi * cpi * dhe2loop - y[29] * dyovery4;

    double dhd1loop = y[30] * (3 * yb2 + ytau2 - 16 * y[3] / 3 - 3 * y[2] - 7 * y[1] / 15) + 6 * y[10] * yb2 + 2 * y[9] * ytau2 + 14 * y[1] / 15 * mm1 + 6 * y[2] * mm2 + 32 * y[3] / 3 * mm3;
    temp = -9 * yb2 * yb2 - 3 * yb2 * ytop2 - 3 * ytau2 * ytau2 + (16 * y[3] - 2 * y[1] / 5) * yb2;
    temp += +6 * y[1] / 5 * ytau2 - 16 * y[3] * y[3] / 9 + 8 * y[3] * y[2] / 9 + 15 * y[2] * y[2] / 2 + y[1] * y[2] + 287 * y[1] * y[1] / 90;
    double dhd2loop = y[30] * temp;
    dhd2loop += -36 * y[10] * yb2 * yb2 - 6 * (y[10] + y[11]) * yb2 * ytop2;
    dhd2loop += -12 * y[9] * ytau2 * ytau2 + (32 * y[3] - 4 * y[1] / 5) * yb2 * y[10] + 12 * y[1] / 5 * ytau2 * y[9] - (32 * y[3] * mm3 - 4 * y[1] / 5 * mm1) * yb2 - 12 * y[1] / 5 * mm1 * ytau2;
    dhd2loop += +64 * y[3] * y[3] / 9 * mm3 - 16 * y[3] * y[2] * (mm3 + mm2) - 16 * y[3] * y[1] / 9 * (mm3 + mm1) - 30 * y[2] * y[2] * mm2 - 2 * y[2] * y[1] * (mm1 + mm2) - 574 * y[1] * y[1] / 45 * mm1;
    temp = -9 * yb2 * yb2 - 3 * yb2 * ytop2 - 3 * ytau2 * ytau2 + (-2 * y[1] / 5 + 16 * y[3]) * yb2 + 6 * y[1] / 5 * ytau2;
    temp += +(7 * nf / 15 + 7.0 / 18) * y[1] * y[1] + (3 * nf - 21.0 / 2) * y[2] * y[2];
    temp += +(16 * nf / 3 - 304.0 / 9) * y[3] * y[3] + y[1] * y[2] + 8 * y[1] * y[3] / 9 + 8 * y[2] * y[3];
    double dyovery5 = cpi * (3 * yb2 + ytau2 - 7 * y[1] / 15 - 3 * y[2] - 16 * y[3] / 3) + pow(cpi, 2) * temp;
    dydx[30] = cpi * dhd1loop + cpi * cpi * dhd2loop - y[30] * dyovery5;

    double dhu1loop = y[31] * (3 * ytop2 - 16 * y[3] / 3 - 3 * y[2] - 13 * y[1] / 15) + 6 * y[11] * ytop2 + 26 * y[1] / 15 * mm1 + 6 * y[2] * mm2 + 32 * y[3] / 3 * mm3;
    temp = -9 * ytop2 * ytop2 - 3 * yb2 * ytop2 + (16 * y[3] + 4 * y[1] / 5) * ytop2 - 16 * y[3] * y[3] / 9;
    temp += +8 * y[3] * y[2] + 15 * y[2] * y[2] / 2 + 136 * y[3] * y[1] / 45 + y[2] * y[1] + 2743 * y[1] * y[1] / 450;
    double dhu2loop = y[31] * temp;
    dhu2loop += -36 * y[11] * ytop2 * ytop2 - 6 * (y[11] + y[10]) * yb2 * ytop2 + (32 * y[3] + 8 * y[1] / 5) * ytop2 * y[11];
    dhu2loop += -(32 * y[3] * mm3 + 8 * y[1] / 5 * mm1) * ytop2 + 64 * y[3] * y[3] / 9 * mm3 - 16 * y[3] * y[2] * (mm3 + mm2);
    dhu2loop += -272 * y[3] * y[1] / 45 * (mm3 + mm1) - 30 * y[2] * y[2] * mm2 - 2 * y[2] * y[1] * (mm1 + mm2) - 5486 * y[1] * y[1] / 225 * mm1;
    temp = -9 * ytop2 * ytop2 - 3 * yb2 * ytop2 + (4 * y[1] / 5 + 16 * y[3]) * ytop2 + (13 * nf / 15 + 403.0 / 450) * y[1] * y[1];
    temp += +(3 * nf - 21.e0 / 2) * y[2] * y[2] + (16 * nf / 3 - 304.0 / 9) * y[3] * y[3] + y[1] * y[2] + 136.0 / 45 * y[1] * y[3] + 8 * y[2] * y[3];
    double dyovery6 = cpi * (3 * ytop2 - 13 * y[1] / 15 - 3 * y[2] - 16 * y[3] / 3) + pow(cpi, 2) * temp;
    dydx[31] = cpi * dhu1loop + cpi * cpi * dhu2loop - y[31] * dyovery6;

}

//
//  Fourth order runge--kutta numerical algorithms solving differential equations by numerical recipes.
//  Needed by the routines above for the rges.
//  ------------------------------------------------------------------------------------------------------------------------------------

void su_rk4(double y[], double dydx[], int n, double x, double h, double yout[], void su_derivs(double, double *, double *))
{
    static int nmax = 31;
    double yt[nmax + 1], dyt[nmax + 1], dym[nmax + 1];
    double hh, h6, xh;

    hh = h / 2.0;
    h6 = h / 6.0;
    xh = x + hh;
    for (int i = 1; i <= n; i++)
        yt[i] = y[i] + hh * dydx[i];

    su_derivs(xh, yt, dyt);
    for (int i = 1; i <= n; i++)
        yt[i] = y[i] + hh * dyt[i];

    su_derivs(xh, yt, dym);
    for (int i = 1; i <= n; i++)
    {
        yt[i] = y[i] + h * dym[i];
        dym[i] = dyt[i] + dym[i];
    }
    su_derivs(x + h, yt, dyt);

    for (int i = 1; i <= n; i++)
        yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2 * dym[i]);
}

//
// Fourth order runge-kutta numerical algorithms solving differential equations by numerical recipes
// Needed by the su_odeint above
// ---------------------------------------------------------------------------------------------------------------------------------------------

double su_rkqc(double y[], double dydx[], int n, double *x, double htry, double eps, double yscal[], double *hdid, double *hnext, int choice)
{
    static int nmax = 31;
    static double fcor = 0.066666666666666667;
    static double one = 1.0;
    static double safety = 0.90;
    static double errcon = 6.0e-6;

    double ysav[nmax + 1], dysav[nmax + 1], ytemp[nmax + 1];

    double pgrow = -0.20;
    double pshrnk = -0.25;
    double xsav = *x;

    double h, hh, errmax;

    void (*su_derivs)(double, double *, double *);

    if (choice == 1)
        su_derivs = su_deriv1;
    else if (choice == 2)
        su_derivs = su_deriv2;

    for (int i = 1; i <= n; i++)
    {
        ysav[i] = y[i];
        dysav[i] = dydx[i];
    }
    h = htry;

    while (1)
    {
        hh = h / 2;
        su_rk4(ysav, dysav, n, xsav, hh, ytemp, su_derivs);
        *x = xsav + hh;

        su_derivs(*x, ytemp, dydx);
        su_rk4(ytemp, dydx, n, *x, hh, y, su_derivs);

        *x = xsav + h;
        if (*x == xsav)
        {
            // pause 'stepsize not significant in rkqc.'
            // write(*,'(a)') 'stepsize not significant in rkqc.'
            iflop = 1;
            return 0;
        }

        su_rk4(ysav, dysav, n, xsav, h, ytemp, su_derivs);
        errmax = 0.0;
        for (int i = 1; i <= n; i++)
        {
            ytemp[i] = y[i] - ytemp[i];
            errmax = (errmax > fabs(ytemp[i] / yscal[i]))? errmax : fabs(ytemp[i] / yscal[i]);
        }
        errmax = errmax / eps;
        if (errmax > one)
        {
            h = safety * h * pow(errmax, pshrnk);

            continue;
        }
        else
        {
            *hdid = h;
            if (errmax > errcon)
                *hnext = safety * h * pow(errmax, pgrow);
            else
                *hnext = 4 * h;
            break;
        }
    }
    for (int i = 1; i <= n; i++)
        y[i] = y[i] + ytemp[i] * fcor;

    return 0;
}

//
// The following routines are for the rge evolution of the parameters  																													su_rkqc(y,		dydx,		nvar,x,		h1,	eps,yscal,	hdid,		hnext,su_derivs)
// ---------------------------------------------------------------------------------------------------------------------

int su_odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1, double hmin, int choice)
{
    //	y,dydx,nvar,x,h,eps,yscal,hdid,hnext,su_derivs
    // -----------------------------------------------------------------
    // this is the main subroutine (from numerical recipes) integrating
    // (coupled) ordinary differential equation

    int nok, nbad;
    static int maxstp = 100000;
    static int nmax = 31;
    static double two = 2.0;
    static double zero = 0.0;
    static double tiny = 1.e-30;
    double xsav;
    double yscal[nmax + 1], y[nmax + 1], dydx[nmax + 1];

    double x = x1;
    double h = h1 * (x2 - x1) / fabs(x2 - x1);
    double dxsav;
    double xp[201],yp[32][201];
    int kmax = 0, kount = 0;

    nok = 0;
    nbad = 0;


    void (*su_derivs)(double, double *, double *);

    if (choice == 1)
        su_derivs = su_deriv1;
    else if (choice == 2)
        su_derivs = su_deriv2;

    for (int i = 1; i <= nvar; i++)
        y[i] = ystart[i];

    xsav = x - dxsav * two;
    for (int nstp = 1; nstp <= maxstp; nstp++)
    {
        su_derivs(x, y, dydx);

        for (int i = 1; i <= nvar; i++)
            yscal[i] = fabs(y[i]) + fabs(h * dydx[i]) + tiny;

        if (kmax > 0)
        {
            if (fabs(x - xsav) > fabs(dxsav))
            {
                if (kount < kmax - 1)
                {
                    kount = kount + 1;
                    xp[kount] = x;
                    for (int i = 1; i <= nvar; i++)
                        yp[i][kount] = y[i];
                    xsav = x;
                }
            }
        }

        if ((x + h - x2) * (x + h - x1) > zero)
            h = x2 - x;

        // New modif to speed up rge integration in the safe zone far from gut:

        double hdid, hnext;

        if (x < log(1.0e14))
            su_rkqc(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, choice);
        else
            su_rkqc(y, dydx, nvar, &x, h1, eps, yscal, &hdid, &hnext, choice);

        if (hdid == h)
            nok = nok + 1;
        else
            nbad = nbad + 1;

        if ((x - x2) * (x2 - x1) >= zero)
        {
            for (int i = 1; i <= nvar; i++)
                ystart[i] = y[i];

            if (kmax != 0)
            {
                kount = kount + 1;
                xp[kount] = x;

                for (int i = 1; i <= nvar; i++)
                    yp[i][kount] = y[i];
            }
            return 0;
        }
        // if(fabs(hnext) < hmin) pause 'stepsize smaller than minimum.'
        iflop = 0;
        if (fabs(hnext) < hmin) // write(*,'(a)') 'stepsize smaller than minimum.'
            iflop = 1;

        h = hnext;
    }
    iflop = 1;
    return 0;
}

//
//  calculates the sfermion masses and the mixing angles at scale=mz
//  for the 3d generation sfermions in the bpmz conventions.
// -------------------------------------------------------------------------------------------------------

void su_sfbpmz(double pizz, double mql, double mur, double mdr, double mel, double mer, double mql1, double mur1, double mdr1, double mel1, double mer1, double al, double at, double ab, double mu, double b_mz, double tb, double rmtau, double rmb, double rmt)
{
    errma2z = 0.0;
    double b_new = atan(tb);
    double cw = sqrt(1.0 - sw2);
    double rmz = sqrt(pow(mz, 2) + pizz);
    double rmw = cw * rmz;
    double mzsave = mz;
    double mwsave = mw;
    mz = rmz;
    mw = rmw;
    double mz2 = pow(mz, 2);
    double ma2;

    if (ipolemz == 0)
    { // running higgs masses in loops at mz
        ma2 = mu * b_mz / sin(b_new) / cos(b_new);
        if (ma2 < 0.0 && irge < irgmax)
            ma2 = 0.10; // temp protection
        if (ma2 < 0.0)
            errma2z = -1.0; // final ma2<0->error flag
        ma = sqrt(fabs(ma2));
    }
    else if (ipolemz == 1)
    { // pole higgs masses in loops at mz
        ma2 = pow(mapole, 2);
        ma = sqrt(fabs(ma2));
    }

    marunz = ma;
    mchrunz = sqrt(fabs(ma2 + pow(mw, 2)));
    double mhht2 = 1.0 / 2 * (ma2 + mz2 + sqrt(pow(ma2 + mz2, 2) - pow(2 * ma * mz * cos(2.0 * b_new), 2)));
    double mht2 = 1.0 / 2 * (ma2 + mz2 - sqrt(pow(ma2 + mz2, 2) - pow(2 * ma * mz * cos(2.0 * b_new), 2)));
    mhrunz = sqrt(fabs(mhht2));
    mlrunz = sqrt(fabs(mht2));

    //	pi = 4*atan(1.0);

    double s2alt = -(ma2 + mz2) * sin(2.0 * b_new);
    double c2alt = -(ma2 - mz2) * cos(2.0 * b_new);
    double t2alt = s2alt / c2alt;

    if (c2alt > 0)
        alpharunz = 0.50 * atan(t2alt);
    else if (c2alt < 0)
        alpharunz = 0.50 * atan(t2alt) - pi / 2.0;
    else
        alpharunz = -pi / 4.0;

    alfa = alpharunz;

    // first two generations:  no mixing included
    // up squarks:

    double mstl2 = pow(mql1, 2) + (0.50 - 2.0 / 3.0 * sw2) * pow(mz, 2) * cos(2.0 * b_new);
    double mstr2 = pow(mur1, 2) + 2.0 / 3.0 * sw2 * pow(mz, 2) * cos(2.0 * b_new);

    double msu1 = sqrt(mstl2);
    double msu2 = sqrt(mstr2);

    // down squarks
    double msbl2 = pow(mql1, 2) + (-0.50 + 1.0 / 3.0 * sw2) * pow(mz, 2) * cos(2.0 * b_new);
    double msbr2 = pow(mdr1, 2) - 1.0 / 3.0 * sw2 * pow(mz, 2) * cos(2.0 * b_new);
    double msd1 = sqrt(msbl2);
    double msd2 = sqrt(msbr2);

    // sleptons
    double msel2 = pow(mel1, 2) + (-0.50 + sw2) * pow(mz, 2) * cos(2.0 * b_new);
    double mser2 = pow(mer1, 2) - sw2 * pow(mz, 2) * cos(2.0 * b_new);
    double msnl2 = pow(mel1, 2) + 0.50 * pow(mz, 2) * cos(2.0 * b_new);
    mse1 = sqrt(msel2);
    mse2 = sqrt(mser2);
    msn1 = sqrt(msnl2);

    // stop parameters
    mstl2 = pow(mql, 2) + (0.50 - 2.0 / 3 * sw2) * pow(mz, 2) * cos(2 * b_new);
    mstr2 = pow(mur, 2) + 2.0 / 3 * sw2 * pow(mz, 2) * cos(2 * b_new);
    double mlrt = at - mu / tb;
    thet = atan(2 * rmt * mlrt / (mstl2 - mstr2)) / 2;
    double ct = cos(thet);
    double st = sin(thet);
    mst1 = sqrt(pow(ct, 2) * (pow(rmt, 2) + mstl2) + pow(st, 2) * (pow(rmt, 2) + mstr2) + 2 * ct * st * rmt * mlrt);
    mst2 = sqrt(pow(st, 2) * (pow(rmt, 2) + mstl2) + pow(ct, 2) * (pow(rmt, 2) + mstr2) - 2 * ct * st * rmt * mlrt);
    mst1sfbp = mst1;
    mst2sfbp = mst2;

    //	printf("\n          .                   .                 %e %e %e %e",rmt,mlrt,mstl2,mstr2);

    if (isnan(mst1)) //!!added protection
    {
        mst1 = 92.0;
        if (irge == irgmax)
            sterr = -1.0;
    }

    if (isnan(mst2))
    {
        mst2 = 92.0;
        if (irge == irgmax)
            ;
        sterr = -1.0;
    }

    // sbottom parameters:

    msbl2 = pow(mql, 2) + (-0.50 + 1.0 / 3 * sw2) * pow(mz, 2) * cos(2 * b_new);
    msbr2 = pow(mdr, 2) - 1.0 / 3 * sw2 * pow(mz, 2) * cos(2 * b_new);
    double mlrb = ab - mu * tb;
    theb = atan(2 * rmb * mlrb / (msbl2 - msbr2)) / 2;
    double cb = cos(theb);
    double sb = sin(theb);
    double msb1 = sqrt(pow(cb, 2) * (pow(rmb, 2) + msbl2) + pow(sb, 2) * (pow(rmb, 2) + msbr2) + 2 * cb * sb * rmb * mlrb);
    double msb2 = sqrt(pow(sb, 2) * (pow(rmb, 2) + msbl2) + pow(cb, 2) * (pow(rmb, 2) + msbr2) - 2 * cb * sb * rmb * mlrb);
    msb1sfbp = msb1;
    msb2sfbp = msb2;

    //	printf("\n\n............................. %e %e\n\n",msb1sfbp,msb2sfbp);

    if (isnan(msb1))
    {
        msb1 = 92.0;
        if (irge == irgmax)
            sberr = -1.0;
    }

    if (isnan(msb2))
    {
        msb2 = 92.0;
        if (irge == irgmax)
            sberr = -1.0;
    }

    //
    //  stau parameters
    //

    msel2 = pow(mel, 2) + (-0.50 + sw2) * pow(mz, 2) * cos(2 * b_new);
    mser2 = pow(mer, 2) - sw2 * pow(mz, 2) * cos(2 * b_new);
    double msntau2 = pow(mel, 2) + 0.50 * pow(mz, 2) * cos(2 * b_new);
    double mlre = al - mu * tb;
    thel = atan(2 * rmtau * mlre / (msel2 - mser2)) / 2;
    cl = cos(thel);
    double sl = sin(thel);
    msta1 = sqrt(pow(cl, 2) * (pow(rmtau, 2) + msel2) + pow(sl, 2) * (pow(rmtau, 2) + mser2) + 2 * cl * sl * rmtau * mlre);
    msta2 = sqrt(pow(sl, 2) * (pow(rmtau, 2) + msel2) + pow(cl, 2) * (pow(rmtau, 2) + mser2) - 2 * cl * sl * rmtau * mlre);

    if (isnan(msta1))
    { // !!added protection
        msta1 = 92.0;
        if (irge == irgmax)
            stauerr = -1.0;
    }

    if (isnan(msta2))
    {
        msta2 = 92.0;
        if (irge == irgmax)
            stauerr = -1.0;
    }

    //   tau sneutrino:
    msntau = sqrt(msntau2);
    mz = mzsave;
    mw = mwsave;

    hml = mlrunz;
    hmh = mhrunz;
    hma = marunz;
    hmch = mchrunz;
    halfa = alpharunz; // SU_higgsrunz
}
