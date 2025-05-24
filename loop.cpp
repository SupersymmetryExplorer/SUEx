#include <stdio.h>
#include "complex.h"
#include "loop.h"
#include "variables.h"

// This header file is for calculating the loop integrals using Passarino-Veltman integrals and tensor reduction.
//
//
// The algorithm is discussed in :
//
// 1) http://www.hephy.at/fileadmin/user_upload/Vortraege/Svit-lecture.pdf
//
// 2) One-loop corrections for e+e- annihilation into l mu+mu-  in the weinberg model
//    G. Passarino and  M. Veltman
//    Nuclear Physics B160  (1979)  151-207
//
// 3) Gauge theories of the strong and electroweak interaction
//		By Manfred Bohm, Manfred Bohm (Dr. rer. nat.), Ansgar Denner, Hans Joos
//
// 4) Techniques for the calculation of electroweak radiative corrections at the one-loop level and results for W-physics at LEP200
//    By A. Denner
//		arXiv:0709.1075v1  [hep-ph]         ... (All the calculations are done here)
//

// 5) Precision corrections in the minimal supersymmetric standard model
//		By D.M. Pierce, Jonathan A. Bagger, K.T. Matchev and Ren-jie Zhang
//		arXiv:hep-ph/9606211v3 					... (All these equations are taken from this paper)

//
// The passarino-veltman one (A) and two points (B0,B1) functions
//  needed for the evalution of the radiative corrections (and also v_loop).
// --------------------------------------------------------------------------------------------------------------------------------------

double su_a(double m)
{

    if (m != 0.0)
        return (m * m * (1.0 - log(m * m / scale / scale)));
    else
        return 0.0;
}

//
// The passarino-veltman one (A) and two points (B0,B1)
//	Functions needed for the evalution of the radiative corrections (and also v_loop).
// The following auxulary functions are required for calculating B0 and B1
// -----------------------------------------------------------------------------------------------------------------------------------
complex xlogx(complex x)
{
    if (cabs(x) == 0)
        return initiallize(0.0, 0.0);
    else
        return cmult(x, clog(x));
}

//
// auxiliary functions used by the B0,B1 two-point functions
// -----------------------------------------------------------------------------------------------------------------------------------

void roots(double p, double m1, double m2, complex *x1, complex *x2, complex *y1, complex *y2, complex *r)
{
    static double eps = 1e-20;
    complex ieps = initiallize(0.0, eps);
    complex onepeps = adddouble(ieps, 1.0);
    complex onemeps = adddouble(mult(-1.0, ieps), 1.0);

    double q;
    complex discriminent;

    discriminent = initiallize(p * (p - 2 * (m1 + m2)) + pow((m1 - m2), 2), 0.0);

    *r = csqrt(discriminent);

    q = p + m1 - m2;

    *x1 = mult(1.0 / 2.0 / p, adddouble(*r, q));
    *x2 = mult(1.0 / 2.0 / p, adddouble(cnegative(*r), q));

    if (cabs(*x2) > cabs(*x1))
    {
        *x1 = cdiv(m1 / p, *x2);
    }
    else if (cabs(*x1) > cabs(*x2))
    {

        *x2 = cdiv(m1 / p, *x1);
    }

    *x1 = add(*x1, mult(cabs(mult(p, *x1)) / p, ieps));
    *x2 = add(*x2, mult(-cabs(mult(p, *x2)) / p, ieps));

    q = p - m1 + m2;
    *y2 = mult(1.0 / 2.0 / p, adddouble(*r, q));
    *y1 = mult(1.0 / 2.0 / p, adddouble(mult(-1, *r), q));

    if (cabs(*y2) > cabs(*y1))
        *y1 = cdiv(m2 / p, *y2);
    else if (cabs(*y1) > cabs(*y2))
        *y2 = cdiv(m2 / p, *y1);

    *y1 = add(*y1, mult(-cabs(mult(p, *y1)) / p, ieps));
    *y2 = add(*y2, mult(cabs(mult(p, *y2)) / p, ieps));
}

//
// Auxiliary functions used by the B0,B1 two-point functions
// -----------------------------------------------------------------------------------------------------------------------------------

complex fpv(int n, complex x, complex y)
{
    static double acc = 1e-12;
    static double eps = 1e-20;
    complex ieps = initiallize(0.0, eps);
    complex onepeps = adddouble(ieps, 1.0);
    complex onemeps = adddouble(mult(-1.0, ieps), 1.0);

    complex fpvv, xm;

    if (cabs(x) < 10)
    {
        if (n == 0)
            fpvv = mult(-1.0, clog(cmult(mult(-1.0, y), cdiv(1.0, x))));
        else if (cabs(x) < acc)
            fpvv = initiallize(-1.0 / n, 0.0);
        else
        {
            fpvv = initiallize(0.0, 0.0);
            xm = initiallize(1.0, 0.0);
            for (int m = 0; m <= n - 1; m++)
            {
                fpvv = add(fpvv, mult(-1.0 / (n - m), xm));
                xm = cmult(xm, x);
            }

            fpvv = add(fpvv, cnegative(cmult(xm, clog(cdivision(cnegative(y), x)))));
        }
    }
    else
    {
        fpvv = initiallize(0.0, 0.0);
        xm = initiallize(1.0, 0.0);
        for (int m = 1; m <= 30; m++)
        {
            xm = cmult(xm, cdiv(1.0, x));
            fpvv = add(fpvv, mult(1.0 / (m + n), xm));
            if (cabs(cmult(xm, cdiv(1.0, fpvv))) < pow(acc, 2))
                return fpvv;
        }
    }

    return fpvv;
}

//
// The main scalar two-point function B0
//	Ref: arXiv:0709.1075v1  [hep-ph]
// -----------------------------------------------------------------------------------------------------------------------------------

complex su_b0(double p, double mm1, double mm2)
{
    static double acc = 1e-12;
    static double eps = 1e-20;
    complex ieps = initiallize(0.0, eps);
    complex onepeps = adddouble(ieps, 1.0);
    complex onemeps = adddouble(mult(-1.0, ieps), 1.0);

    complex x1, x2, y1, y2, r;
    double minacc;
    complex b0;

    double m1 = pow(mm1, 2);
    double m2 = pow(mm2, 2);
    double divergence = 0.0;
    double mudim2 = pow(scale, 2);
    minacc = acc * (m1 + m2);

    // General case
    if (fabs(p) > minacc)
    {
        roots(p, m1, m2, &x1, &x2, &y1, &y2, &r);

        if (cabs(y1) > .50 && cabs(y2) > .50)
        {
            b0 = add(initiallize(-log(m2 / mudim2), 0.0), mult(-1.0, add(fpv(1, x1, y1), fpv(1, x2, y2))));
        }
        else if (cabs(x1) < 10 && cabs(x2) < 10)
        {
            complex b01 = add(xlogx(mult(-1.0, x1)), xlogx(mult(-1.0, x2)));
            complex b02 = add(xlogx(y1), xlogx(y2));
            b0 = add(add(initiallize(2, 0.0), mult(-1.0, clog(mult(p / mudim2, onemeps)))), add(b01, mult(-1.0, b02)));
        }
        else if (cabs(x1) > .50 && cabs(x2) > .50)
        {
            b0 = add(initiallize(-log(m1 / mudim2), 0.0), mult(-1, add(fpv(1, y1, x1), fpv(1, y2, x2))));
        }
        else
        {
            b0 = initiallize(999.0e300, 0.0);
        }
    }

    // Zero momentum
    else if (fabs(m1 - m2) > minacc)
    {
        x2 = mult(m1 / (m1 - m2), onemeps);
        y2 = mult(m2 / (m2 - m1), onemeps);
        if (cabs(y2) > .50)
            b0 = add(initiallize(-log(m2 / mudim2), 0.0), mult(-1, fpv(1, x2, y2)));
        else
            b0 = add(initiallize(-log(m1 / mudim2), 0.0), mult(-1, fpv(1, y2, x2)));
    }

    else
        b0 = initiallize(-log(m2 / mudim2), 0.0);

    return adddouble(b0, divergence);
}

//
// One-loop corrections for e+e- annihilation into l mu+mu-  in the weinberg model
//    G. Passarino and  M. Veltman
//    Nuclear Physics B160  (1979)  151-207
// -----------------------------------------------------------------------------------------------------------------------------------

complex su_b1(double s, double mi, double mj)
{
    complex su_b1v;

    if (mi == mj)
        su_b1v = mult(1.0 / 2, su_b0(s, mi, mj));
    else
    {

        complex sxz;
        sxz = su_b0(s, mi, mj);

        su_b1v = mult(1.0 / 2, add(initiallize(su_a(mj) / s - su_a(mi) / s, 0.0), mult((1.0 + pow(mi, 2) / s - pow(mj, 2) / s), su_b0(s, mi, mj)))); //.. eq(D6)
    }
    return su_b1v;
}

//
// Precision corrections in the minimal supersymmetric standard model
//		By D.M. Pierce, Jonathan A. Bagger, K.T. Matchev and Ren-jie Zhang      Eq(B10), Eq(B14)
//		arXiv:9606/9606211v3  [hep-ph]
// -----------------------------------------------------------------------------------------------------------------------------------

complex su_bt22(double s, double mi, double mj)
{
    complex su_bt22v;
    double ans;

    su_bt22v = adddouble(mult((pow(mi, 2) / s + pow(mj, 2) / s - 1.0 / 2), su_b0(s, mi, mj)), su_a(mi) / s / 2 + su_a(mj) / s / 2);
    ans = (pow(mj, 2) / s - pow(mi, 2) / s) / 2;

    su_bt22v = add(su_bt22v, mult(-ans, adddouble(mult((pow(mj, 2) / s - pow(mi, 2) / s), su_b0(s, mi, mj)), -su_a(mj) / s + su_a(mi) / s)));
    su_bt22v = adddouble(su_bt22v, +pow(mi, 2) / s + pow(mj, 2) / s - 1.0 / 3 - 3 * su_a(mi) / s / 2 - 3 * su_a(mj) / s / 2);
    su_bt22v = mult(s / 6.0, su_bt22v);

    return su_bt22v;
}

//
// Precision corrections in the minimal supersymmetric standard model
//		By D.M. Pierce, Jonathan A. Bagger, K.T. Matchev and Ren-jie Zhang      Eq(B13)
//		arXiv:9606/9606211v3  [hep-ph]
// -----------------------------------------------------------------------------------------------------------------------------------

complex su_bh(double s, double mi, double mj)
{
    complex temp, su_bhv;

    temp = initiallize(0.0, 0.0);
    temp = adddouble(temp, su_a(mi) / s / 2 + su_a(mj) / s / 2);

    temp = add(temp, mult((pow(mi, 2) / s + pow(mj, 2) / s - 1.0 / 2), su_b0(s, mi, mj)));
    temp = add(temp, mult((pow(mj, 2) - pow(mi, 2)) / s / 2, adddouble(mult((-pow(mj, 2) + pow(mi, 2)) / s, su_b0(s, mi, mj)), su_a(mj) / s - su_a(mi) / s)));
    temp = adddouble(temp, +pow(mi, 2) / s + pow(mj, 2) / s - 1.0 / 3);
    su_bhv = mult(4.0 * s / 6, temp);
    su_bhv = add(su_bhv, mult(s, adddouble(mult((1.0 - pow(mi, 2) / s - pow(mj, 2) / s), su_b0(s, mi, mj)), -su_a(mi) / s - su_a(mj) / s)));
    return su_bhv;
}

//
// Precision corrections in the minimal supersymmetric standard model
//		By D.M. Pierce, Jonathan A. Bagger, K.T. Matchev and Ren-jie Zhang      Eq(B11)
//		arXiv:9606/9606211v3  [hep-ph]
// -----------------------------------------------------------------------------------------------------------------------------------

complex su_bf(double s, double mi, double mj)
{
    complex su_bfv = mult(s, adddouble(mult(-1 * (2.0 + 2 * pow(mi, 2) / s - pow(mj, 2) / s), su_b0(s, mi, mj)), su_a(mi) / s - 2 * su_a(mj) / s));
    return (su_bfv);
}

//
// Precision corrections in the minimal supersymmetric standard model
//		By D.M. Pierce, Jonathan A. Bagger, K.T. Matchev and Ren-jie Zhang      Eq(B12)
//		arXiv:9606/9606211v3  [hep-ph]
// -----------------------------------------------------------------------------------------------------------------------------------

complex su_bg(double s, double mi, double mj)
{
    complex su_bgv;
    su_bgv = mult(s, adddouble(mult((1.0 - pow(mi, 2) / s - pow(mj, 2) / s), su_b0(s, mi, mj)), -su_a(mi) / s - su_a(mj) / s));
    return su_bgv;
}

//
// Calculation for Higgs.
// -----------------------------------------------------------------------------------------------------------------------------------

complex f0_hdec(double m1, double m2, double qsq)
{
    complex cd, cr, cq2, ieps, cbet, cxx;
    double m1sq = m1 * m1;
    double m2sq = m2 * m2;

    ieps = initiallize(1.0, 1.0e-12);
    cq2 = mult(qsq, ieps);
    cd = cdiv((m1sq - m2sq), cq2);
    cr = csqrt(add(cpow(adddouble(cd, 1.0), 2.0), cdiv(-4.0 * m1sq, cq2)));

    complex ta;
    double tb;
    if (qsq == 0.0)
        return initiallize(0.0, 0.0);
    else
    {

        if (fabs(m1 - m2) < 1.0e-5)
        {
            return adddouble(cmult(cr, clog(cdivision(adddouble(cr, 1.0), adddouble(cr, -1.0)))), -2.0);
        }
        else
        {
            cbet = csqrt(adddouble(cnegative(cdiv(4 * m1 * m2, adddouble(cq2, -(m1 - m2) * (m1 - m2)))), 1.0));
            cxx = cdivision(adddouble(cbet, -1.0), adddouble(cbet, +1.0));
            ta = mult((qsq - (m1 - m2) * (m1 - m2)) / qsq, cmult(cbet, clog(cxx)));
            tb = log(m2sq / m1sq) * (((qsq + m2sq - m1sq) / 2.0 / qsq) - m2sq / (m2sq - m1sq));

            return cnegative(adddouble(ta, 1.0 - tb));
        }
    }
    return initiallize(0.0, 0.0);
}
