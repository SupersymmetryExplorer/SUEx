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

// Stable massive-loop limits and corrections for |s| << (|mi| + |mj|)^2.
static bool loop_small_momentum(double s, double mi, double mj)
{
    const long double mass_sum = fabsl(mi) + fabsl(mj);
    return mass_sum > 0.0L &&
           fabsl(s) <= 1.0e-3L * mass_sum * mass_sum;
}

static long double loop_b0_zero(long double a, long double b)
{
    const long double log_scale = 2.0L * logl(fabsl(scale));
    if (a == 0.0L)
        return 1.0L - logl(b) + log_scale;
    if (b == 0.0L)
        return 1.0L - logl(a) + log_scale;

    const long double sum = a + b;
    const long double q = (a - b) / sum;
    if (fabsl(q) < 0.25L) {
        long double power = q * q;
        long double correction = 0.0L;
        for (int n = 1; n <= 16; ++n) {
            correction += power / ((2.0L * n) * (2.0L * n + 1.0L));
            power *= q * q;
        }
        return -logl(sum / 2.0L) + log_scale + correction;
    }
    return 1.0L + log_scale - (a * logl(a) - b * logl(b)) / (a - b);
}

static long double loop_b0_slope_zero(long double a, long double b)
{
    if (a == 0.0L)
        return 0.5L / b;
    if (b == 0.0L)
        return 0.5L / a;

    const long double sum = a + b;
    const long double difference = a - b;
    const long double q = difference / sum;
    if (fabsl(q) < 0.25L) {
        long double power = 1.0L;
        long double result = 0.0L;
        for (int n = 0; n <= 16; ++n) {
            result += power / ((2.0L * n + 1.0L) * (2.0L * n + 3.0L));
            power *= q * q;
        }
        return result / sum;
    }
    return sum / (2.0L * difference * difference)
         - a * b * (logl(a) - logl(b)) /
           (difference * difference * difference);
}

template <typename Function>
static long double loop_parameter_integral(const Function& function)
{
    const long double fl = function(0.0L);
    const long double fm = function(0.5L);
    const long double fr = function(1.0L);
    const long double fq1 = function(0.25L);
    const long double fq3 = function(0.75L);
    const long double magnitude = fmaxl(fmaxl(fabsl(fl), fabsl(fr)),
        fmaxl(fabsl(fm), fmaxl(fabsl(fq1), fabsl(fq3))));
    if (magnitude == 0.0L)
        return 0.0L;

    auto integrate = [&](auto&& self, long double left, long double right,
                         long double left_value, long double mid_value,
                         long double right_value, long double whole,
                         long double tolerance, int depth) -> long double {
        const long double middle = (left + right) / 2.0L;
        const long double f1 = function((left + middle) / 2.0L);
        const long double f2 = function((middle + right) / 2.0L);
        const long double first = (middle - left) *
            (left_value + 4.0L * f1 + mid_value) / 6.0L;
        const long double second = (right - middle) *
            (mid_value + 4.0L * f2 + right_value) / 6.0L;
        const long double error = first + second - whole;
        if (fabsl(error) <= 15.0L * tolerance)
            return first + second + error / 15.0L;
        if (depth == 0)
            return NAN;
        return self(self, left, middle, left_value, f1, mid_value,
                    first, tolerance / 2.0L, depth - 1)
             + self(self, middle, right, mid_value, f2, right_value,
                    second, tolerance / 2.0L, depth - 1);
    };
    return integrate(integrate, 0.0L, 1.0L, fl, fm, fr,
                     (fl + 4.0L * fm + fr) / 6.0L,
                     2.0e-13L * magnitude, 48);
}

static long double loop_momentum_ratio(long double t, long double a,
                                       long double b, long double s)
{
    // Removable endpoints for a single massless internal particle.
    if (t == 0.0L && b == 0.0L)
        return s / a;
    if (t == 1.0L && a == 0.0L)
        return s / b;
    return s * t * (1.0L - t) / (t * a + (1.0L - t) * b);
}

static complex loop_stable_b0(double s, double mi, double mj)
{
    if (!loop_small_momentum(s, mi, mj))
        return su_b0(s, mi, mj);
    const long double a = (long double)mi * mi;
    const long double b = (long double)mj * mj;
    long double value = loop_b0_zero(a, b);
    if (s != 0.0) {
        value += loop_parameter_integral([&](long double t) {
            return -log1pl(-loop_momentum_ratio(t, a, b, s));
        });
    }
    return initiallize(static_cast<double>(value), 0.0);
}

static complex loop_small_b1(double s, double mi, double mj)
{
    const long double a = (long double)mi * mi;
    const long double b = (long double)mj * mj;
    long double value = (loop_b0_zero(a, b) +
                         (a - b) * loop_b0_slope_zero(a, b)) / 2.0L;
    if (s != 0.0) {
        value += loop_parameter_integral([&](long double t) {
            return -(1.0L - t) *
                log1pl(-loop_momentum_ratio(t, a, b, s));
        });
    }
    return initiallize(static_cast<double>(value), 0.0);
}

static complex loop_small_bt22(double s, double mi, double mj)
{
    const long double a = (long double)mi * mi;
    const long double b = (long double)mj * mj;
    const long double difference = a - b;
    long double value = difference * difference *
                         loop_b0_slope_zero(a, b) / 4.0L;
    if (s != 0.0) {
        const long double log_scale = 2.0L * logl(fabsl(scale));
        value += loop_parameter_integral([&](long double t) {
            const long double mass = t * a + (1.0L - t) * b;
            if (mass == 0.0L)
                return 0.0L;
            const long double u = loop_momentum_ratio(t, a, b, s);
            // u + (1-u)*log(1-u), expanded without cancellation.
            const long double remainder = u * u * (0.5L + u *
                (1.0L / 6.0L + u * (1.0L / 12.0L + u *
                (1.0L / 20.0L + u * (1.0L / 30.0L + u *
                (1.0L / 42.0L + u / 56.0L))))));
            return 0.5L * (s * t * (1.0L - t) *
                (logl(mass) - log_scale) - mass * remainder);
        });
    }
    return initiallize(static_cast<double>(value), 0.0);
}

complex su_b1(double s, double mi, double mj)
{
    if (loop_small_momentum(s, mi, mj))
        return loop_small_b1(s, mi, mj);
    if (mi == mj)
        return mult(0.5, su_b0(s, mi, mj));
    return mult(0.5 / s,
        adddouble(mult(s + mi * mi - mj * mj, su_b0(s, mi, mj)),
                  su_a(mj) - su_a(mi)));
}

//
// Precision corrections in the minimal supersymmetric standard model
//		By D.M. Pierce, Jonathan A. Bagger, K.T. Matchev and Ren-jie Zhang      Eq(B10), Eq(B14)
//		arXiv:9606/9606211v3  [hep-ph]
// -----------------------------------------------------------------------------------------------------------------------------------

complex su_bt22(double s, double mi, double mj)
{
    if (s == 0.0 && mi == 0.0 && mj == 0.0)
        return initiallize(0.0, 0.0);
    if (loop_small_momentum(s, mi, mj))
        return loop_small_bt22(s, mi, mj);
    const double a = mi * mi;
    const double b = mj * mj;
    const double difference = a - b;
    const complex b0 = su_b0(s, mi, mj);
    const complex divided_term = mult(-difference / (2.0 * s),
        adddouble(mult(difference, b0), -su_a(mi) + su_a(mj)));
    return mult(1.0 / 6.0, adddouble(
        add(mult(a + b - s / 2.0, b0), divided_term),
        a + b - s / 3.0 - su_a(mi) - su_a(mj)));
}

//
// Precision corrections in the minimal supersymmetric standard model
//		By D.M. Pierce, Jonathan A. Bagger, K.T. Matchev and Ren-jie Zhang      Eq(B13)
//		arXiv:9606/9606211v3  [hep-ph]
// -----------------------------------------------------------------------------------------------------------------------------------

complex su_bh(double s, double mi, double mj)
{
    if (s == 0.0 && mi == 0.0 && mj == 0.0)
        return initiallize(0.0, 0.0);
    return add(mult(4.0, su_bt22(s, mi, mj)),
        mult(s - mi * mi - mj * mj, loop_stable_b0(s, mi, mj)));
}

//
// Precision corrections in the minimal supersymmetric standard model
//		By D.M. Pierce, Jonathan A. Bagger, K.T. Matchev and Ren-jie Zhang      Eq(B11)
//		arXiv:9606/9606211v3  [hep-ph]
// -----------------------------------------------------------------------------------------------------------------------------------

complex su_bf(double s, double mi, double mj)
{
    if (s == 0.0 && mi == 0.0 && mj == 0.0)
        return initiallize(0.0, 0.0);
    return adddouble(
        mult(-(2.0 * s + 2.0 * mi * mi - mj * mj),
             loop_stable_b0(s, mi, mj)),
        su_a(mi) - 2.0 * su_a(mj));
}

//
// Precision corrections in the minimal supersymmetric standard model
//		By D.M. Pierce, Jonathan A. Bagger, K.T. Matchev and Ren-jie Zhang      Eq(B12)
//		arXiv:9606/9606211v3  [hep-ph]
// -----------------------------------------------------------------------------------------------------------------------------------

complex su_bg(double s, double mi, double mj)
{
    if (s == 0.0 && mi == 0.0 && mj == 0.0)
        return initiallize(0.0, 0.0);
    return adddouble(
        mult(s - mi * mi - mj * mj, loop_stable_b0(s, mi, mj)),
        -su_a(mi) - su_a(mj));
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
