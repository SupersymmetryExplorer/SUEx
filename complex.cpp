#include <math.h>
#include <stdlib.h>
#include <complex>
#include <limits>

#include "complex.h"

const double pi = 3.14159265359; 		


complex mult(double x,complex y)
{
	complex ans;
	ans.re = x*y.re;
	ans.im = x*y.im;
	return(ans);
	}	

complex initiallize(double re,double im)
{
	complex ans;
	ans.re = re;
	ans.im = im;
	return ans;
	}		

double carg(complex x)
{
    if (x.re == 0.0 && x.im == 0.0)
        return 0.0;

    return atan2(x.im, x.re);
	}
	
complex adddouble(complex x, double y)
{
	x.re = x.re + y;
	return x;
	}

complex add(complex x, complex y)
{
	x.re = x.re + y.re;
	x.im = x.im + y.im;
	return x;
	}

complex clog(complex z)
{
    if (z.re == 0.0 && z.im == 0.0)
        return initiallize(-INFINITY, 0.0);

    const auto result = std::log(std::complex<long double>(z.re, z.im));
    return initiallize(static_cast<double>(result.real()),
                      static_cast<double>(result.imag()));
}

complex cmult(complex x, complex y)
{
	complex ans;
	ans.re = x.re*y.re - x.im*y.im;
	ans.im = x.re*y.im + x.im*y.re;
	return ans;
	}	
	
double cabs(complex x)
{
    return std::hypot(x.re, x.im);
}
	

	
complex csqrt(complex x)
{
    if (x.re == 0.0 && x.im == 0.0)
        return initiallize(0.0, 0.0);

    const auto result = std::sqrt(std::complex<long double>(x.re, x.im));
    return initiallize(static_cast<double>(result.real()),
                      static_cast<double>(result.imag()));
}

complex cdiv(double x, complex z)
{
    return cdivision(initiallize(x, 0.0), z);
}

complex cdivision(complex x, complex y)
{
    // Wider intermediates also protect division near the double range limits.
    const std::complex<long double> numerator(x.re, x.im);
    const std::complex<long double> denominator(y.re, y.im);
    const auto result = numerator / denominator;
    return initiallize(static_cast<double>(result.real()),
                      static_cast<double>(result.imag()));
}

complex cnegative(complex x)
{
	x.re *=-1.0;
	x.im *=-1.0;
	return x;
	}
		
double creal(complex z)
{
	return z.re;
	}

complex cpow(complex z, double x)
{
    // Preserve the existing zero-base conventions without evaluating log(0).
    if (x == 0.0)
        return initiallize(1.0, 0.0);
    if (z.re == 0.0 && z.im == 0.0) {
        if (x > 0.0)
            return initiallize(0.0, 0.0);
        if (x < 0.0)
            return initiallize(INFINITY, NAN);
    }

    const std::complex<long double> base(z.re, z.im);
    std::complex<long double> result;
    // Integer powers avoid introducing phase-rounding errors into exact zeros.
    if (std::isfinite(x) && x == std::trunc(x) &&
        x > std::numeric_limits<int>::min() &&
        x <= std::numeric_limits<int>::max()) {
        result = std::pow(base, static_cast<int>(x));
    } else {
        result = std::pow(base, static_cast<long double>(x));
    }
    return initiallize(static_cast<double>(result.real()),
                      static_cast<double>(result.imag()));
}
