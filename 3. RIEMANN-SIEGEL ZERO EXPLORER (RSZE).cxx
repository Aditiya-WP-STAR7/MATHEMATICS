#include <iostream>
#include <vector>
#include <cmath>
#include <mpfr.h>

using namespace std;

// Precision: 256 bits ≈ 77 decimal digits
static const int PREC = 256;

/*
 θ(t) = Im(log Γ(1/4 + it/2)) − (t/2) log π
 */
void theta(mpfr_t result, mpfr_t t) {
    mpfr_t a, b, lgamma_real, lgamma_imag, logpi, temp;
    mpfr_inits2(PREC, a, b, lgamma_real, lgamma_imag, logpi, temp, (mpfr_ptr) 0);

    // a = 1/4
    mpfr_set_d(a, 0.25, MPFR_RNDN);

    // b = t / 2
    mpfr_div_ui(b, t, 2, MPFR_RNDN);

    // Compute log Γ(a + i b)
    mpfr_lngamma_complex(lgamma_real, lgamma_imag, a, b, MPFR_RNDN);

    // temp = (t/2) * log(pi)
    mpfr_const_pi(logpi, MPFR_RNDN);
    mpfr_log(logpi, logpi, MPFR_RNDN);
    mpfr_mul(temp, b, logpi, MPFR_RNDN);

    // result = Im(log Γ) − temp
    mpfr_sub(result, lgamma_imag, temp, MPFR_RNDN);

    mpfr_clears(a, b, lgamma_real, lgamma_imag, logpi, temp, (mpfr_ptr) 0);
}

/*
 Riemann–Siegel Z(t)
 */
void Z_function(mpfr_t Z, mpfr_t t) {
    mpfr_t th, sum, n, limit, term, logn, phase, sqrt_n;
    mpfr_inits2(PREC, th, sum, n, limit, term, logn, phase, sqrt_n, (mpfr_ptr) 0);

    theta(th, t);
    mpfr_set_zero(sum, 1);

    // limit = floor(sqrt(t / (2π)))
    mpfr_const_pi(limit, MPFR_RNDN);
    mpfr_mul_ui(limit, limit, 2, MPFR_RNDN);
    mpfr_div(limit, t, limit, MPFR_RNDN);
    mpfr_sqrt(limit, limit, MPFR_RNDN);
    unsigned long N = mpfr_get_ui(limit, MPFR_RNDN);

    for (unsigned long i = 1; i <= N; ++i) {
        mpfr_set_ui(n, i, MPFR_RNDN);

        mpfr_log(logn, n, MPFR_RNDN);
        mpfr_mul(phase, t, logn, MPFR_RNDN);
        mpfr_sub(phase, th, phase, MPFR_RNDN);

        mpfr_cos(term, phase, MPFR_RNDN);
        mpfr_sqrt(sqrt_n, n, MPFR_RNDN);
        mpfr_div(term, term, sqrt_n, MPFR_RNDN);

        mpfr_add(sum, sum, term, MPFR_RNDN);
    }

    mpfr_mul_ui(Z, sum, 2, MPFR_RNDN);

    mpfr_clears(th, sum, n, limit, term, logn, phase, sqrt_n, (mpfr_ptr) 0);
}

int main() {
    cout << "=============================================\n";
    cout << " Riemann–Siegel Zeta Zero Explorer (RSZE)\n";
    cout << " High-Precision Analytic Number Theory Tool\n";
    cout << "=============================================\n\n";

    while (true) {
        cout << "Enter t (imaginary part of zero, approx) or -1 to exit:\n> ";

        double t_input;
        cin >> t_input;
        if (t_input < 0) break;

        mpfr_t t, Z;
        mpfr_inits2(PREC, t, Z, (mpfr_ptr) 0);
        mpfr_set_d(t, t_input, MPFR_RNDN);

        Z_function(Z, t);

        cout << "\nZ(t) = ";
        mpfr_out_str(stdout, 10, 40, Z, MPFR_RNDN);
        cout << "\n";

        if (mpfr_sgn(Z) == 0)
            cout << ">>> EXACT ZERO DETECTED <<<\n";
        else
            cout << ">>> Sign indicates proximity to a non-trivial zero <<<\n";

        cout << "\n---------------------------------------------\n\n";

        mpfr_clears(t, Z, (mpfr_ptr) 0);
    }

    cout << "Session ended. Mathematics never ends.\n";
    return 0;
}
