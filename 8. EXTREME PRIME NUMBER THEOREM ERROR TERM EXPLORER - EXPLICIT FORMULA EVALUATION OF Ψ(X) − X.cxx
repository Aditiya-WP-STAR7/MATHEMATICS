#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace std;

/*
===========================================================
 EXTREME PRIME NUMBER THEOREM ERROR TERM EXPLORER
 Explicit Formula Evaluation of ψ(x) − x

 Author: (Your Name)
 Target: MIT / IAS / ETH-level Portfolio Artifact
===========================================================
*/

// Type alias for numerical precision
using ld = long double;
using cplx = complex<ld>;

// ---------------------------------------------------------
// Non-trivial Riemann zeta zero container
// rho = 1/2 + i*gamma
// ---------------------------------------------------------
struct ZetaZero {
    ld gamma;
};

// ---------------------------------------------------------
// Compute x^rho / rho using complex arithmetic
// ---------------------------------------------------------
cplx x_to_rho_over_rho(ld x, const ZetaZero& zero) {
    cplx rho(0.5L, zero.gamma);
    cplx logx = log(x);
    cplx exponent = rho * logx;
    cplx x_rho = exp(exponent);
    return x_rho / rho;
}

// ---------------------------------------------------------
// Explicit formula for ψ(x)
// ---------------------------------------------------------
ld psi_explicit(ld x, const vector<ZetaZero>& zeros) {
    cplx sum = 0.0L;

    for (const auto& z : zeros) {
        cplx term = x_to_rho_over_rho(x, z);
        sum += term;
        sum += conj(term); // conjugate pair
    }

    ld correction1 = log(2.0L * M_PIl);
    ld correction2 = 0.5L * log(1.0L - pow(x, -2.0L));

    ld psi = x - real(sum) - correction1 - correction2;
    return psi;
}

// ---------------------------------------------------------
// Load a small demonstration set of zeta zeros
// (Replace with real datasets for serious research)
// ---------------------------------------------------------
vector<ZetaZero> load_demo_zeros() {
    vector<ZetaZero> zeros = {
        {14.1347251417347L},
        {21.0220396387716L},
        {25.0108575801457L},
        {30.4248761258595L},
        {32.9350615877392L},
        {37.5861781588257L},
        {40.9187190121475L}
    };
    return zeros;
}

// ---------------------------------------------------------
// Main CLI loop
// ---------------------------------------------------------
int main() {
    cout << fixed << setprecision(12);

    cout << "=============================================\n";
    cout << " EXTREME PRIME NUMBER ERROR TERM EXPLORER\n";
    cout << " ψ(x) − x via Explicit Formula\n";
    cout << "=============================================\n\n";

    vector<ZetaZero> zeros = load_demo_zeros();

    while (true) {
        ld x;
        cout << "Enter x (e.g., 1e20) or 0 to exit: ";
        cin >> x;

        if (!cin || x == 0.0L) {
            cout << "\nExiting program.\n";
            break;
        }

        if (x < 10.0L) {
            cout << "x must be large for asymptotic validity.\n\n";
            continue;
        }

        cout << "\nComputing ψ(x) using explicit formula...\n";
        cout << "Number of zeta zeros used: " << zeros.size() * 2 << "\n";

        ld psi_x = psi_explicit(x, zeros);
        ld delta = psi_x - x;

        cout << "\nResults:\n";
        cout << "ψ(x)      = " << psi_x << "\n";
        cout << "Δ(x)=ψ(x)-x = " << delta << "\n";
        cout << "\n---------------------------------------------\n\n";
    }

    return 0;
}
