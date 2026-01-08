#include <iostream>
#include <cmath>
#include <gmpxx.h>
#include <iomanip>

using namespace std;

/*
  Fourier expansion of the modular j-invariant:

  j(tau) = 1/q + 744 + sum_{n=1}^∞ c(n) q^n
  where q = exp(2πiτ)

  At τ = i√163:
  q = exp(-2π√163) ≈ 2.6 × 10⁻¹⁸

  This explains the near-integer phenomenon.
*/

// Known integer Fourier coefficients (first few)
const long long C[] = {
    196884,
    21493760,
    864299970,
    20245856256,
    333202640600
};

mpz_class compute_exact_core() {
    mpz_class base("640320");
    mpz_class cube = base * base * base;
    return -cube;
}

void compute_j_invariant(int terms) {
    cout << "\n[ Computing j(i√163) using Fourier expansion ]\n";

    long double sqrt163 = sqrt((long double)163.0);
    long double q = expl(-2.0L * M_PI * sqrt163);

    cout << fixed << setprecision(30);
    cout << "q = exp(-2π√163) = " << q << "\n\n";

    // Start with floating evaluation
    long double j_float = (1.0L / q) + 744.0L;

    long double qn = q;
    for (int i = 0; i < terms && i < 5; ++i) {
        j_float += C[i] * qn;
        qn *= q;
    }

    cout << "Floating approximation:\n";
    cout << "j(i√163) ≈ " << j_float << "\n\n";

    // Exact algebraic integer core
    mpz_class exact_core = compute_exact_core();
    mpz_class exact_value = exact_core + 744;

    cout << "Exact integer expression:\n";
    cout << "j(i√163) = -640320^3 + 744\n";
    cout << "            = " << exact_value << "\n\n";

    mpz_class difference = exact_value - mpz_class(j_float);
    cout << "Difference (Exact - Floating) ≈ " << difference << "\n";
}

int main() {
    cout << "===============================================\n";
    cout << "  Modular j-Invariant CLI (MIT-Level Program)\n";
    cout << "  Exact Arithmetic with GMP\n";
    cout << "===============================================\n";

    while (true) {
        int terms;
        cout << "\nEnter number of Fourier terms (1–5 recommended): ";
        cin >> terms;

        compute_j_invariant(terms);

        char choice;
        cout << "\nCompute again? (y/n): ";
        cin >> choice;

        if (choice != 'y' && choice != 'Y') {
            cout << "\nProgram terminated. Stay legendary.\n";
            break;
        }
    }

    return 0;
}
