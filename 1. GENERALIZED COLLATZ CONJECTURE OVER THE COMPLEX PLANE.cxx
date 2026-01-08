#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <iomanip>
#include <limits>

using namespace std;

// ==========================
// Generalized Collatz on ℂ
// ==========================

// Determine "parity" via modulus-based criterion
bool isEvenComplex(const complex<double>& z) {
    long long m = static_cast<long long>(floor(abs(z)));
    return (m % 2 == 0);
}

// Collatz iteration in complex plane
complex<double> collatzComplex(const complex<double>& z) {
    if (isEvenComplex(z)) {
        return z / 2.0;
    } else {
        return 3.0 * z + complex<double>(1.0, 0.0);
    }
}

// Orbit simulation
void simulateOrbit(
    complex<double> z0,
    int maxIterations,
    double divergenceThreshold
) {
    complex<double> z = z0;

    cout << "\nIteration | Re(z)               Im(z)               | |z|\n";
    cout << "------------------------------------------------------------------\n";

    for (int i = 0; i < maxIterations; ++i) {
        double modulus = abs(z);

        cout << setw(9) << i << " | "
             << setw(18) << fixed << setprecision(10) << z.real() << " "
             << setw(18) << fixed << setprecision(10) << z.imag() << " | "
             << setw(12) << modulus << "\n";

        if (modulus > divergenceThreshold) {
            cout << "\n⚠ Orbit diverged beyond threshold.\n";
            return;
        }

        z = collatzComplex(z);
    }

    cout << "\n✓ Simulation completed without divergence.\n";
}

// ==========================
// Main CLI Program
// ==========================
int main() {
    cout << "====================================================\n";
    cout << " Generalized Collatz Conjecture over ℂ (CLI Program)\n";
    cout << " Research-Grade Numerical Experiment\n";
    cout << "====================================================\n";

    while (true) {
        double realPart, imagPart;
        int iterations;
        double threshold;

        cout << "\nEnter real part of initial z: ";
        cin >> realPart;

        cout << "Enter imaginary part of initial z: ";
        cin >> imagPart;

        cout << "Enter maximum iterations: ";
        cin >> iterations;

        cout << "Enter divergence threshold (e.g. 1e6): ";
        cin >> threshold;

        complex<double> z0(realPart, imagPart);

        simulateOrbit(z0, iterations, threshold);

        char choice;
        cout << "\nRun another simulation? (y/n): ";
        cin >> choice;

        if (choice != 'y' && choice != 'Y') {
            cout << "\nExiting program. Thank you for exploring ℂ.\n";
            break;
        }
    }

    return 0;
}
