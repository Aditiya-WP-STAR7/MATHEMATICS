#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

/*
    ============================================================
    MULTIDIMENSIONAL MONTE CARLO INTEGRATION (7D)
    Radiative Correction Approximation to Fine-Structure Constant
    ------------------------------------------------------------
    Author  : Aditiya Widodo Putra
    Purpose : High-order Feynman-like integral estimation
    Method  : Importance Sampling Monte Carlo
    ============================================================
*/

// Physical constants (SI units, normalized where appropriate)
constexpr double ALPHA_0 = 1.0 / 137.035999084; // CODATA
constexpr double PI = 3.141592653589793;
constexpr int DIM = 7;

// Probability density function for importance sampling
double gaussianPDF(double x, double sigma) {
    return exp(-x * x / (2 * sigma * sigma)) / (sqrt(2 * PI) * sigma);
}

// Feynman-like integrand (simplified but physically inspired)
double feynmanIntegrand(const vector<double>& k) {
    double sum_sq = 0.0;
    double interaction = 1.0;

    for (double v : k) {
        sum_sq += v * v;
        interaction *= cos(v);
    }

    // Mimics loop momentum suppression + interaction vertex
    return interaction * exp(-sum_sq);
}

// Monte Carlo integration with importance sampling
double monteCarloIntegral(
    long long samples,
    double& variance_out
) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dist(0.0, 1.0);

    double sum = 0.0;
    double sum_sq = 0.0;

    for (long long i = 0; i < samples; ++i) {
        vector<double> k(DIM);
        double weight = 1.0;

        for (int d = 0; d < DIM; ++d) {
            k[d] = dist(gen);
            weight *= gaussianPDF(k[d], 1.0);
        }

        double value = feynmanIntegrand(k) / weight;
        sum += value;
        sum_sq += value * value;
    }

    double mean = sum / samples;
    variance_out = (sum_sq / samples) - (mean * mean);
    return mean;
}

int main() {
    cout << fixed << setprecision(10);

    cout << "\n=============================================\n";
    cout << "  Quantum Monte Carlo: Fine-Structure Constant\n";
    cout << "=============================================\n";
    cout << "Dimension       : 7D\n";
    cout << "Method          : Importance Sampling\n";
    cout << "Base alpha      : " << ALPHA_0 << "\n\n";

    while (true) {
        long long samples;
        cout << "Enter number of Monte Carlo samples (e.g. 1e6): ";
        cin >> samples;

        if (!cin || samples <= 0) {
            cout << "Invalid input. Please enter a positive integer.\n";
            cin.clear();
            cin.ignore(10000, '\n');
            continue;
        }

        double variance = 0.0;
        double delta_alpha = monteCarloIntegral(samples, variance);

        double alpha_effective = ALPHA_0 + delta_alpha * 1e-4;
        double std_error = sqrt(variance / samples);

        cout << "\n===== RESULTS =====\n";
        cout << "Radiative Correction (Δα) : " << delta_alpha << "\n";
        cout << "Estimated α_eff          : " << alpha_effective << "\n";
        cout << "Standard Error           : ±" << std_error << "\n";
        cout << "Relative Precision       : "
             << fabs(std_error / delta_alpha) * 100.0 << " %\n";

        cout << "\nCompute another integral? (y/n): ";
        char choice;
        cin >> choice;

        if (choice != 'y' && choice != 'Y') {
            cout << "\nExiting program. Scientific computation complete.\n";
            break;
        }

        cout << "\n---------------------------------------------\n";
    }

    return 0;
}
