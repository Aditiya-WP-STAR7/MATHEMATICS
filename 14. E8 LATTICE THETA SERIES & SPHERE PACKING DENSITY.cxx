#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace std;

/*
    E8 Lattice Theta Series Approximation
    ------------------------------------
    Θ_E8(q) = Σ_{v ∈ E8} q^{||v||^2}

    E8 definition:
    - Vectors in Z^8 with even coordinate sum
    - OR half-integer vectors with odd coordinate sum
*/

const int DIM = 8;

/* Compute squared Euclidean norm */
double normSquared(const vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += x * x;
    return s;
}

/* Generate integer lattice vectors */
void generateIntegerVectors(
    int idx,
    vector<double>& v,
    int bound,
    vector<vector<double>>& results
) {
    if (idx == DIM) {
        int sum = 0;
        for (double x : v) sum += static_cast<int>(x);
        if (sum % 2 == 0) results.push_back(v);
        return;
    }
    for (int k = -bound; k <= bound; ++k) {
        v[idx] = k;
        generateIntegerVectors(idx + 1, v, bound, results);
    }
}

/* Generate half-integer lattice vectors */
void generateHalfIntegerVectors(
    int idx,
    vector<double>& v,
    int bound,
    vector<vector<double>>& results
) {
    if (idx == DIM) {
        int sum = 0;
        for (double x : v) sum += static_cast<int>(2 * x);
        if (sum % 4 == 2) results.push_back(v);
        return;
    }
    for (int k = -bound; k <= bound; ++k) {
        v[idx] = k + 0.5;
        generateHalfIntegerVectors(idx + 1, v, bound, results);
    }
}

/* Theta series computation */
double thetaE8(double q, double maxNorm, int bound) {
    vector<vector<double>> lattice;
    vector<double> v(DIM);

    generateIntegerVectors(0, v, bound, lattice);
    generateHalfIntegerVectors(0, v, bound, lattice);

    double theta = 0.0;

    for (const auto& vec : lattice) {
        double nsq = normSquared(vec);
        if (nsq <= maxNorm)
            theta += pow(q, nsq);
    }

    return theta;
}

/* Exact E8 packing density */
double exactPackingDensity() {
    return pow(M_PI, 4) / 384.0;
}

int main() {
    cout << fixed << setprecision(10);

    while (true) {
        double q;
        double maxNorm;
        int bound;

        cout << "\n===== E8 Lattice Theta Series CLI =====\n";
        cout << "Enter q (0 < q < 1): ";
        cin >> q;

        cout << "Enter max norm squared cutoff: ";
        cin >> maxNorm;

        cout << "Enter coordinate bound (recommended 2 or 3): ";
        cin >> bound;

        if (q <= 0.0 || q >= 1.0) {
            cout << "Invalid q. Must be in (0,1).\n";
            continue;
        }

        cout << "\nComputing Θ_E8(q)... please wait...\n";

        double theta = thetaE8(q, maxNorm, bound);

        cout << "\n===== Results =====\n";
        cout << "Theta_E8(q) ≈ " << theta << endl;
        cout << "Exact E8 Packing Density = "
             << exactPackingDensity() << endl;

        cout << "\nCompute again? (y/n): ";
        char choice;
        cin >> choice;
        if (choice != 'y' && choice != 'Y')
            break;
    }

    cout << "\nProgram terminated. Stay legendary.\n";
    return 0;
}
