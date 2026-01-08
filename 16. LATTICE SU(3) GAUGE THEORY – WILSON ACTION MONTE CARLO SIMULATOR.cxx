#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <cmath>
#include <iomanip>

using namespace std;

/*
===========================================================
 LATTICE SU(3) GAUGE THEORY
 Monte Carlo Simulation with Wilson Action
 Author: (Your Name)
 Purpose: Demonstrate professional-level lattice QCD numerics
===========================================================
*/

static constexpr int Nc = 3;
using Complex = complex<double>;

/* ================= SU(3) MATRIX ================= */

struct SU3 {
    Complex m[Nc][Nc];

    static SU3 identity() {
        SU3 U;
        for (int i = 0; i < Nc; ++i)
            for (int j = 0; j < Nc; ++j)
                U.m[i][j] = (i == j) ? 1.0 : 0.0;
        return U;
    }
};

/* ================= MATRIX OPERATIONS ================= */

SU3 dagger(const SU3& U) {
    SU3 R;
    for (int i = 0; i < Nc; ++i)
        for (int j = 0; j < Nc; ++j)
            R.m[i][j] = conj(U.m[j][i]);
    return R;
}

SU3 operator*(const SU3& A, const SU3& B) {
    SU3 R;
    for (int i = 0; i < Nc; ++i)
        for (int j = 0; j < Nc; ++j) {
            R.m[i][j] = 0.0;
            for (int k = 0; k < Nc; ++k)
                R.m[i][j] += A.m[i][k] * B.m[k][j];
        }
    return R;
}

double real_trace(const SU3& U) {
    double tr = 0.0;
    for (int i = 0; i < Nc; ++i)
        tr += real(U.m[i][i]);
    return tr;
}

/* ================= RANDOM SU(3) GENERATOR ================= */

SU3 random_su3(mt19937& rng, double eps) {
    normal_distribution<double> N(0.0, eps);

    SU3 U = SU3::identity();

    for (int i = 0; i < Nc; ++i)
        for (int j = 0; j < Nc; ++j)
            U.m[i][j] += Complex(N(rng), N(rng));

    return U;
}

/* ================= LATTICE ================= */

struct Lattice {
    int L;
    vector<SU3> link;

    Lattice(int L_) : L(L_), link(4 * pow(L, 4), SU3::identity()) {}

    int index(int x, int y, int z, int t, int mu) const {
        x = (x + L) % L;
        y = (y + L) % L;
        z = (z + L) % L;
        t = (t + L) % L;
        return (((t * L + z) * L + y) * L + x) * 4 + mu;
    }

    SU3& U(int x, int y, int z, int t, int mu) {
        return link[index(x, y, z, t, mu)];
    }
};

/* ================= WILSON ACTION ================= */

double wilson_action(const Lattice& lat, double beta) {
    double S = 0.0;
    int L = lat.L;

    for (int x = 0; x < L; ++x)
        for (int y = 0; y < L; ++y)
            for (int z = 0; z < L; ++z)
                for (int t = 0; t < L; ++t)
                    for (int mu = 0; mu < 4; ++mu)
                        for (int nu = mu + 1; nu < 4; ++nu) {

                            SU3 U1 = lat.link[lat.index(x, y, z, t, mu)];
                            SU3 U2 = lat.link[lat.index(x + (mu == 0), y + (mu == 1), z + (mu == 2), t + (mu == 3), nu)];
                            SU3 U3 = dagger(lat.link[lat.index(x + (nu == 0), y + (nu == 1), z + (nu == 2), t + (nu == 3), mu)]);
                            SU3 U4 = dagger(lat.link[lat.index(x, y, z, t, nu)]);

                            SU3 P = U1 * U2 * U3 * U4;
                            S += 1.0 - real_trace(P) / Nc;
                        }

    return beta * S;
}

/* ================= MONTE CARLO ================= */

void metropolis(Lattice& lat, double beta, int sweeps) {
    mt19937 rng(random_device{}());
    uniform_real_distribution<double> U(0.0, 1.0);

    double S = wilson_action(lat, beta);

    for (int sweep = 0; sweep < sweeps; ++sweep) {
        for (auto& link : lat.link) {
            SU3 old = link;
            SU3 proposal = random_su3(rng, 0.05);
            link = proposal * old;

            double S_new = wilson_action(lat, beta);
            double dS = S_new - S;

            if (dS < 0 || U(rng) < exp(-dS)) {
                S = S_new;
            } else {
                link = old;
            }
        }

        cout << "Sweep " << sweep + 1 << " | Action = " << fixed << setprecision(6) << S << endl;
    }
}

/* ================= MAIN LOOP ================= */

int main() {
    cout << "\n=== SU(3) Lattice Gauge Theory Simulator ===\n";

    while (true) {
        int L, sweeps;
        double beta;

        cout << "\nEnter lattice size L (e.g. 4): ";
        cin >> L;

        cout << "Enter beta (e.g. 5.7): ";
        cin >> beta;

        cout << "Enter Monte Carlo sweeps: ";
        cin >> sweeps;

        Lattice lat(L);
        metropolis(lat, beta, sweeps);

        char again;
        cout << "\nRun another simulation? (y/n): ";
        cin >> again;
        if (again != 'y' && again != 'Y') break;
    }

    cout << "\nSimulation finished. Exiting professionally.\n";
    return 0;
}
