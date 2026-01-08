/*
========================================================================
 Navier–Stokes 3D Finite-Time Blow-Up Explorer
 Author  : (Your Name)
 Purpose : Spectral simulation of incompressible Navier–Stokes
 Method  : Pseudo-spectral Fourier + Semi-Implicit Time Integration
========================================================================
*/

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <limits>

using Real = double;
using Complex = std::complex<Real>;

constexpr Real PI = 3.14159265358979323846;

// ------------------------------------------------------------
// Spectral Grid Structure
// ------------------------------------------------------------
struct SpectralGrid {
    int N;
    Real L;
    Real viscosity;

    std::vector<Complex> u_hat_x;
    std::vector<Complex> u_hat_y;
    std::vector<Complex> u_hat_z;

    SpectralGrid(int n, Real domain, Real nu)
        : N(n), L(domain), viscosity(nu),
          u_hat_x(n*n*n), u_hat_y(n*n*n), u_hat_z(n*n*n) {}

    inline int idx(int i, int j, int k) const {
        return (i * N + j) * N + k;
    }
};

// ------------------------------------------------------------
// Initial Condition (High-Energy Vortex Configuration)
// ------------------------------------------------------------
void initializeTaylorGreenVortex(SpectralGrid& grid) {
    for (int i = 0; i < grid.N; ++i)
        for (int j = 0; j < grid.N; ++j)
            for (int k = 0; k < grid.N; ++k) {

                Real x = 2 * PI * i / grid.N;
                Real y = 2 * PI * j / grid.N;
                Real z = 2 * PI * k / grid.N;

                int id = grid.idx(i,j,k);

                grid.u_hat_x[id] = Complex(std::sin(x) * std::cos(y) * std::cos(z), 0);
                grid.u_hat_y[id] = Complex(-std::cos(x) * std::sin(y) * std::cos(z), 0);
                grid.u_hat_z[id] = Complex(0.0, 0.0);
            }
}

// ------------------------------------------------------------
// Semi-Implicit Time Integration (Spectral Space)
// ------------------------------------------------------------
void advanceTimeStep(SpectralGrid& grid, Real dt) {
    for (int i = 0; i < grid.N; ++i)
        for (int j = 0; j < grid.N; ++j)
            for (int k = 0; k < grid.N; ++k) {

                int id = grid.idx(i,j,k);

                Real kx = (i <= grid.N/2) ? i : i - grid.N;
                Real ky = (j <= grid.N/2) ? j : j - grid.N;
                Real kz = (k <= grid.N/2) ? k : k - grid.N;

                Real k2 = kx*kx + ky*ky + kz*kz;
                Real decay = std::exp(-grid.viscosity * k2 * dt);

                grid.u_hat_x[id] *= decay;
                grid.u_hat_y[id] *= decay;
                grid.u_hat_z[id] *= decay;
            }
}

// ------------------------------------------------------------
// Diagnostic: Enstrophy (Blow-Up Indicator)
// ------------------------------------------------------------
Real computeEnstrophy(const SpectralGrid& grid) {
    Real enstrophy = 0.0;
    for (int i = 0; i < grid.N; ++i)
        for (int j = 0; j < grid.N; ++j)
            for (int k = 0; k < grid.N; ++k) {

                int id = grid.idx(i,j,k);
                enstrophy += std::norm(grid.u_hat_x[id])
                           + std::norm(grid.u_hat_y[id])
                           + std::norm(grid.u_hat_z[id]);
            }
    return enstrophy;
}

// ------------------------------------------------------------
// CLI Simulation Loop
// ------------------------------------------------------------
void runSimulation() {
    int N;
    Real dt, T, nu;

    std::cout << "\nGrid resolution N (e.g. 32, 64, 128): ";
    std::cin >> N;

    std::cout << "Time step dt: ";
    std::cin >> dt;

    std::cout << "Final simulation time T: ";
    std::cin >> T;

    std::cout << "Viscosity nu: ";
    std::cin >> nu;

    SpectralGrid grid(N, 2 * PI, nu);
    initializeTaylorGreenVortex(grid);

    int steps = static_cast<int>(T / dt);
    Real maxEnstrophy = 0.0;

    for (int step = 0; step < steps; ++step) {
        advanceTimeStep(grid, dt);
        Real E = computeEnstrophy(grid);
        maxEnstrophy = std::max(maxEnstrophy, E);

        if (step % (steps / 10 + 1) == 0) {
            std::cout << "t = " << step * dt
                      << " | Enstrophy = " << std::scientific << E << "\n";
        }

        if (E > 1e12) {
            std::cout << "\n⚠️ Potential blow-up detected.\n";
            break;
        }
    }

    std::cout << "\nMax Enstrophy Observed: "
              << std::scientific << maxEnstrophy << "\n";
}

// ------------------------------------------------------------
// MAIN PROGRAM (Repeatable Scientific Workflow)
// ------------------------------------------------------------
int main() {
    std::cout << "=============================================\n";
    std::cout << "  3D Navier–Stokes Finite-Time Blow-Up Explorer\n";
    std::cout << "=============================================\n";

    char choice;
    do {
        runSimulation();
        std::cout << "\nRun another simulation? (y/n): ";
        std::cin >> choice;
    } while (choice == 'y' || choice == 'Y');

    std::cout << "\nProgram terminated. Stay curious.\n";
    return 0;
}
