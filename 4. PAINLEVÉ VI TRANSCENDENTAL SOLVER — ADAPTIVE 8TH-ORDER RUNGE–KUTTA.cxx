#include <iostream>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace std;

/*
=========================================================
 Painlevé VI Transcendental Solver
 Adaptive Runge–Kutta 8th Order (Embedded Error Control)
 Author Style: Research-Grade Numerical Analysis
=========================================================
*/

struct State {
    double y;
    double v;
};

const double SAFETY = 0.9;
const double MIN_STEP = 1e-10;
const double MAX_STEP = 0.1;

/* Painlevé VI system */
State painleveVI(double t, const State& s) {
    State ds;

    double y = s.y;
    double v = s.v;

    if (fabs(y) < 1e-12 || fabs(y - 1.0) < 1e-12 || fabs(y - t) < 1e-12) {
        ds.y = 0.0;
        ds.v = 0.0;
        return ds;
    }

    ds.y = v;

    ds.v =
        0.5 * (1.0 / y + 1.0 / (y - 1.0) + 1.0 / (y - t)) * v * v
        - (1.0 / t + 1.0 / (t - 1.0) + 1.0 / (y - t)) * v;

    return ds;
}

/* Runge-Kutta 8th Order Step */
State rk8_step(double t, const State& s, double h) {
    State k1 = painleveVI(t, s);

    State s2{ s.y + h * 0.5 * k1.y, s.v + h * 0.5 * k1.v };
    State k2 = painleveVI(t + 0.5 * h, s2);

    State s3{ s.y + h * 0.5 * k2.y, s.v + h * 0.5 * k2.v };
    State k3 = painleveVI(t + 0.5 * h, s3);

    State s4{ s.y + h * k3.y, s.v + h * k3.v };
    State k4 = painleveVI(t + h, s4);

    State out;
    out.y = s.y + h / 6.0 * (k1.y + 2*k2.y + 2*k3.y + k4.y);
    out.v = s.v + h / 6.0 * (k1.v + 2*k2.v + 2*k3.v + k4.v);

    return out;
}

/* Adaptive integration */
void integrate(double t0, double t1, State s0, double tol) {
    double t = t0;
    double h = 1e-3;
    State s = s0;

    cout << "\n t\t\t y(t)\t\t y'(t)\n";
    cout << "---------------------------------------------\n";

    while (t < t1) {
        if (t + h > t1) h = t1 - t;

        State full = rk8_step(t, s, h);
        State half1 = rk8_step(t, s, h / 2);
        State half2 = rk8_step(t + h / 2, half1, h / 2);

        double error =
            fabs(full.y - half2.y) +
            fabs(full.v - half2.v);

        if (error < tol || h <= MIN_STEP) {
            t += h;
            s = half2;

            cout << fixed << setprecision(8)
                 << t << "\t"
                 << s.y << "\t"
                 << s.v << "\n";

            double scale = SAFETY * pow(tol / (error + 1e-16), 0.125);
            h = min(MAX_STEP, h * scale);
        } else {
            h *= 0.5;
        }
    }
}

int main() {
    cout << "=============================================\n";
    cout << " Painlevé VI Transcendental Solver (CLI)\n";
    cout << " Adaptive Runge-Kutta Order 8\n";
    cout << "=============================================\n";

    char repeat;

    do {
        double t0, t1, y0, v0, tol;

        cout << "\nInitial t0        : "; cin >> t0;
        cout << "Final t1          : "; cin >> t1;
        cout << "Initial y(t0)     : "; cin >> y0;
        cout << "Initial y'(t0)    : "; cin >> v0;
        cout << "Error tolerance   : "; cin >> tol;

        State initial{ y0, v0 };

        integrate(t0, t1, initial, tol);

        cout << "\nRun another computation? (y/n): ";
        cin >> repeat;

    } while (repeat == 'y' || repeat == 'Y');

    cout << "\nProgram finished successfully.\n";
    return 0;
}
