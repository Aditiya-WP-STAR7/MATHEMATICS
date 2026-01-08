#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

constexpr int DIM = 6;

// ================= VECTOR & MATRIX UTILITIES =================

typedef vector<double> Vec;
typedef vector<Vec> Mat;

Vec zeros() {
    return Vec(DIM, 0.0);
}

Mat identity() {
    Mat I(DIM, Vec(DIM, 0.0));
    for (int i = 0; i < DIM; ++i) I[i][i] = 1.0;
    return I;
}

// ================= GRAM–SCHMIDT QR DECOMPOSITION =================

void qrDecomposition(const Mat& A, Mat& Q, Mat& R) {
    Q = identity();
    R = Mat(DIM, Vec(DIM, 0.0));

    vector<Vec> v = A;

    for (int i = 0; i < DIM; ++i) {
        for (int j = 0; j < i; ++j) {
            double dot = 0.0;
            for (int k = 0; k < DIM; ++k)
                dot += v[i][k] * Q[j][k];

            R[j][i] = dot;
            for (int k = 0; k < DIM; ++k)
                v[i][k] -= dot * Q[j][k];
        }

        double norm = 0.0;
        for (double x : v[i]) norm += x * x;
        norm = sqrt(norm);

        R[i][i] = norm;

        for (int k = 0; k < DIM; ++k)
            Q[i][k] = v[i][k] / norm;
    }
}

// ================= HÉNON 6D MAP =================

Vec henon6D(const Vec& X, double a, double b, double c,
            double d, double e, double f, double eps) {
    Vec Y(DIM);

    Y[0] = 1.0 - a * X[0] * X[0] + X[1] + eps * X[2];
    Y[1] = b * X[0];
    Y[2] = c * X[1];
    Y[3] = d * X[2];
    Y[4] = e * X[3];
    Y[5] = f * X[4];

    return Y;
}

// ================= JACOBIAN MATRIX =================

Mat jacobian(const Vec& X, double a, double b, double c,
             double d, double e, double f, double eps) {
    Mat J(DIM, Vec(DIM, 0.0));

    J[0][0] = -2.0 * a * X[0];
    J[0][1] = 1.0;
    J[0][2] = eps;

    J[1][0] = b;
    J[2][1] = c;
    J[3][2] = d;
    J[4][3] = e;
    J[5][4] = f;

    return J;
}

// ================= LYAPUNOV EXPONENTS =================

Vec lyapunovSpectrum(Vec X, int iterations,
                     double a, double b, double c,
                     double d, double e, double f, double eps) {

    Mat Q = identity();
    Vec lambda(DIM, 0.0);

    for (int i = 0; i < iterations; ++i) {
        Mat J = jacobian(X, a, b, c, d, e, f, eps);

        Mat A(DIM, Vec(DIM, 0.0));
        for (int r = 0; r < DIM; ++r)
            for (int k = 0; k < DIM; ++k)
                for (int c2 = 0; c2 < DIM; ++c2)
                    A[r][k] += J[r][c2] * Q[c2][k];

        Mat R;
        qrDecomposition(A, Q, R);

        for (int j = 0; j < DIM; ++j)
            lambda[j] += log(abs(R[j][j]));

        X = henon6D(X, a, b, c, d, e, f, eps);
    }

    for (double& l : lambda)
        l /= iterations;

    return lambda;
}

// ================= MAIN CLI =================

int main() {
    cout << fixed << setprecision(10);

    while (true) {
        Vec X(DIM);
        double a, b, c, d, e, f, eps;
        int iterations;

        cout << "\n=== 6D Generalized Hénon Chaos Analyzer ===\n";
        cout << "Enter initial conditions (x y z u v w): ";
        for (double& x : X) cin >> x;

        cout << "Enter parameters a b c d e f epsilon: ";
        cin >> a >> b >> c >> d >> e >> f >> eps;

        cout << "Number of iterations: ";
        cin >> iterations;

        Vec lyap = lyapunovSpectrum(X, iterations, a, b, c, d, e, f, eps);

        cout << "\nLyapunov Exponents:\n";
        for (int i = 0; i < DIM; ++i)
            cout << "λ" << i + 1 << " = " << lyap[i] << endl;

        cout << "\nRun another computation? (y/n): ";
        char choice;
        cin >> choice;
        if (choice != 'y' && choice != 'Y') break;
    }

    cout << "\nProgram terminated. Stay chaotic.\n";
    return 0;
}
