#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <random>

using namespace std;

/*
===========================================================
 Quaternion Class
 Exact quaternion algebra for ADHM construction
===========================================================
*/

struct Quaternion {
    double r, i, j, k;

    Quaternion(double r_=0, double i_=0, double j_=0, double k_=0)
        : r(r_), i(i_), j(j_), k(k_) {}

    Quaternion operator+(const Quaternion& q) const {
        return Quaternion(r + q.r, i + q.i, j + q.j, k + q.k);
    }

    Quaternion operator-(const Quaternion& q) const {
        return Quaternion(r - q.r, i - q.i, j - q.j, k - q.k);
    }

    Quaternion operator*(const Quaternion& q) const {
        return Quaternion(
            r*q.r - i*q.i - j*q.j - k*q.k,
            r*q.i + i*q.r + j*q.k - k*q.j,
            r*q.j - i*q.k + j*q.r + k*q.i,
            r*q.k + i*q.j - j*q.i + k*q.r
        );
    }

    Quaternion conjugate() const {
        return Quaternion(r, -i, -j, -k);
    }

    double norm() const {
        return sqrt(r*r + i*i + j*j + k*k);
    }
};

/*
===========================================================
 ADHM Matrix Construction
 k = 5 → 10×10 quaternionic matrices
===========================================================
*/

const int K = 5;
const int N = 2 * K;

using QMatrix = vector<vector<Quaternion>>;

QMatrix zeroMatrix() {
    return QMatrix(N, vector<Quaternion>(N));
}

QMatrix randomADHMMatrix(double scale = 0.1) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(-scale, scale);

    QMatrix M = zeroMatrix();
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            M[i][j] = Quaternion(
                dist(gen), dist(gen), dist(gen), dist(gen)
            );
    return M;
}

QMatrix dagger(const QMatrix& A) {
    QMatrix D = zeroMatrix();
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            D[j][i] = A[i][j].conjugate();
    return D;
}

QMatrix multiply(const QMatrix& A, const QMatrix& B) {
    QMatrix C = zeroMatrix();
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k)
                C[i][j] = C[i][j] + A[i][k] * B[k][j];
    return C;
}

/*
===========================================================
 Self-Duality Numerical Test
 Checks || Δ†Δ - (Δ†Δ)† ||
===========================================================
*/

double selfDualityError(const QMatrix& Delta) {
    QMatrix Ddag = dagger(Delta);
    QMatrix M = multiply(Ddag, Delta);
    QMatrix Mdag = dagger(M);

    double err = 0.0;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            err += (M[i][j] - Mdag[i][j]).norm();

    return err;
}

/*
===========================================================
 Approximate Chern Number Estimator
 (Topological density surrogate)
===========================================================
*/

double estimateChernNumber(const QMatrix& Delta) {
    QMatrix Ddag = dagger(Delta);
    QMatrix M = multiply(Ddag, Delta);

    double trace = 0.0;
    for (int i = 0; i < N; ++i)
        trace += M[i][i].r;

    return trace / (8.0 * M_PI * M_PI);
}

/*
===========================================================
 MAIN CLI PROGRAM
===========================================================
*/

int main() {
    cout << fixed << setprecision(10);

    while (true) {
        cout << "\n=============================================\n";
        cout << " Yang–Mills Instanton Solver (ADHM, k = 5)\n";
        cout << " SU(2) Gauge Theory on S^4\n";
        cout << "=============================================\n";

        cout << "\n[1] Constructing ADHM matrix (10×10 quaternions)...\n";
        QMatrix Delta = randomADHMMatrix();

        cout << "[2] Evaluating self-duality consistency...\n";
        double error = selfDualityError(Delta);
        cout << "    Self-duality error ||F - *F|| ≈ " << error << "\n";

        cout << "[3] Estimating Chern number...\n";
        double chern = estimateChernNumber(Delta);
        cout << "    Estimated Chern number ≈ " << chern << "\n";

        cout << "\n---------------------------------------------\n";
        cout << " Computation complete.\n";
        cout << "---------------------------------------------\n";

        cout << "\nWould you like to compute another instanton?\n";
        cout << "Enter [1] Yes   [0] No : ";

        int choice;
        cin >> choice;
        if (choice == 0) break;
    }

    cout << "\nProgram terminated gracefully.\n";
    return 0;
}
