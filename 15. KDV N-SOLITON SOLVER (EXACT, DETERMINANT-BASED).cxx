#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>

using namespace std;

/* ===============================
   HIGH PRECISION TYPE
   =============================== */
using Real = long double;

/* ===============================
   DENSE MATRIX CLASS
   =============================== */
class Matrix {
public:
    int n;
    vector<vector<Real>> a;

    Matrix(int n_) : n(n_), a(n_, vector<Real>(n_, 0.0L)) {}

    static Matrix identity(int n) {
        Matrix I(n);
        for (int i = 0; i < n; ++i) I.a[i][i] = 1.0L;
        return I;
    }

    Real determinant() const {
        Matrix m = *this;
        Real det = 1.0L;

        for (int i = 0; i < n; ++i) {
            int pivot = i;
            for (int j = i + 1; j < n; ++j)
                if (fabsl(m.a[j][i]) > fabsl(m.a[pivot][i]))
                    pivot = j;

            if (fabsl(m.a[pivot][i]) < 1e-18L)
                return 0.0L;

            if (pivot != i) {
                swap(m.a[pivot], m.a[i]);
                det = -det;
            }

            det *= m.a[i][i];
            Real inv = 1.0L / m.a[i][i];

            for (int j = i + 1; j < n; ++j) {
                Real factor = m.a[j][i] * inv;
                for (int k = i; k < n; ++k)
                    m.a[j][k] -= factor * m.a[i][k];
            }
        }
        return det;
    }
};

/* ===============================
   BUILD A MATRIX
   =============================== */
Matrix buildA(const vector<Real>& p,
              const vector<Real>& c,
              Real x, Real t) {

    int N = p.size();
    Matrix A(N);

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            Real exponent =
                (p[i] + p[j]) * x -
                (p[i]*p[i]*p[i] + p[j]*p[j]*p[j]) * t;

            A.a[i][j] =
                (c[i] * c[j]) / (p[i] + p[j]) *
                expl(exponent);
        }
    return A;
}

/* ===============================
   LOG DET(I + A)
   =============================== */
Real logDetIA(const vector<Real>& p,
              const vector<Real>& c,
              Real x, Real t) {

    Matrix A = buildA(p, c, x, t);
    Matrix I = Matrix::identity(p.size());

    for (int i = 0; i < A.n; ++i)
        for (int j = 0; j < A.n; ++j)
            A.a[i][j] += I.a[i][j];

    Real det = A.determinant();
    if (det <= 0)
        throw runtime_error("Non-positive determinant encountered.");

    return logl(det);
}

/* ===============================
   SECOND DERIVATIVE (5-POINT)
   =============================== */
Real secondDerivative(const vector<Real>& p,
                      const vector<Real>& c,
                      Real x, Real t) {

    const Real h = 1e-4L;

    Real f1 = logDetIA(p, c, x - 2*h, t);
    Real f2 = logDetIA(p, c, x - h, t);
    Real f3 = logDetIA(p, c, x, t);
    Real f4 = logDetIA(p, c, x + h, t);
    Real f5 = logDetIA(p, c, x + 2*h, t);

    return (-f5 + 16*f4 - 30*f3 + 16*f2 - f1) / (12*h*h);
}

/* ===============================
   MAIN SOLITON SOLVER
   =============================== */
Real kdvSolution(const vector<Real>& p,
                 const vector<Real>& c,
                 Real x, Real t) {
    return -2.0L * secondDerivative(p, c, x, t);
}

/* ===============================
   CLI APPLICATION
   =============================== */
int main() {
    cout << fixed << setprecision(12);
    cout << "\n=== KdV N-Soliton Exact Solver ===\n";
    cout << "Integrable Systems | Determinant Method\n";

    while (true) {
        int N;
        cout << "\nEnter number of solitons N: ";
        cin >> N;

        vector<Real> p(N), c(N);
        for (int i = 0; i < N; ++i) {
            cout << "p[" << i << "], c[" << i << "]: ";
            cin >> p[i] >> c[i];
        }

        Real x, t;
        cout << "Enter x and t: ";
        cin >> x >> t;

        try {
            Real u = kdvSolution(p, c, x, t);
            cout << "\nExact KdV Solution u(x,t) = " << u << "\n";
        } catch (const exception& e) {
            cout << "Error: " << e.what() << "\n";
        }

        char again;
        cout << "\nCompute another case? (y/n): ";
        cin >> again;
        if (again != 'y' && again != 'Y')
            break;
    }

    cout << "\nProgram terminated professionally.\n";
    return 0;
}
