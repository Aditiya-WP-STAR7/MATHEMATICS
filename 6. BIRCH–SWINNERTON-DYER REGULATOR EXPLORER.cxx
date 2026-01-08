#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>

using namespace std;

/*
 ============================================================
  Birch–Swinnerton-Dyer Regulator Explorer
  Elliptic Curve: y^2 = x^3 + 7823
  Author: Research-Grade Numerical Prototype
 ============================================================
*/

static const long double A = 0.0L;
static const long double B = 7823.0L;

/* ---------------- Elliptic Curve Point ---------------- */

struct ECPoint {
    long double x;
    long double y;
    bool infinity;

    ECPoint() : x(0), y(0), infinity(true) {}
    ECPoint(long double _x, long double _y) : x(_x), y(_y), infinity(false) {}
};

/* ---------------- Utility ---------------- */

long double naiveHeight(const ECPoint& P) {
    if (P.infinity) return 0.0L;
    return logl(max(1.0L, fabsl(P.x)));
}

/* ---------------- Elliptic Curve Arithmetic ---------------- */

ECPoint add(const ECPoint& P, const ECPoint& Q) {
    if (P.infinity) return Q;
    if (Q.infinity) return P;

    if (fabsl(P.x - Q.x) < 1e-18L && fabsl(P.y + Q.y) < 1e-18L)
        return ECPoint(); // Point at infinity

    long double lambda;

    if (fabsl(P.x - Q.x) < 1e-18L) {
        // Point doubling
        lambda = (3 * P.x * P.x + A) / (2 * P.y);
    } else {
        // Point addition
        lambda = (Q.y - P.y) / (Q.x - P.x);
    }

    long double xr = lambda * lambda - P.x - Q.x;
    long double yr = lambda * (P.x - xr) - P.y;

    return ECPoint(xr, yr);
}

ECPoint multiply(ECPoint P, int k) {
    ECPoint result;
    while (k > 0) {
        if (k & 1) result = add(result, P);
        P = add(P, P);
        k >>= 1;
    }
    return result;
}

/* ---------------- Canonical Height Approximation ---------------- */

long double canonicalHeight(ECPoint P, int iterations = 15) {
    long double height = 0.0L;
    ECPoint Q = P;

    for (int n = 1; n <= iterations; ++n) {
        Q = add(Q, Q); // 2^n P
        long double h = naiveHeight(Q);
        long double scaled = h / powl(4.0L, n);
        height = scaled;
    }

    return height;
}

/* ---------------- Main CLI Program ---------------- */

int main() {
    cout << fixed << setprecision(18);

    while (true) {
        cout << "\n==============================================\n";
        cout << " Birch–Swinnerton-Dyer Regulator Explorer\n";
        cout << " Elliptic Curve: y^2 = x^3 + 7823\n";
        cout << "==============================================\n";

        long double x, y;
        cout << "Enter x-coordinate of point P: ";
        cin >> x;

        cout << "Enter y-coordinate of point P: ";
        cin >> y;

        long double lhs = y * y;
        long double rhs = x * x * x + B;

        if (fabsl(lhs - rhs) > 1e-6L) {
            cout << "\n[ERROR] The point is NOT on the elliptic curve.\n";
            cout << "Please enter a valid point.\n";
            continue;
        }

        ECPoint P(x, y);

        cout << "\nComputing canonical height approximation...\n";

        long double h_hat = canonicalHeight(P);

        cout << "\n----------------------------------------------\n";
        cout << "Approximate Canonical Height (Néron–Tate):\n";
        cout << "ĥ(P) ≈ " << h_hat << "\n";
        cout << "----------------------------------------------\n";

        cout << "\nThis value contributes directly to the\n";
        cout << "Birch–Swinnerton-Dyer regulator R.\n";

        char choice;
        cout << "\nWould you like to compute another point? (y/n): ";
        cin >> choice;

        if (choice != 'y' && choice != 'Y') {
            cout << "\nExiting BSD Regulator Explorer.\n";
            cout << "Stay curious. Stay mathematical.\n";
            break;
        }
    }

    return 0;
}
