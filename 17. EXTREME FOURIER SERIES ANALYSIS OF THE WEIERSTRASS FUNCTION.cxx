#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>

using namespace std;

using ld = long double;
using cd = complex<ld>;

const ld PI = acosl(-1.0L);

/* ===========================
   FFT IMPLEMENTATION
   =========================== */
void fft(vector<cd>& a, bool invert) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j |= bit;
        if (i < j) swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        ld ang = 2 * PI / len * (invert ? -1 : 1);
        cd wlen(cosl(ang), sinl(ang));
        for (int i = 0; i < n; i += len) {
            cd w(1);
            for (int j = 0; j < len / 2; j++) {
                cd u = a[i + j];
                cd v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        for (cd& x : a) x /= n;
    }
}

/* ===========================
   WEIERSTRASS FUNCTION
   =========================== */
ld weierstrass(ld x, ld a, int b, int terms) {
    ld sum = 0.0L;
    for (int n = 0; n < terms; n++) {
        sum += powl(a, n) * cosl(powl(b, n) * PI * x);
    }
    return sum;
}

/* ===========================
   FOURIER COEFFICIENT VIA
   OSCILLATORY INTEGRAL
   =========================== */
ld fourierCoefficient(int k, ld a, int b, int terms, int samples) {
    ld integral = 0.0L;
    ld dx = 2.0L / samples;

    for (int i = 0; i < samples; i++) {
        ld x = -1.0L + i * dx;
        ld fx = weierstrass(x, a, b, terms);
        integral += fx * cosl(PI * k * x);
    }

    return integral * dx;
}

/* ===========================
   MAIN PROGRAM LOOP
   =========================== */
int main() {
    cout << fixed << setprecision(15);

    while (true) {
        ld a;
        int b, terms, samples, maxK;

        cout << "\n=== Extreme Fourier Analysis of Weierstrass Function ===\n";
        cout << "Enter parameter a (0 < a < 1): ";
        cin >> a;
        cout << "Enter integer b (>1): ";
        cin >> b;
        cout << "Number of Weierstrass terms: ";
        cin >> terms;
        cout << "Numerical integration samples: ";
        cin >> samples;
        cout << "Max Fourier mode k: ";
        cin >> maxK;

        cout << "\nComputing Fourier coefficients...\n";

        vector<ld> coeff(maxK + 1);
        for (int k = 0; k <= maxK; k++) {
            coeff[k] = fourierCoefficient(k, a, b, terms, samples);
            cout << "a_" << k << " = " << coeff[k] << "\n";
        }

        cout << "\nComputation complete.\n";
        cout << "Run another computation? (y/n): ";
        char choice;
        cin >> choice;
        if (choice != 'y' && choice != 'Y') break;
    }

    cout << "\nProgram terminated. Keep exploring the impossible.\n";
    return 0;
}
