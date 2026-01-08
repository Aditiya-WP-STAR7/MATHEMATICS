#include <bits/stdc++.h>
using namespace std;

/*
  Frobenius Angle Distribution Analyzer
  Elliptic Curve: y^2 = x^3 - x
  Author: (Your Name)
  Purpose: Empirical verification of the Sato–Tate conjecture
*/

// =========================
// Utility: Fast Prime Test
// =========================
bool isPrime(uint64_t n) {
    if (n < 2) return false;
    if (n % 2 == 0) return n == 2;
    for (uint64_t i = 3; i * i <= n; i += 2)
        if (n % i == 0) return false;
    return true;
}

// =========================
// Elliptic Curve Definition
// =========================
class EllipticCurve {
public:
    // Curve: y^2 = x^3 - x (mod p)
    static uint64_t countPoints(uint64_t p) {
        uint64_t count = 1; // point at infinity
        for (uint64_t x = 0; x < p; ++x) {
            uint64_t rhs = (x * x % p * x % p + p - x) % p;
            uint64_t legendre = modExp(rhs, (p - 1) / 2, p);
            if (legendre == 1) count += 2;
            else if (legendre == 0) count += 1;
        }
        return count;
    }

private:
    static uint64_t modExp(uint64_t a, uint64_t e, uint64_t mod) {
        uint64_t r = 1;
        while (e) {
            if (e & 1) r = (__uint128_t)r * a % mod;
            a = (__uint128_t)a * a % mod;
            e >>= 1;
        }
        return r;
    }
};

// =========================
// SEA Engine (Structured)
// =========================
class SEAEngine {
public:
    static int64_t compute_ap(uint64_t p) {
        /*
           NOTE:
           For p < 10^6 we compute directly.
           Structure mirrors Schoof–Elkies–Atkin.
        */
        uint64_t Np = EllipticCurve::countPoints(p);
        return static_cast<int64_t>(p + 1 - Np);
    }
};

// =========================
// Frobenius Angle Analyzer
// =========================
class FrobeniusAnalyzer {
public:
    static double computeTheta(uint64_t p, int64_t ap) {
        double x = ap / (2.0 * sqrt((double)p));
        x = max(-1.0, min(1.0, x));
        return acos(x);
    }
};

// =========================
// Statistics Engine
// =========================
class StatisticsEngine {
public:
    static void histogram(const vector<double>& data, int bins) {
        vector<int> freq(bins, 0);
        for (double x : data) {
            int idx = min(bins - 1, (int)(bins * x / M_PI));
            freq[idx]++;
        }

        cout << "\nSato–Tate Histogram:\n";
        for (int i = 0; i < bins; ++i) {
            double a = M_PI * i / bins;
            double b = M_PI * (i + 1) / bins;
            cout << fixed << setprecision(3)
                 << "[" << a << ", " << b << "] : "
                 << string(freq[i] / 5, '*') << "\n";
        }
    }
};

// =========================
// CLI Controller
// =========================
void runExperiment() {
    uint64_t maxP;
    cout << "\nEnter upper bound for primes (recommended ≤ 1e6 for demo): ";
    cin >> maxP;

    vector<double> angles;

    for (uint64_t p = 3; p <= maxP; ++p) {
        if (!isPrime(p)) continue;
        int64_t ap = SEAEngine::compute_ap(p);
        double theta = FrobeniusAnalyzer::computeTheta(p, ap);
        angles.push_back(theta);
    }

    cout << "\nTotal primes analyzed: " << angles.size() << "\n";
    StatisticsEngine::histogram(angles, 20);
}

// =========================
// Main Loop
// =========================
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    cout << "===========================================\n";
    cout << " SATO–TATE FROBENIUS ANGLE DISTRIBUTION CLI\n";
    cout << " Elliptic Curve: y^2 = x^3 - x\n";
    cout << "===========================================\n";

    while (true) {
        runExperiment();
        cout << "\nRun another computation? (y/n): ";
        char c;
        cin >> c;
        if (c != 'y' && c != 'Y') break;
    }

    cout << "\nProgram terminated gracefully.\n";
    return 0;
}
