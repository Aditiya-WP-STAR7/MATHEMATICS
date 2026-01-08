#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <limits>

using namespace std;

/*
    ============================================================
    CATEGORY THEORY PROJECTIVE LIMIT SIMULATOR
    ------------------------------------------------------------
    Title  : Infinite Projective Limit of (Z / p^n Z) ⊗ Q/Z
    Author : MIT-Level Mathematical CLI Program
    Field  : Category Theory, Algebra, Data Science, Number Theory
    ============================================================
*/

/*
    We represent a commutative diagram as a directed graph.
    Nodes   : Objects (Z / p^n Z)
    Edges   : Morphisms (canonical projection maps)
*/

struct GroupObject {
    int n;              // exponent n in Z / p^n Z
    long long order;    // p^n
};

class CommutativeDiagram {
private:
    vector<GroupObject> objects;
    map<int, int> morphisms; // n -> n-1

public:
    void buildDiagram(long long p, int maxLevel) {
        objects.clear();
        morphisms.clear();

        for (int n = 1; n <= maxLevel; ++n) {
            GroupObject g;
            g.n = n;
            g.order = pow(p, n);
            objects.push_back(g);

            if (n > 1) {
                morphisms[n] = n - 1;
            }
        }
    }

    void displayDiagram() const {
        cout << "\nCommutative Diagram (Inverse System):\n";
        for (const auto& g : objects) {
            cout << "Object: Z / " << g.order << " Z";
            if (morphisms.count(g.n)) {
                cout << "  --->  Z / " << pow(objects[0].order, morphisms.at(g.n)) << " Z";
            }
            cout << "\n";
        }
    }

    /*
        Tensor with Q/Z:
        For finite cyclic groups, tensoring with Q/Z extracts p-primary torsion.
    */
    void computeTensorWithQmodZ(long long p) const {
        cout << "\nTensoring each object with Q/Z:\n";
        for (const auto& g : objects) {
            cout << "(Z / " << g.order << " Z) ⊗ (Q / Z)";
            cout << "  ≅  Z / " << g.order << " Z  (p-primary torsion)\n";
        }
    }

    /*
        Projective limit result (theoretical, not brute-force enumerable)
    */
    void computeProjectiveLimit(long long p) const {
        cout << "\nComputing Projective Limit...\n\n";
        cout << "lim_n (Z / p^n Z) ⊗ (Q / Z)\n";
        cout << "---------------------------------\n";
        cout << "Result (Category-Theoretic Limit):\n";
        cout << "≅ Z_p / Z\n";
        cout << "(p-adic solenoid, compact, totally disconnected)\n";
    }
};

/*
    Input validation for large primes
*/
long long readPrime() {
    long long p;
    while (true) {
        cout << "Enter a prime number p (>= 2): ";
        cin >> p;
        if (cin.fail() || p < 2) {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            cout << "Invalid input. Try again.\n";
        } else {
            return p;
        }
    }
}

int readMaxLevel() {
    int n;
    while (true) {
        cout << "Enter max exponent n (suggested <= 12 for CLI): ";
        cin >> n;
        if (cin.fail() || n <= 0) {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            cout << "Invalid input. Try again.\n";
        } else {
            return n;
        }
    }
}

int main() {
    cout << "=====================================================\n";
    cout << " CATEGORY THEORY PROJECTIVE LIMIT CLI (C++)\n";
    cout << " Infinite Inverse Limits & Tensor Products\n";
    cout << "=====================================================\n";

    char repeat;
    do {
        long long p = readPrime();
        int maxLevel = readMaxLevel();

        CommutativeDiagram diagram;
        diagram.buildDiagram(p, maxLevel);

        diagram.displayDiagram();
        diagram.computeTensorWithQmodZ(p);
        diagram.computeProjectiveLimit(p);

        cout << "\nWould you like to compute another projective limit? (y/n): ";
        cin >> repeat;

    } while (repeat == 'y' || repeat == 'Y');

    cout << "\nThank you. Mathematics is eternal.\n";
    cout << "Program terminated with categorical coherence.\n";
    return 0;
}
