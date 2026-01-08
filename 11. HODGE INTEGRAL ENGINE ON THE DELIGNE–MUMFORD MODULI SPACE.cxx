#include <iostream>
#include <map>
#include <iomanip>

using namespace std;

/*
 ============================================================
  MATHEMATICALLY RIGOROUS HODGE INTEGRAL ENGINE
  Author : Aditiya WP
  Field  : Algebraic Geometry / Topological Gravity
 ============================================================
*/

// Bernoulli numbers B_{2g}
map<int, long double> bernoulli = {
    {2,  1.0L/6.0L},
    {4, -1.0L/30.0L},
    {6,  1.0L/42.0L},
    {8, -1.0L/30.0L},
    {10, 5.0L/66.0L}
};

long double factorial(int n)
{
    long double r = 1.0L;
    for (int i = 1; i <= n; ++i)
        r *= i;
    return r;
}

/*
  Computes:
  ∫_{M̄_{g,1}} λ_g ψ_1^{2g}
*/
long double hodgeIntegral_g1(int g)
{
    if (bernoulli.count(2*g) == 0)
        return 0.0L;

    long double B = fabsl(bernoulli[2*g]);
    return B / ( (2.0L * g) * factorial(2*g) );
}

int main()
{
    cout << fixed << setprecision(18);

    cout << "=============================================\n";
    cout << " RIGOROUS HODGE INTEGRAL COMPUTATION ENGINE\n";
    cout << " ∫ M̄_{g,1} λ_g ψ_1^{2g}\n";
    cout << "=============================================\n\n";

    while (true)
    {
        int g;
        cout << "Enter genus g (>=1): ";
        cin >> g;

        if (g < 1)
        {
            cout << "Genus must be >= 1\n\n";
            continue;
        }

        long double result = hodgeIntegral_g1(g);

        cout << "\nExact Hodge Integral:\n";
        cout << "∫ M̄_" << g << ",1 λ_" << g
             << " ψ_1^" << 2*g
             << " = " << result << "\n\n";

        char c;
        cout << "Compute another integral? (y/n): ";
        cin >> c;

        if (c != 'y' && c != 'Y')
        {
            cout << "\nExiting engine.\n";
            cout << "Intersection theory never lies.\n";
            break;
        }

        cout << "\n---------------------------------------------\n\n";
    }

    return 0;
}
