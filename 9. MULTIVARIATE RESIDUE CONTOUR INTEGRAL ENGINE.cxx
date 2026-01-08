#include <iostream>
#include <complex>
#include <cmath>
#include <iomanip>

using namespace std;

static const double PI = acos(-1.0);

// ===============================
// Core Integrand
// ===============================
complex<double> integrand(double theta1, double theta2)
{
    complex<double> i(0.0, 1.0);

    complex<double> z1 = exp(i * theta1);
    complex<double> z2 = exp(i * theta2);

    complex<double> denominator =
        pow(z1, 3.0) + pow(z2, 3.0) + complex<double>(1.0, 0.0)
        - complex<double>(3.0, 0.0) * z1 * z2;

    // dz1 dz2 = (i z1)(i z2) dθ1 dθ2 = - z1 z2 dθ1 dθ2
    complex<double> jacobian = -z1 * z2;

    return jacobian / denominator;
}

// ===============================
// Numerical Integration on T²
// ===============================
complex<double> computeContourIntegral(int resolution)
{
    double dtheta = 2.0 * PI / resolution;
    complex<double> sum(0.0, 0.0);

    for (int i = 0; i < resolution; ++i)
    {
        for (int j = 0; j < resolution; ++j)
        {
            double t1 = i * dtheta;
            double t2 = j * dtheta;
            sum += integrand(t1, t2);
        }
    }

    return sum * dtheta * dtheta;
}

// ===============================
// CLI Interface
// ===============================
void runProgram()
{
    cout << "\n=============================================\n";
    cout << " MULTIVARIATE COMPLEX RESIDUE INTEGRATOR\n";
    cout << " Numerical Contour Integration on T²\n";
    cout << "=============================================\n";

    while (true)
    {
        int resolution;
        cout << "\nEnter angular resolution (e.g. 200, 400, 800): ";
        cin >> resolution;

        if (resolution <= 0)
        {
            cout << "Invalid resolution.\n";
            continue;
        }

        cout << "\nComputing integral...\n";

        complex<double> result = computeContourIntegral(resolution);

        cout << fixed << setprecision(10);
        cout << "\nApproximate Integral Value:\n";
        cout << "Real Part      : " << real(result) << "\n";
        cout << "Imaginary Part : " << imag(result) << "\n";
        cout << "Magnitude      : " << abs(result) << "\n";

        char choice;
        cout << "\nCompute another integral? (y/n): ";
        cin >> choice;

        if (choice != 'y' && choice != 'Y')
            break;
    }

    cout << "\nProgram terminated. Stay mathematical.\n";
}

// ===============================
// Main
// ===============================
int main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    runProgram();
    return 0;
}
