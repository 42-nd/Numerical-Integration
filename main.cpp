#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>
#include "Integration.h"
#define I_exact 1.2869387

using namespace std;
double function_to_integrate(double x) {
    return exp(-x / 2) + x;
}

void buildGridsAndCalculate(const function<double(double)>& f, double a, double b) {
    vector<double> h_values = {0.1, 0.05, 0.025}; 
    NumericalIntegrator integrator(a, b, 0.1);

    vector<double> trapezoidal_results;
    cout << "\nTrapezoidal Results:\n" << scientific;
    for (double h : h_values) {
        integrator.set_h(h);
        double trapezoidal_result = integrator.trapezoidal(f);
        trapezoidal_results.push_back(trapezoidal_result);
        double error = trapezoidal_result - I_exact;
        cout << "h = " << h << " : " << trapezoidal_result << " Error: "  << error << endl;
    }
    cout << "\nTrapezoidal Results:\n";
    cout << setw(8) << "h"
        << setw(12) << "k"
        << setw(15) << "Error"
        << setw(15) << "Error_norm"
        << setw(15) << "Runge"
        << setw(15) << "Richardson"
        << setw(15) << "Err_Richardson" << endl;
    for (int i = 1; i < trapezoidal_results.size(); i++) {
        double k = abs(log2(1 + (trapezoidal_results[i - 1] - trapezoidal_results[i]) / (I_exact - trapezoidal_results[i - 1])));
        double err = I_exact - trapezoidal_results[i];
        double err_norm = (I_exact - trapezoidal_results[i]) / (I_exact - trapezoidal_results[i - 1]);
        double runge = (trapezoidal_results[i - 1] - trapezoidal_results[i]) / (pow(2, k) - 1);
        double richardson = trapezoidal_results[i - 1] + runge;
        double err_richardson = I_exact - richardson;

        cout << setw(8) << h_values[i]
            << setw(12) << k
            << setw(15) << err
            << setw(15) << err_norm
            << setw(15) << runge
            << setw(15) << richardson
            << setw(15) << err_richardson << endl;


    }

    double gauss3_result = integrator.gauss3(f);
    double gauss3_error = gauss3_result - I_exact;

    cout << "\nGauss-3:\n";
    cout << "Result: " << gauss3_result << " Error: " << scientific << gauss3_error << endl;
}

int main() {
    double a = 0.0, b = 1.0;
    buildGridsAndCalculate(function_to_integrate, a, b);
    return 0;
}