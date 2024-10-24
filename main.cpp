#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>
#include "Integration.h"
#define I_exact 1.2869386805747331527924009300176390931161637290256260886342156825 // взято из wolfram alpha

using namespace std;

double function_to_integrate(double x) {
    return exp(-x / 2) + x;
}

void buildGridsAndCalculate(const function<double(double)>& f, double a, double b) {
    vector<double> h_values = { 0.1, 0.05, 0.025 };
    NumericalIntegrator integrator(a, b, 0.1);

    // --- Трапеции ---
    vector<double> trapezoidal_results;
    cout << "\nTrapezoidal Results:\n" << scientific;

    for (double h : h_values) {
        integrator.set_h(h);
        double trapezoidal_result = integrator.trapezoidal(f);
        trapezoidal_results.push_back(trapezoidal_result);
        double error = trapezoidal_result - I_exact;
        cout << "h = " << h << " : " << trapezoidal_result << " Error: " << error << endl;
    }
    cout << "\nTrapezoidal Results Table:\n";
    cout << setw(15) << "h"
        << setw(15) << "k"
        << setw(15) << "Error"
        << setw(15) << "Error_norm"
        << setw(15) << "Runge"
        << setw(15) << "Richardson"
        << setw(15) << "Err_Richardson" << endl;
    for (int i = 1; i < trapezoidal_results.size(); i++) {
        double k = abs(log2(abs(1 + ((trapezoidal_results[i - 1] - trapezoidal_results[i]) / (I_exact - trapezoidal_results[i - 1])))));
        double err = I_exact - trapezoidal_results[i];
        double err_norm = (I_exact - trapezoidal_results[i]) / (I_exact - trapezoidal_results[i - 1]);
        double runge = (trapezoidal_results[i - 1] - trapezoidal_results[i]) / (pow(2, k) - 1);
        double richardson = trapezoidal_results[i - 1] + runge;
        double err_richardson = I_exact - richardson;

        cout << setw(15) << h_values[i]
            << setw(15) << k
            << setw(15) << err
            << setw(15) << err_norm
            << setw(15) << runge
            << setw(15) << richardson
            << setw(15) << err_richardson << endl;
    }

    // --- Гаусс-3 ---
    vector<double> gauss3_results;
    cout << "\nGauss-3 Results:\n" << scientific;
    for (double h : h_values) {
        integrator.set_h(h);
        double gauss3_result = integrator.gauss3(f);
        gauss3_results.push_back(gauss3_result);
        double error = gauss3_result - I_exact;
        cout << "h = " << h << " : " << gauss3_result << " Error: " << error << endl;
    }

    cout << "\nGauss-3 Results Table:\n";
    cout << setw(15) << "h"
        << setw(15) << "k"
        << setw(15) << "Error"
        << setw(15) << "Error_norm"
        << setw(15) << "Runge"
        << setw(15) << "Richardson"
        << setw(15) << "Err_Richardson" << endl;
    for (int i = 1; i < gauss3_results.size(); i++) {
        double k = abs(log2(abs(1 + (gauss3_results[i - 1] - gauss3_results[i]) / (I_exact - gauss3_results[i - 1]))));
        double err = I_exact - gauss3_results[i];
        double err_norm = (I_exact - gauss3_results[i]) / (I_exact - gauss3_results[i - 1]);
        double runge = (gauss3_results[i - 1] - gauss3_results[i]) / (pow(2, k) - 1);
        double richardson = gauss3_results[i - 1] + runge;
        double err_richardson = I_exact - richardson;

        cout << setw(15) << h_values[i]
            << setw(15) << k
            << setw(15) << err
            << setw(15) << err_norm
            << setw(15) << runge
            << setw(15) << richardson
            << setw(15) << err_richardson << endl;
    }
}

int main() {
    double a = 0.0, b = 1.0;
    buildGridsAndCalculate(function_to_integrate, a, b);
    return 0;
}
