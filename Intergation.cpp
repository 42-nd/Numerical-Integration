#include "Integration.h"
using namespace std;

double NumericalIntegrator::trapezoidal(const function<double(double)>& f) {
    double result = h * f(a) + h * f(b);

    for (int i = 1; i < grid.size() - 1; ++i) {
        double x = a + i * h;
        result += (h * 2) * f(x);
    }
    return result * 0.5;
}

double NumericalIntegrator::gauss3(const function<double(double)>& f) {
    const double gauss_points[] = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
    const double gauss_weights[] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

    double midpoint = (a + b) / 2;
    double half_length = (b - a) / 2;
    double result = 0.0;

    for (int i = 0; i < 3; ++i) {
        result += gauss_weights[i] * f(midpoint + half_length * gauss_points[i]);
    }
    return result * half_length;
}

void NumericalIntegrator::buildGrid() {
    double step = a + h;
    grid.clear();
    grid.push_back(a);
    while (step < b) {
        grid.push_back(step);
        step += h;
    }
    grid.push_back(b);
}

double NumericalIntegrator::get_h() {
    return h;
}

void NumericalIntegrator::set_h(double new_h) {
    h = new_h;
    buildGrid();
}