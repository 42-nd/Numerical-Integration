#include "Integration.h"
using namespace std;

double NumericalIntegrator::trapezoidal(const function<double(double)>& f) {
    double result = h * f(a) + h * f(b);
    for (int i = 1; i < grid.size()-1; i++) {
        double x = grid[i];
        result += (h + h) * f(x);
    }
    return result * 0.5;
}

double NumericalIntegrator::gauss3(const function<double(double)>& f) {
    const double gauss_points[] = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
    const double gauss_weights[] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };


    double result = 0.0;
    for (int k = 1; k < grid.size(); k++) {
        double x_left = grid[k - 1];
        double x_right = grid[k];
        double half_h = (x_right - x_left) / 2.0;

        for (int i = 0; i < 3; i++) {
            double x = (x_left + x_right) / 2.0 + half_h * gauss_points[i];
            result += gauss_weights[i] * f(x);
        }
    }

    return result * (h / 2.0);
}

void NumericalIntegrator::buildGrid() {
    double step = a + h;
    grid.clear();
    grid.push_back(a);
    while (abs(step - b) >= 1e-6) {
        grid.push_back(step);
        step += h;
    }
    if (grid.back() != b) {
        grid.push_back(b);
    }

}

double NumericalIntegrator::get_h() {
    return h;
}

void NumericalIntegrator::set_h(double new_h) {
    h = new_h;
    buildGrid();
}