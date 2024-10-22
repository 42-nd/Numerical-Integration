#pragma once
#ifndef integration_h
#define integration_h

#include <vector>
#include <cmath>
#include <functional>
#include <iostream>
using namespace std;
class NumericalIntegrator {
public:
    NumericalIntegrator(double a, double b, double h) : a(a), b(b), h(h) {
        buildGrid();
    }

    // ����� ��� �������������� ������� ��������
    double trapezoidal(const function<double(double)>& f);
    // ����� ��� �������������� � ������� ���������� ������-3
    double gauss3(const function<double(double)>& f);
    double get_h();
    void set_h(double new_h);
private:
    double a, b, h;
    vector<double> grid;

    void buildGrid();
};
#endif