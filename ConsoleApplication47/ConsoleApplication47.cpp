#include <iostream>
#include <cmath>

using namespace std;

const int MAX_NODES = 100; // Максимальное количество узлов

struct CubicPolynomial {
    double a, b, c, d;
};

void solveTridiagonal(const double a[], const double b[], const double c[], const double d[], double x[], int n) {
    double ac[MAX_NODES], bc[MAX_NODES], dc[MAX_NODES];
    ac[0] = a[0];
    bc[0] = b[0];
    dc[0] = d[0];

    for (int i = 1; i < n; i++) {
        double m = ac[i - 1] / bc[i - 1];
        bc[i] = b[i] - m * c[i - 1];
        dc[i] = d[i] - m * dc[i - 1];
        ac[i] = a[i];
    }

    x[n - 1] = dc[n - 1] / bc[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = (dc[i] - c[i] * x[i + 1]) / bc[i];
    }
}

void computeSplineCoefficients(const double x[], const double y[], CubicPolynomial polynomials[], int n) {
    double h[MAX_NODES];
    for (int i = 0; i < n - 1; i++) {
        h[i] = x[i + 1] - x[i];
    }

    double a[MAX_NODES], b[MAX_NODES], c[MAX_NODES], d[MAX_NODES];
    for (int i = 0; i < n - 2; i++) {
        a[i] = h[i];
        b[i] = 2 * (h[i] + h[i + 1]);
        c[i] = h[i + 1];
        d[i] = 3 * ((y[i + 2] - y[i + 1]) / h[i + 1] - (y[i + 1] - y[i]) / h[i]);
    }

    double m[MAX_NODES];
    solveTridiagonal(a, b, c, d, m, n - 1);

    double secondDerivatives[MAX_NODES];
    secondDerivatives[0] = 0;
    for (int i = 0; i < n - 1; i++) {
        secondDerivatives[i + 1] = m[i];
    }
    secondDerivatives[n] = 0;

    for (int i = 0; i < n - 1; i++) {
        polynomials[i].a = (secondDerivatives[i + 1] - secondDerivatives[i]) / (6 * h[i]);
        polynomials[i].b = secondDerivatives[i] / 2;
        polynomials[i].c = (y[i + 1] - y[i]) / h[i] - h[i] * (2 * secondDerivatives[i] + secondDerivatives[i + 1]) / 6;
        polynomials[i].d = y[i];
    }
}

double interpolate(CubicPolynomial polynomials[], const double x[], int n, double t) {
    if (t < x[0] || t > x[n]) {
        return 0;
    }

    int i = 0;
    while (i < n && t > x[i + 1]) {
        ++i;
    }

    double dx = t - x[i];
    return polynomials[i].a * pow(dx, 3) + polynomials[i].b * pow(dx, 2) + polynomials[i].c * dx + polynomials[i].d;
}

int main() {
    setlocale(LC_ALL, "Ru");
    double x[MAX_NODES] = { 0, 1, 2, 3, 4 };
    double y[MAX_NODES] = { 1, 3, 2, 4, 0 };
    CubicPolynomial polynomials[MAX_NODES];

    computeSplineCoefficients(x, y, polynomials, 5);

    for (double t = 0.0; t <= 4.0; t += 0.5) {
        double interpolatedValue = interpolate(polynomials, x, 5, t);
        cout << "t = " << t << ", интерполированное значение = " << interpolatedValue << endl;
    }

    return 0;
}