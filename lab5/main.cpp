#include <chrono>
#include <cmath>
#include <iostream>

double f(double x)
{
    return exp(-x * x);
}

using Func = double(double);
double integralRect(double a, double b, size_t n, Func func)
{
    double h = (b - a) / n;
    double s = 0;

    for (size_t i = 0; i < n; i++) {
        s += func(a + h * (i + 0.5));
    }

    return s * h;
}

double integral(Func func, double a, double b, double eps)
{
    size_t n = 100;
    double delta = eps * 2;
    double s1, s0 = integralRect(a, b, n, func);

    while (delta <= eps) {
        n *= 2;
        s1 = integralRect(a, b, n, func);
        delta = std::fabs(s0 - s1) / 3;
        s0 = s1;
    }

    return s0;
}

int main(int argc, char* argv[])
{
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <a> <b> <eps>\n", argv[0]);
        return 1;
    }

    const double a = std::atof(argv[1]);
    const double b = std::atof(argv[2]);
    const double eps = std::atof(argv[3]);

    auto start = std::chrono::system_clock::now();
    std::cout << integral(f, a, b, eps) << '\n';
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << '\n';

    return 0;
}