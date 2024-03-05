#include <cmath>
#include <iostream>
#include <tuple>

double f(double x)
{
    return 3 - 1.0 / 2.0 * sqrt(x) - exp(-1.0 / 2.0 * x * x);
}

auto solution(double a, double b, double eps)
{
    double x = (a + b) / 2;
    double x_new = x - (std::pow(f(x), 2)) / (f(x + f(x)) - f(x));
    size_t i = 1;

    while (fabs(x_new - x) > eps) {
        x = x_new;
        x_new = x - (std::pow(f(x), 2)) / (f(x + f(x)) - f(x));
        i++;
    }

    return std::make_tuple(x_new, i);
}

int main(int argc, char* argv[])
{
    if (argc < 4) {
        return 1;
    }
    double a = std::atof(argv[1]);
    double b = std::atof(argv[2]);
    double eps = std::atof(argv[3]);

    auto [x, i] = solution(a, b, eps);
    std::cout << "x = " << x << " iter = " << i << '\n';

    return 0;
}