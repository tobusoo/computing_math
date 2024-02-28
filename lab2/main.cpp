#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <tuple>

void print_maxtrix(double** a, double* x, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n - 1; j++) {
            std::cout << a[i][j] << "x" << j + 1 << " + ";
            // std::cout << a[i][j] << ' ';
        }
        // std::cout << a[i][n - 1] << ' ';
        std::cout << a[i][n - 1] << "x" << n << " = ";
        if (x == NULL)
            std::cout << '\n';
        else
            std::cout << x[i] << '\n';
    }
    std::cout << '\n';
}

bool is_diagonal_predominance_h(double** a, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        double sum = 0;
        for (size_t j = 0; j < n; j++)
            sum += std::fabs(a[i][j]);
        sum -= std::fabs(a[i][i]);
        if (sum > std::fabs(a[i][i]))
            return false;
    }

    return true;
}

bool is_diagonal_predominance(double** a, double* x, size_t n)
{
    if (is_diagonal_predominance_h(a, n)) {
        return true;
    }

    for (size_t i = 0; i < n; i++) {
        double sum = 0;
        for (size_t j = 0; j < n; j++)
            sum += std::fabs(a[i][j]);
        sum -= std::fabs(a[i][i]);
        if (sum > std::fabs(a[i][i])) {
            bool swap = false;
            // std::cout << i << '\n';
            size_t l = i;
            for (size_t k = i + 1; k < n; k++) {
                sum = 0;
                for (size_t j = 0; j < n; j++)
                    sum += std::fabs(a[k][j]);
                sum -= std::fabs(a[k][l]);
                // std::cout << sum << ' ' << std::fabs(a[k][l]) << '\n';
                if (sum <= std::fabs(a[k][l])) {
                    std::swap(a[i], a[k]);
                    std::swap(x[i], x[k]);
                    swap = true;
                    break;
                }
            }
            if (!swap) {
                // print_maxtrix(a, x, n);
                return 0;
            }
            // print_maxtrix(a, x, n);
        }
    }

    // print_maxtrix(a, x, n);

    return is_diagonal_predominance_h(a, n);
}

bool is_converge(std::vector<double>& x1, std::vector<double>& x2, size_t n, double eps)
{
    double distance = 0;
    for (size_t i = 0; i < n; i++) {
        distance += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
    distance = std::sqrt(distance);

    return distance < eps;
}

std::tuple<std::vector<double>, size_t> seidel(double** a, double* b, size_t n, double eps)
{
    std::vector<double> x(n, 0);
    std::vector<double> prev(n, 0);
    size_t iter = 0;

    do {
        for (size_t i = 0; i < n; i++) {
            prev[i] = x[i];
        }
        for (size_t i = 0; i < n; i++) {
            double temp = 0;
            for (size_t j = 0; j < n; j++) {
                if (j != i)
                    temp += a[i][j] * x[j];
            }
            x[i] = (b[i] - temp) / a[i][i];
        }
        iter++;
    } while (!is_converge(x, prev, n, eps));

    return std::make_tuple(x, iter);
}

int main(int argc, char* argv[])
{
    std::ifstream file(argv[1]);
    double** a;
    double* b;
    size_t n = 0;
    file >> n;

    try {
        a = new double*[n];
        for (size_t i = 0; i < n; i++) {
            a[i] = new double[n];
        }
        b = new double[n];
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
    }

    double num;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            file >> num;
            a[i][j] = num;
        }
        file >> num;
        b[i] = num;
    }

    print_maxtrix(a, b, n);
    if (!is_diagonal_predominance(a, b, n)) {
        std::cout << "The SLAE doesn't converge" << '\n';
        return 0;
    }
    print_maxtrix(a, b, n);

    double eps = 1e-5;
    auto [x, iter] = seidel(a, b, n, eps);

    std::cout << iter << " iterations; eps = " << eps << "\nSolution:\n";
    for (size_t i = 0; i < x.size(); i++) {
        std::cout << x[i] << ' ';
    }
    std::cout << '\n';

    return 0;
}
