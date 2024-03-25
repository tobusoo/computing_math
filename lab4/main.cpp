#include <SFML/Graphics.hpp>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

class L {
private:
    std::vector<double> xi;
    std::vector<double> yi;

    double c(double x, size_t k)
    {
        double temp = 1;
        for (size_t i = 0; i < xi.size(); i++) {
            if (i != k)
                temp *= (x - xi[i]) / (xi[k] - xi[i]);
        }

        return temp;
    }

public:
    L(std::vector<double>& xx, std::vector<double>& yy) : xi(xx), yi(yy)
    {
        assert(xx.size() == yy.size());
    }

    double operator()(double x)
    {
        double y = 0;

        for (size_t i = 0; i < yi.size(); i++) {
            y += yi[i] * c(x, i);
        }

        return y;
    }
};

int main()
{
    double temp;
    std::vector<double> x;
    std::vector<double> y;
    std::fstream file("in.txt");

    for (size_t i = 0; i < 4; i++) {
        file >> temp;
        x.push_back(temp);
    }

    for (size_t i = 0; i < 4; i++) {
        file >> temp;
        y.push_back(temp);
    }

    L l(x, y);
    std::cout << "L(2) = " << l(2) << '\n';

    return 0;
}