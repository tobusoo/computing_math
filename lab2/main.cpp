#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>

struct FractionalNum {
    //   numerator
    // -------------
    //  denominator

    long long numerator;
    long long denominator;

    FractionalNum(long long numerator = 1, long long denominator = 1)
        : numerator(numerator), denominator(denominator)
    {
        assert(denominator != 0);
        to_short();
    }

    void check_sign()
    {
        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }
    }

    void mul(const FractionalNum& b)
    {
        numerator *= b.numerator;
        denominator *= b.denominator;
        to_short();
    }

    FractionalNum operator*(const FractionalNum& other)
    {
        FractionalNum temp(numerator, denominator);
        temp.mul(other);

        return temp;
    }

    FractionalNum& operator*=(const FractionalNum& other)
    {
        this->mul(other);

        return *this;
    }

    void div(const FractionalNum& b)
    {
        if (b.numerator == 0)
            throw std::runtime_error("division by zero");
        numerator *= b.denominator;
        denominator *= b.numerator;
        to_short();
    }

    FractionalNum operator/(const FractionalNum& other)
    {
        FractionalNum temp(numerator, denominator);
        temp.div(other);

        return temp;
    }

    FractionalNum& operator/=(const FractionalNum& other)
    {
        this->div(other);

        return *this;
    }

    void sum(const FractionalNum& b)
    {
        numerator = numerator * b.denominator + b.numerator * denominator;
        denominator = denominator * b.denominator;
        to_short();
    }

    FractionalNum operator+(const FractionalNum& other)
    {
        FractionalNum temp(numerator, denominator);
        temp.sum(other);

        return temp;
    }

    FractionalNum& operator+=(const FractionalNum& other)
    {
        this->sum(other);

        return *this;
    }

    void sub(const FractionalNum& b)
    {
        numerator = numerator * b.denominator - b.numerator * denominator;
        denominator = denominator * b.denominator;
        to_short();
    }

    FractionalNum operator-(const FractionalNum& other)
    {
        FractionalNum temp(numerator, denominator);
        temp.sub(other);

        return temp;
    }

    FractionalNum& operator-=(const FractionalNum& other)
    {
        this->sub(other);

        return *this;
    }

    FractionalNum fabs()
    {
        return FractionalNum(std::fabs(numerator), std::fabs(denominator));
    }

    void to_short()
    {
        assert(denominator != 0);
        if (numerator == 0) {
            denominator = 1;
            return;
        }

        long long nod = std::gcd(numerator, denominator);
        numerator /= nod;
        denominator /= nod;
        check_sign();
    }

    FractionalNum& operator=(const FractionalNum& other)
    {
        if (this != &other) {
            numerator = other.numerator;
            denominator = other.denominator;
        }

        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, FractionalNum& num)
    {
        num.to_short();

        if (num.denominator == 1)
            os << num.numerator;
        else
            os << num.numerator << '/' << num.denominator;

        return os;
    }

    bool operator>(const FractionalNum& b) const
    {
        double temp = numerator / denominator;
        double temp2 = b.numerator / b.denominator;
        return temp > temp2;
    }

    bool operator<(const FractionalNum& b) const
    {
        return !(*this > b);
    }

    friend std::istream& operator>>(std::istream& os, FractionalNum& num)
    {
        long long numer = 1, denomin = 1;
        char ch = 0;

        os >> numer >> ch;
        if (ch != '/') {
            num.numerator = numer;
            num.denominator = 1;
            os.putback(ch);
        } else {
            os >> denomin;
            if (denomin == 0) {
                throw std::runtime_error("Could not read the number");
            }
            num.numerator = numer;
            num.denominator = denomin;
        }

        num.to_short();
        return os;
    }
};

void print_maxtrix(FractionalNum** a, FractionalNum* x, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n - 1; j++) {
            // std::cout << a[i][j] << "x" << j + 1 << " + ";
            std::cout << a[i][j] << ' ';
        }
        std::cout << a[i][n - 1] << ' ';
        // std::cout << a[i][n - 1] << "x" << n << " = ";
        std::cout << x[i] << '\n';
    }
    std::cout << '\n';
}

int seidel(FractionalNum** a, FractionalNum* x, long long n)
{
    return 0;
}

int main(int argc, char* argv[])
{
    std::ifstream file(argv[1]);
    FractionalNum** a;
    FractionalNum* x;
    size_t n = 0;
    file >> n;
    try {
        a = new FractionalNum*[n];
        for (size_t i = 0; i < n; i++) {
            a[i] = new FractionalNum[n];
        }
        x = new FractionalNum[n];
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
    }

    FractionalNum num;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            file >> num;
            a[i][j] = num;
        }
        file >> num;
        x[i] = num;
    }

    print_maxtrix(a, x, n);

    int solutions = seidel(a, x, n);
    if (solutions == 1) {
        for (size_t i = 0; i < n; i++) {
            std::cout << x[i] << ' ';
        }
        std::cout << '\n';
    } else if (solutions == -1)
        printf("Infinity solutions\n");
    else
        printf("Zero solutions\n");

    return 0;
}