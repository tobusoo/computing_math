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

int rsly_gauss(FractionalNum** a, FractionalNum* x, long long n)
{
    long long imax, i, j, k;
    FractionalNum amax, temp_value;
    std::vector<int> where(n, -1);

    for (k = 0; k < n; k++) {
        imax = k;
        amax = a[k][k].fabs();
        for (i = k + 1; i < n; i++) {
            temp_value = a[i][k].fabs();
            if (temp_value > amax) {
                amax = temp_value;
                imax = i;
            }
        }
        if (k != imax) {
            std::swap(a[k], a[imax]);
            std::swap(x[k], x[imax]);
        }

        if (a[k][k].numerator == 0) {
            continue;
        }
        temp_value = FractionalNum(1) / a[k][k];
        where[k] = k;

        for (i = k; i < n; i++) {
            a[k][i] *= temp_value;
        }
        x[k] *= temp_value;

        for (i = k + 1; i < n; i++) {
            for (j = k + 1; j < n; j++) {
                a[i][j] -= a[i][k] * a[k][j];
            }
            x[i] -= a[i][k] * x[k];

            // Зануление столбцов
            // a[i][k] -= a[i][k] * a[k][k]; // Всегда равно 0
            a[i][k] = 0;
        }

        print_maxtrix(a, x, n);
    }

    bool is_inf_or_zero_sol = false;
    std::vector<FractionalNum> x_before(x, x + n);
    for (i = n - 1; i >= 0; i--) {
        if (where[i] == -1) {
            x[i] = FractionalNum(0);
            is_inf_or_zero_sol = true;
            continue;
        }
        for (j = i + 1; j < n; j++) {
            x[i] -= a[i][j] * x[j];
        }
    }

    if (is_inf_or_zero_sol) {
        for (i = 0; i < n; i++) {
            FractionalNum sum(0);
            for (j = 0; j < n; j++) {
                sum += x[j] * a[i][j];
            }
            sum -= x_before[i];
            temp_value = sum.fabs();
            if (temp_value.numerator != 0)
                return 0;
        }

        for (i = 0; i < n; ++i) {
            if (where[i] == -1)
                return -1; // INFINITY SOLUTION
        }
    }

    return 1;
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

    int solutions = rsly_gauss(a, x, n);
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
