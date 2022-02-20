#ifndef AP_HW1_H
#define AP_HW1_H
#include <iostream>
#include <vector>
#include <random>

using Matrix = std::vector<std::vector<double>>;

namespace algebra
{
    Matrix zeros(size_t n, size_t m);
    Matrix ones(size_t n, size_t m);
    Matrix random(size_t n, size_t m, double min, double max);
    void show(const Matrix& matrix);
    Matrix multiply(const Matrix& matrix, double c);
    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2);
}


#endif //AP_HW1_H
