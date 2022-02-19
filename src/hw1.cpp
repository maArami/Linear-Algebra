#include "hw1.h"
#include <iostream>
#include <vector>

using Matrix = std::vector<std::vector<double>>;

Matrix algebra::zeros(size_t n, size_t m)
{

    Matrix mat(n , std::vector<double> (m,0));
    
    return mat;


}


