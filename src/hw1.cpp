#include "hw1.h"
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using Matrix = std::vector<std::vector<double>>;
//------------------------------zeros----------------------------------------
Matrix algebra::zeros(size_t n, size_t m)
{
    if(n<=0 || m<=0)
        throw std::logic_error("invalied input");
    else
    {
    Matrix mat(n, std::vector<double>(m, 0));
    return mat;
    }
}
//------------------------------ones----------------------------------------
Matrix algebra::ones(size_t n, size_t m)
{
    if(n<=0 || m<=0)
        throw std::logic_error("invalied input");
    else
    {   
    Matrix mat(n, std::vector<double>(m, 1));
    return mat;
    }
}
//------------------------------random---------------------------------------
Matrix  algebra::random(size_t n, size_t m, double min, double max)
{
    if(min>=max)
    {
        throw std::logic_error("invalied input");
    }
    else
    {
        Matrix mat{algebra::zeros(n,m)};
        static std::default_random_engine e{};
        static std::uniform_real_distribution<double> d{min, max};


        for (size_t i {}; i < n; i++)
            for (size_t j {}; j < m; j++)
                mat[i][j]= d(e);
        return mat;
    }  
}
//------------------------------show---------------------------------------
void algebra::show(const Matrix& matrix)
{
    size_t n { matrix.size() };
    size_t m { matrix[0].size() };

    std::cout << std::setprecision(3);
    for (size_t i {}; i < n; ++i) {
        for (size_t j {}; j < m; ++j)
            std::cout << matrix[i][j] << "  ";
        std::cout << std::endl;
    }
}
//------------------------------multiply1---------------------------------------
Matrix algebra::multiply(const Matrix& matrix, double c)
{
    size_t n { matrix.size() };
    size_t m { matrix[0].size() };

    Matrix mat1 { algebra::zeros(n,m) };

    for (size_t i {}; i < n; ++i)
        for (size_t j {}; j < m; ++j)
            mat1[i][j] = matrix[i][j] * c;
    return mat1;
}
//------------------------------multiply2---------------------------------------
Matrix algebra::multiply(const Matrix& matrix1, const Matrix& matrix2)
{
    if(matrix1.size()*matrix2.size() ==0)
    {
        Matrix rslt{};
        return rslt;
    }
    
    size_t n1 { matrix1.size() };
    size_t m1 { matrix1[0].size() };
    size_t n2 { matrix2.size() };
    size_t m2 { matrix2[0].size() };  

    if(m1!=n2 )
    {
        throw std::logic_error("you can't Multiply these matrix");
    }
    else
    {    
    Matrix rslt{algebra::zeros(n1,m2)};
    double tmp{};

    for(size_t i1{};i1<n1;++i1)
        for(size_t j2{};j2<m2;++j2)
        {
            for(size_t i2{};i2<n2;++i2)
                tmp += matrix1[i1][i2]*matrix2[i2][j2];
        rslt[i1][j2] = tmp;
        tmp = 0;
        }
     return rslt;    
    }
}