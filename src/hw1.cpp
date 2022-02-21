#include "hw1.h"
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <cmath>

using Matrix = std::vector<std::vector<double>>;
//------------------------------zeros----------------------------------------
Matrix algebra::zeros(size_t n, size_t m)
{
    if (n <= 0 || m <= 0)
        throw std::logic_error("invalied input");
    else {
        Matrix mat(n, std::vector<double>(m, 0));
        return mat;
    }
}
//------------------------------ones----------------------------------------
Matrix algebra::ones(size_t n, size_t m)
{
    if (n <= 0 || m <= 0)
        throw std::logic_error("invalied input");
    else {
        Matrix mat(n, std::vector<double>(m, 1));
        return mat;
    }
}
//------------------------------random---------------------------------------
Matrix algebra::random(size_t n, size_t m, double min, double max)
{
    if (min >= max) {
        throw std::logic_error("invalied input");
    } else {
        Matrix mat { algebra::zeros(n, m) };
        static std::default_random_engine e {};
        static std::uniform_real_distribution<double> d { min, max };

        for (size_t i {}; i < n; i++)
            for (size_t j {}; j < m; j++)
                mat[i][j] = d(e);
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
            std::cout << std::setw(8)<<matrix[i][j] ;
        std::cout << std::endl;
    }
}
//------------------------------multiply1---------------------------------------
Matrix algebra::multiply(const Matrix& matrix, double c)
{
    if (matrix.size() == 0) {
        Matrix rslt {};
        return rslt;
    }

    size_t n { matrix.size() };
    size_t m { matrix[0].size() };

    Matrix mat1 { algebra::zeros(n, m) };

    for (size_t i {}; i < n; ++i)
        for (size_t j {}; j < m; ++j)
            mat1[i][j] = matrix[i][j] * c;
    return mat1;
}
//------------------------------multiply2---------------------------------------
Matrix algebra::multiply(const Matrix& matrix1, const Matrix& matrix2)
{
    if (matrix1.size() * matrix2.size() == 0) {
        Matrix rslt {};
        return rslt;
    }

    size_t n1 { matrix1.size() };
    size_t m1 { matrix1[0].size() };
    size_t n2 { matrix2.size() };
    size_t m2 { matrix2[0].size() };

    if (m1 != n2) {
        throw std::logic_error("you can't Multiply these matrix");
    } else {
        Matrix rslt { algebra::zeros(n1, m2) };
        double tmp {};

        for (size_t i1 {}; i1 < n1; ++i1)
            for (size_t j2 {}; j2 < m2; ++j2) {
                for (size_t i2 {}; i2 < n2; ++i2)
                    tmp += matrix1[i1][i2] * matrix2[i2][j2];
                rslt[i1][j2] = tmp;
                tmp = 0;
            }
        return rslt;
    }
}
//------------------------------sum1---------------------------------------
Matrix algebra::sum(const Matrix& matrix, double c)
{

    if (matrix.size() == 0) {
        Matrix rslt {};
        return rslt;
    }

    size_t n { matrix.size() };
    size_t m { matrix[0].size() };
    Matrix mat1 { algebra::zeros(n, m) };

    for (size_t i {}; i < n; ++i)
        for (size_t j {}; j < m; ++j)
            mat1[i][j] = matrix[i][j] + c;
    return mat1;
}
//------------------------------sum2---------------------------------------
Matrix algebra::sum(const Matrix& matrix1, const Matrix& matrix2)
{
    if (matrix1.size()==0 && matrix2.size() == 0) {
        Matrix rslt {};
        return rslt;
    }
    if(matrix1.size()==0 || matrix2.size() == 0)
    {
        throw std::logic_error("you can't sum these matrix");
    }

    size_t n1 { matrix1.size() };
    size_t m1 { matrix1[0].size() };
    size_t n2 { matrix2.size() };
    size_t m2 { matrix2[0].size() };
    if (m1 != m2 && n1 != n2) {
        throw std::logic_error("you can't sum these matrix");
    } else {
        Matrix rslt { algebra::zeros(n1, m1) };
        for (size_t i1 {}; i1 < n1; ++i1)
            for (size_t j1 {}; j1 < m1; ++j1)
                rslt[i1][j1] = matrix1[i1][j1] + matrix2[i1][j1];
        return rslt;

    }
}
//------------------------------transpose---------------------------------------
Matrix algebra::transpose(const Matrix& matrix)
{
    if (matrix.size()==0) {
        Matrix rslt {};
        return rslt;
    }

    size_t n { matrix.size() };
    size_t m { matrix[0].size() };
    Matrix rslt {zeros(m,n)};

    for (size_t i {}; i < n; ++i)
        for (size_t j {}; j < m; ++j)
            rslt[j][i] = matrix[i][j];
    return rslt;
}

//------------------------------minor---------------------------------------
Matrix algebra::minor(const Matrix& matrix, size_t n, size_t m)
{   
    if(matrix.size()==0)
    {
        Matrix rslt{};
        return rslt;
    }

    
        //n = n-1;
        //m = m-1;
        size_t row{matrix.size()};
        size_t col{matrix[0].size()};
        Matrix rslt{algebra::zeros(row-1,col-1)};
        size_t ii{};
        size_t jj{};

        for (size_t i {}; i < row; ++i)
            if(i != n)
            {
                for (size_t j {}; j < col; ++j)
                    if(j != m)
                    {
                        rslt[ii][jj] = matrix[i][j];
                        jj++;
                    }
                ii++;
                jj=0;
            }
        return rslt;
}

//------------------------------deteminan---------------------------------------
double algebra::determinant(const Matrix& matrix)
{
    if(matrix.size()==0)
    {
        double rslt{1};
        return rslt;
    }

    size_t n { matrix.size() };
    size_t m { matrix[0].size() };
    if(n==1 && m==1)
    {
        double rslt{matrix[0][0]};
        return rslt;
    }
    if(n!=m)
        throw std::logic_error("this matrix hasn't determinan");
    if(n==2)
    {
        double rslt{matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]};
        return rslt;
    }
    else
    {
        double rslt{};
        for(size_t i{};i<n;++i)
        {
            rslt = rslt + (std::pow(-1,1+i+1))*determinant(algebra::minor(matrix,i,0))*matrix[i][0];
        }
        return rslt;
    }
}
//------------------------------invers---------------------------------------
Matrix algebra::inverse(const Matrix& matrix)
{
    if(matrix.size()==0)
    {
        Matrix rslt{};
        return rslt;
    }
    size_t n { matrix.size() };
    size_t m { matrix[0].size() };
    if(n!=m || algebra::determinant(matrix)==0)
        throw std::logic_error("this matrix hasn't inverse");
    

    Matrix adj{algebra::zeros(n,n)};
    Matrix c{algebra::zeros(n,n)};
    Matrix rslt{};
    for(size_t i{};i<n;i++)
        for(size_t j{};j<n;j++)
            c[i][j]=(std::pow(-1,i+j+2))*algebra::determinant(algebra::minor(matrix,i,j));
    adj = algebra::transpose(c);
    double a{algebra::determinant(matrix)};
    adj = algebra::multiply(adj,1/a);
        return adj;
}
//------------------------------concatenate--------------------------------------------
Matrix algebra::concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis)
{
    if(matrix1.size() == 0)
        return matrix2;
    if(matrix2.size() == 0)
        return matrix1;

    size_t n1 { matrix1.size() };
    size_t m1 { matrix1[0].size() };
    size_t n2 { matrix2.size() };
    size_t m2 { matrix2[0].size() };

    if((axis == 0 && m1 != m2) || (axis == 1 && n1 != n2))
        throw std::logic_error("invalied input.......");
    
    if(axis==0)
    {
        Matrix rslt{algebra::zeros(n1+n2,m1)};

        for(size_t i{};i<n1+n2;i++)
            for(size_t j{};j<m1;j++)
            {   if(i<n1)
                    rslt[i][j] = matrix1[i][j];
                else
                    rslt[i][j] = matrix2[i-n1][j];
            }
        return rslt;
    }
    else
    {
        Matrix rslt{algebra::zeros(n1,m1+m2)};

        for(size_t i{};i<n1;i++)
            for(size_t j{};j<m1+m2;j++)
            {   if(j<m1)
                    rslt[i][j] = matrix1[i][j];
                else
                    rslt[i][j] = matrix2[i][j-m1];
            }
        return rslt;       
    }    
 
}
//------------------------------ero_swap--------------------------------------------
Matrix algebra::ero_swap(const Matrix& matrix, size_t r1, size_t r2)
{
    if(matrix.size() == 0)
        return matrix;

    size_t n { matrix.size() };
    size_t m { matrix[0].size() }; 
    if(r1>n-1 || r2>n-1)
        throw std::logic_error("the row that you entered is larger than the size of matrix");
    
    Matrix rslt{algebra::zeros(n,m)};
    size_t tmp{};
    for(size_t i{};i<n;i++)
    {
        tmp = i;
        if(i==r1)
            tmp=r2;
        else if(i==r2)
            tmp=r1;
        for(size_t j{};j<m;j++)
            rslt[i][j] = matrix[tmp][j];   
    }
    return rslt;
}
//------------------------------ero_multiply--------------------------------------------
Matrix algebra::ero_multiply(const Matrix& matrix, size_t r, double c)
{
    if(matrix.size()==0)
        return matrix;
    size_t n { matrix.size() };
    size_t m { matrix[0].size() }; 
    if(r>n-1)
        throw std::logic_error("the row that you entered is larger than the size of matrix");
    
    Matrix rslt{algebra::zeros(n,m)};
    for(size_t i{};i<n;i++)
        for(size_t j{};j<m;j++)
        {
            if(i == r)
                rslt[i][j] = matrix[i][j]*c;
            else
                rslt[i][j] = matrix[i][j];
        }
    return rslt;
    
}
//------------------------------ero_sum-----------------------------------------------
Matrix algebra::ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2)
{
    if(matrix.size()==0)
        return matrix;
    size_t n { matrix.size() };
    size_t m { matrix[0].size() }; 
    if(r1>n-1 || r2>n-1)
        throw std::logic_error("the row that you entered is larger than the size of matrix");
    
    Matrix rslt{algebra::zeros(n,m)};
    for(size_t i{};i<n;i++)
        for(size_t j{};j<m;j++)
        {
            if(i == r2)
                rslt[i][j] = matrix[r1][j]*c + matrix[r2][j];
            else
                rslt[i][j] = matrix[i][j];
        }
    return rslt;
    
}