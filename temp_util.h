#ifndef UTILCFD_H
#define UTILCFD_H

// std lib
#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<string>
#include<utility>
#include<algorithm>
#include<memory>
#include<numeric>
#include<functional>
#include<chrono>
// third party
#include<mshio/mshio.h>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include"HumidAirProp.h"
#include"sqlite3.h"

using namespace Eigen;

std::unique_ptr<bool> is_true(new bool(true));
std::unique_ptr<bool> is_false(new bool(false));
std::unique_ptr<double> pi(new double(4*std::atan(1.0)));
typedef std::chrono::milliseconds milli;
typedef Eigen::Vector3d coor;
typedef Eigen::Triplet<double> sparse_input;
template<typename U>
struct make
{
    typedef std::vector<U> vec;
    typedef std::vector<std::vector<U>> std_mat;
    typedef std::map<int, U> map_int;
    typedef std::map<std::string, U> map_str;
    typedef Eigen::SparseMatrix<U, RowMajor> sp_mat;
    typedef std::pair<std::string, U> unique;
};
double square_sum (double x, double y)
{
    return x + y * y;
};
double sqrt_sum(coor v)
{
    return pow(pow(v(0), 2) + pow(v(1), 2) + pow(v(2), 2), 0.5);
};
struct axes
{
    make<double>::sp_mat x;
    make<double>::sp_mat y;
    make<double>::sp_mat z;
    axes() {};
    axes(make<double>::sp_mat& in): x(in), y(in), z(in) {};
    axes(const axes& other)
    {
        this->x = other.x;
        this->y = other.y;
        this->z = other.z;
    };
    axes& operator()(int row, int col, coor& in)
    {
        this->x.coeffRef(row, col) = in(0);
        this->y.coeffRef(row, col) = in(1);
        this->z.coeffRef(row, col) = in(2);
    };
    axes& operator=(int x)
    {
        return *this;
    };
    axes& operator+=(const axes& other)
    {
        this->x = this->x + other.x;
        this->y = this->y + other.y;
        this->z = this->z + other.z;
        return *this;
    };
    axes& operator+=(const axes& other)
    {
        this->x = this->x + other.x;
        this->x = this->x + other.x;
        this->x = this->x + other.x;
        return *this;
    };
    coor& axes_to_coor(int row, int col)
    {
        make<double>::vec ax_vec;
        ax_vec.push_back(this->x.coeffRef(row, col));
        ax_vec.push_back(this->y.coeffRef(row, col));
        ax_vec.push_back(this->z.coeffRef(row, col));
        coor ax_coor(ax_vec.data());
        return ax_coor;
    };
    double& axes_to_val(int row, int col)
    {
        make<double>::vec ax_vec;
        ax_vec.push_back(this->x.coeffRef(row, col));
        ax_vec.push_back(this->y.coeffRef(row, col));
        ax_vec.push_back(this->z.coeffRef(row, col));
        coor ax_coor(ax_vec.data());
        double ax_val = sqrt_sum(ax_coor);
        return ax_val;
    };
};

#endif