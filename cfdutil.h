#ifndef CFDUTIL_H
#define CFDUTIL_H

#include<iostream>
#include<vector>
#include<map>
#include<string>
#include<utility>
#include<algorithm>
#include<memory>
#include<numeric>
#include<functional>
#include<mshio/mshio.h>
#include<Eigen/Core>
#include<Eigen/Sparse>

using namespace Eigen;

std::unique_ptr<bool> is_true(new bool(true));
std::unique_ptr<bool> is_false(new bool(false));
std::unique_ptr<double> pi(new double(4*std::atan(1.0)));
typedef Eigen::Vector3d coor;
typedef Eigen::Triplet<double> trip;
typedef unsigned int uint;
template<typename U>
struct make
{
    typedef std::vector<U> vec;
    typedef std::vector<std::vector<U>> std_mat;
    typedef std::map<uint, U> map_uint;
    typedef std::map<std::string, U> map_str;
    typedef Eigen::SparseMatrix<U, RowMajor> sp_mat;
};
double square_sum (double x, double y)
{
    return x + y * y;
};
double sqrt_sum(make<double>::vec v)
{
    return pow(std::accumulate(v.begin(), v.end(), 0.0, square_sum), 0.5);
};
struct axes
{
    make<double>::sp_mat x;
    make<double>::sp_mat y;
    make<double>::sp_mat z;
    axes() {};
    axes(make<double>::sp_mat x_in, make<double>::sp_mat y_in, make<double>::sp_mat z_in):
                     x(x_in), y(y_in), z(z_in) {};
    axes(const axes& other)
    {
        this->x = other.x;
        this->y = other.y;
        this->z = other.z;
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
    coor& axes_to_coor(uint row, uint col)
    {
        make<double>::vec ax_vec;
        ax_vec.push_back(this->x.coeffRef(row, col));
        ax_vec.push_back(this->y.coeffRef(row, col));
        ax_vec.push_back(this->z.coeffRef(row, col));
        coor ax_coor(ax_vec.data());
        return ax_coor;
    };
    double& axes_to_val(uint row, uint col)
    {
        make<double>::vec ax_vec;
        ax_vec.push_back(this->x.coeffRef(row, col));
        ax_vec.push_back(this->y.coeffRef(row, col));
        ax_vec.push_back(this->z.coeffRef(row, col));
        double ax_val = sqrt_sum(ax_vec);
        return ax_val;
    };
};

#endif