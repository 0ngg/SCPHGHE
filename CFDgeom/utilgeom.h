/*
utilgeom
    1|	Surface area and volume computation
    2|  Check-in-boundary
    3|  Match vector elements
    4| 
*/

#ifndef UTILGEOM_H
#define UTILGEOM_H

#include<iostream>
#include<vector>
#include<map>
#include<string>
#include<utility>
#include<cmath>
#include<algorithm>
#include<numeric>
#include<functional>
#include<Eigen/Dense>

typedef unsigned int uint;
typedef Eigen::Vector3d coor;

namespace cfdgeom
{

std::string check_tag(std::map<uint, std::string> tag_map, uint id)
{
    if(std::find(tag_map.begin(), tag_map.end(), id) != tag_map.end())
    {
        return tag_map[id];
    };

    return "none";
}; // check_tag

bool match_vec(std::vector<uint> vec1, std::vector<uint> vec2)
{
    uint target = vec1.size();
    for(auto i = vec1.begin(); i != vec1.end(); i++)
    {
        {
        if(ctd != target && *i == *j)
        {
            ctd++;
        };
        else if(ctd == target)
        {
            return true;
        };
        }; // if else
    }; // for

    return false;
}; // match_vec

std::map<uint, std::pair<coor, coor>> find_diagonals(std::vector<coor> lcoor)
{
    std::pair<coor, coor> diag_1;
    std::pair<coor, coor> diag_2;
    float length = 0;
    coor vec << 0.0d 0.0d 0.0d;

    for(uint i = 0; i < sizeof(lcoor)-1; i++)
    {
        for(uint j = i; j < sizeof(lcoor); j++)
        {
            // get length
            for(uint k = 0; k < 3; k++)
            {
                vec[k] = *i[k] + *j[k];
                vec[k] = pow(vec[k], 2);
            };

            double rhs = std::accumulate(vec.begin(), vec.end(), 0.0);
            rhs = pow(rhs, 0.5);

            if(rhs > length)
            {
                length = rhs;
                diag_1.first = lcoor[i];
                diag_1.second = lcoor[j];

                std::vector<uint> temp_d2;
                for(uint k = 0; k < sizeof(lcoor); k++)
                {
                    if(k != diag_1.first && k != diag_1.second)
                    {
                        temp_d2.push_back(k);
                    };

                    diag_2.first = lcoor(temp_d2[0]);
                    diag_2.second = lcoor(temp_d2[1]);
                };
            }; // if
        }; // for
    }; // for

    std::map<uint, std::pair<coor, coor>> get;
    get.insert({0, diag_1});
    get.insert({1, diag_2});

    return get;
}; // find_diagonals

std::pair<std::pair<coor, coor>, double> info_2d(std::vector<coor> lcoor)
{
   // triangle and quadrilaterals only for now
    std::pair<std::pair<coor, coor>, double> get;

    switch(sizeof(lcoor))
    {   
        case 3:
        coor tcentroid = std::accumulate(lcoor.begin(), lcoor.end(), Eigen::Vector3d::Zero());
        tcentroid = tcentroid/3;

        // area = 0.5 * |ab x ac|
        double tarea = 0.5 * ((lcoor[1] - lcoor[0]).cross(lcoor[2] - lcoor[0])).norm();
        coor tnormal = (lcoor[1] - lcoor[0]).cross(lcoor[2] - lcoor[0]) / (lcoor[1] - lcoor[0]).cross(lcoor[2] - lcoor[0]).norm();

        case 4:
        std::map<uint, std::pair<coor, coor>> diag = find_diagonals(lcoor);

        // centroid = sum of all (tri centroid * area) / (sum of area)
        // area = sum of all triangle area
        coor c1 = ( diag[0].first + diag[0].second + diag[1].first  ) / 3;
        coor c2 = ( diag[0].first + diag[0].second + diag[1].second  ) / 3
        double a1 = 0.5 * ((diag[1].first - diag[0].first).cross(diag[1].first - diag[0].second)).norm();
        double a2 = 0.5 * ((diag[1].second - diag[0].first).cross(diag[1].second - diag[0].second)).norm();

        coor tcentroid = ( c1*a1 + c2*a2 ) / (a1 + a2);
        double tarea = a1 + a2;

        coor tnormal =  (diag[0].first - diag[1].first).cross(diag[0].second - diag[1].first) /
                        (diag[0].first - diag[1].first).cross(diag[0].second - diag[1].first).norm();

        default:
        cout << "shape not supported. Abort.." << endl;
        return;
    }; // switch

    std::pair<coor, coor> arnorm;
    arnorm.first = tcentroid;
    arnorm.second = tnormal;

    get.first = arnorm;
    get.second = tarea;

    return get
}; // info_2d

std::pair<coor, double> info_3d(std::vector<uint> lface, std::map<uint, finfo>* aface)
{
    std::pair<coor, double> get;
    coor tcentroid = Eigen::Vector3d::Zero();
    double ctd_area = 0;

    // face-weighted centroid
    for(auto i = lface.begin(); i != lface.end(); i++)
    {
        tcentroid = tcentroid + (*aface[*i]->fcentroid * (*aface[*i]->farea));
        ctd_area += *aface[*i]->farea;
    };

    tcentroid = tcentroid/ctd_area;

    // volume computation
    double tvol = 0;
    for(auto i = lface.begin(); i != lface.end(); i++)
    {
        coor cf = *aface[*i] - centroid;
        {
        if(*aface[*i]->fnormal.dot(cf) >= 0)
        {
            coor nf = *aface[*i];
        };
        else:
        {
            coor nf = *aface[*i] * -1;
        };
        };

        tvol = tvol + (cf.dot(nf)) * *aface[*i]->farea;
    };

    tvol = tvol/3;

    get.first = tcentroid;
    get.second = tvol;

    return get
}; // info_3d

}; // namespace cfdgeom

#endif