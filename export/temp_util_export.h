#include"temp_util.h"
#include"temp_struct_scheme.h"
#include"temp_struct_solver.h"

#ifndef CFDEXPORT_H
#define CFDEXPORT_H

template <typename V>
std::string vec_to_str(make<V>::vec vec_in)
{
    std::string get("'");
    for(auto i = vec_in.begin(); i != vec_in.end(); i++)
    {
        if(*i != vec_in.end())
        {
            get += std::to_string(*i + ",");
        }
        else
        {
            get += std::to_string(*i);
        };
    };
    get += "'";
    return get;
};

#endif