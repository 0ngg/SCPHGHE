#ifndef SCHEMECFD_H
#define SCHEMECFD_H

#include"temp_util.h"

namespace cfdscheme
{
// minfo
bool is_unique(std::string s)
{
    return !s.empty() && std::find_if(s.begin(), s.end(),
    [](char c) {return std::isdigit(c);}) != s.end(); 
};
int is_neighbor(make<int>::vec v1, make<int>::vec v2)
{
    make<int>::vec::iterator check = std::find_if(v1.begin(), v1.end(), [](make<int>::vec v2,
    int i) {return std::find(v2.begin(), v2.end(), i);});
    if(check != v1.end())
    {
        return *check;
    }
    else
    {
        return -1;
    };
};
};
#endif