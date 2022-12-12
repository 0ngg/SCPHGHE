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
make<std::pair<std::string, int>>::vec is_boundary(int fid, make<make<make<int>::vec>::map_int>::map_str& fid_, make<std::string>::vec& reserve)
{
    make<std::pair<std::string, int>>::vec get = make<std::pair<std::string, int>>::vec();
    for(std::pair<std::string, make<make<int>::vec>::map_int> entry1 : fid_)
    {
        if(std::find(reserve.begin(), reserve.end(), entry1.first) != reserve.end())
        {
            for(std::pair<int, make<int>::vec> entry2 : fid_[entry1.first])
            {
                for(auto i = entry2.second.begin(); i != entry2.second.end(); i++)
                {
                    if(*i == fid)
                    {
                        std::pair<std::string, int> temp(entry1.first, entry2.first);
                        get.push_back(temp);
                    };
                };
            };
        };
    };
    return get;
};
};
#endif