#ifndef LINEARCFD_H
#define LINEARCFD_H

#include"temp_util.h"
#include"temp_struct_scheme.h"

namespace cfdlinear
{
template <class V>
void append_template(make<std::string>::vec which, V* eq, cfdscheme::scheme& scheme_)
{
    // refs
    minfo& mesh_ = scheme_.mesh;
    make<make<double>::sp_mat>::map_str& fc_ = mesh_.fc;
    make<make<double>::sp_mat>::map_str& cc_ = mesh_.cc;
    // iter
    cfdscheme::vinfo value_(which, mesh_);
    make<make<double>::map_int>::map_str gamma_; // cell
    make<make<double>::sp_mat>::map_str lhs_cc_;
    make<make<double>::sp_mat>::map_str lhs_fc_;
    make<make<double>::sp_mat>::map_str rhs_cc_;
    make<make<double>::sp_mat>::map_str rhs_fc_;
    for(auto i = which.begin(); i != which.end(); i++)
    {
        lhs_cc_.insert({*i, make<double>::sp_mat(cc_[*i])});
        rhs_cc_.insert({*i, make<double>::sp_mat(cc_[*i].rows(), 1)});
        if(*i.compare("s2s") != 0)
        {
            gamma_.insert({*i, make<double>::map_int(scheme_.prop.rho["fluid"])});
            lhs_fc_.insert({*i, make<double>::sp_mat(fc_[*i])});
            rhs_fc_.insert({*i, make<double>::sp_mat(fc_[*i])});
        };
    };
    eq->value = value_; 
    eq->lhs_cc = lhs_cc_; eq->rhs_cc = rhs_cc_;
    if(std::find(which.begin(), which.end(), "s2s") == which.end())
    {
        eq->gamma = gamma_; eq->lhs_fc = lhs_fc_; eq->rhs_fc = rhs_fc_;
    };
};
template <class V>
V linear_face(V valC__, V valF__, double gc__, bool is_gamma)
{
    if(is_gamma)
    {
        return (gc__ * valF__ + (1 - gc__) * valC__) / (valC__ * valF__);
    }
    else
    {
        return gc__ * valC__ + (1 - gc__) * valF__;
    };
};
/*
double quick_face()
{

};
coor quick_grad()
{
    
};
*/
double calc_fluid_prop(std::string what, double P, double T)
{

};
};
#endif