#include"temp_util.h"
#include"temp_util_linear.h"
#include"temp_struct_linear.h"

#ifdef LINEARCFD_H
namespace cfdlinear
{
// pcorrect
void pcorrect::make_linear(cfdscheme::scheme& scheme_)
{
    // fluid only
    append_template(make<std::string>::vec{"fluid"}, this, scheme_);
};
void pcorrect::calc_source(cfdscheme::scheme& scheme_)
{
    // inlet fluid lhs_fc, outlet fluid lhs_fc
    
};
void pcorrect::calc_lhs_cc(cfdscheme::scheme& scheme_)
{

};
void pcorrect::calc_lhs_fc(cfdscheme::scheme&)
{

};
void pcorrect::calc_rhs_cc(cfdscheme::scheme&)
{

};
void pcorrect::calc_rhs_fc(cfdscheme::scheme&)
{

};
// momentum
void momentum::make_linear(cfdscheme::scheme& scheme_)
{
    // fluid only
    append_template(make<std::string>::vec{"fluid"}, this, scheme_);
};
void momentum::calc_source(cfdscheme::scheme&)
{
    // ns fluid lhs_fc, ns fluid rhs_fc
    // inlet fluid lhs_fc, outlet fluid lhs_fc
};
void momentum::calc_lhs_cc(cfdscheme::scheme&)
{

};
void momentum::calc_lhs_fc(cfdscheme::scheme&)
{

};
void momentum::calc_rhs_cc(cfdscheme::scheme&)
{

};
void momentum::calc_rhs_fc(cfdscheme::scheme&)
{

};
void momentum::calc_wall(cfdscheme::scheme&)
{

};
// energy
void energy::make_linear(cfdscheme::scheme& scheme_)
{
    // all but s2s
    // refs
    make<make<double>::sp_mat>::map_str& cc_ = scheme_.mesh.cc;
    make<std::string>::vec append_arg;
    for(std::pair<std::string, make<double>::sp_mat> entry : cc_)
    {
        if(entry.first.compare("s2s") != 0)
        {
            append_arg.push_back(entry.first);
        };
    };
    append_template(append_arg, this, scheme_);
};
void energy::calc_source(cfdscheme::scheme&)
{
    // temp temp(unique)/fluid-solid lhs_fc, temp temp(unique)/fluid-solid rhs_cc
    // irr irr(unique)/fluid-solid rhs_fc, s2s s2s(unique)/fluid-solid rhs_fc
    // hamb lhs_fc, hamb rhs_fc
    
};
void energy::calc_lhs_cc(cfdscheme::scheme&)
{

};
void energy::calc_lhs_fc(cfdscheme::scheme&)
{

};
void energy::calc_rhs_cc(cfdscheme::scheme&)
{

};
void energy::calc_rhs_fc(cfdscheme::scheme&)
{

};
void energy::calc_wall(cfdscheme::scheme&)
{

};
// turb_k
void turb_k::make_linear(cfdscheme::scheme& scheme_)
{
    // fluid only
    append_template(make<std::string>::vec{"fluid"}, this, scheme_);
};
void turb_k::calc_source(cfdscheme::scheme&)
{
    // inlet lhs_fc, inlet rhs_fc
    // outlet lhs_fc, outlet rhs_fc
};
void turb_k::calc_lhs_cc(cfdscheme::scheme&)
{

};
void turb_k::calc_lhs_fc(cfdscheme::scheme&)
{

};
void turb_k::calc_rhs_cc(cfdscheme::scheme&)
{

};
void turb_k::calc_rhs_fc(cfdscheme::scheme&)
{

};
void turb_k::calc_wall(cfdscheme::scheme&)
{

};
// turb_e
void turb_e::make_linear(cfdscheme::scheme& scheme_)
{
    // fluid only
    append_template(make<std::string>::vec{"fluid"}, this, scheme_);
};
void turb_e::calc_source(cfdscheme::scheme&)
{
    // inlet lhs_fc, inlet rhs_fc
    // outlet lhs_fc, outlet rhs_fc
};
void turb_e::calc_lhs_cc(cfdscheme::scheme&)
{

};
void turb_e::calc_lhs_fc(cfdscheme::scheme&)
{

};
void turb_e::calc_rhs_cc(cfdscheme::scheme&)
{

};
void turb_e::calc_rhs_fc(cfdscheme::scheme&)
{

};
void turb_e::calc_wall(cfdscheme::scheme&)
{

};
// s2s
void s2s::make_linear(cfdscheme::scheme& scheme_)
{
    // s2s only
    append_template(make<std::string>::vec{"s2s"}, this, scheme_);
};
void s2s::calc_lhs_cc(cfdscheme::scheme&)
{

};
void s2s::calc_rhs_cc(cfdscheme::scheme&)
{

};
};
#endif