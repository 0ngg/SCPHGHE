#include"temp_util_linear.h"

#ifdef LINEARCFD_H
namespace cfdlinear
{
class linear
{
    public:
    cfdscheme::vinfo value;
    make<make<make<double>::map_int>::map_str>::map_str gamma; // cc
    make<make<double>::sp_mat>::map_str lhs_cc;
    make<make<double>::sp_mat>::map_str lhs_fc;
    make<make<double>::sp_mat>::map_str rhs_cc;
    make<make<double>::sp_mat>::map_str rhs_fc;
    linear() {};
    private:
    virtual void make_linear();
    virtual void update_linear();
    virtual void calc_gamma();
    virtual void calc_lhs();
    virtual void calc_rhs();
    virtual void calc_bound_lhs();
    virtual void calc_bound_rhs();
};
class pcorrect : public linear
{
    public:
    pcorrect() {};
    private:
    void make_linear(cfdscheme::scheme& const);
    void update_linear(cfdscheme::scheme& const, momentum& const, momentum& const, momentum& const);
    void update_correction(cfdscheme::scheme& const, momentum& const, momentum& const, momentum& const);
    void calc_lhs(cfdscheme::scheme& const, momentum& const, momentum& const, momentum& const);
    void calc_rhs(cfdscheme::scheme& const, momentum& const, momentum& const, momentum& const);
    void calc_bound_lhs(int, int, std::string, int, cfdscheme::scheme& const, double);
};
class momentum : public linear
{
    public:
    int axis;
    momentum() {};
    private:
    void make_linear(cfdscheme::scheme& const, int);
    void update_linear(cfdscheme::scheme& const, turb_k& const, momentum& const, momentum& const);
    void calc_gamma(cfdscheme::scheme& const);
    void calc_wall(cfdscheme::scheme&, momentum& const, momentum& const);
    void calc_lhs(cfdscheme::scheme& const, momentum& const, momentum& const);
    void calc_rhs(cfdscheme::scheme& const, turb_k& const, momentum& const, momentum& const);
    void calc_bound_lhs(int, int, std::string, int, cfdscheme::scheme& const, momentum& const, momentum& const);
    void calc_bound_rhs(int, int, std::string, int, cfdscheme::scheme& const, momentum& const, momentum& const);
};
class turb_k : public linear
{
    turb_k() {};
    private:
    void make_linear(cfdscheme::scheme& const);
    void update_linear(cfdscheme::scheme& const, turb_e& const, momentum& const, momentum& const, momentum& const);
    void calc_gamma(cfdscheme::scheme& const);
    void calc_wall(cfdscheme::scheme& const, turb_e& const);
    void calc_lhs(cfdscheme::scheme& const, momentum& const, momentum& const, momentum& const);
    void calc_rhs(cfdscheme::scheme& const, turb_e& const, momentum& const, momentum& const, momentum& const);
    void calc_bound_lhs(int, int, std::string, int, cfdscheme::scheme& const, momentum& const, momentum& const, momentum& const);
    void calc_bound_rhs(int, int, std::string, int, cfdscheme::scheme& const, momentum& const, momentum& const, momentum& const);
};
class turb_e : public linear
{
    turb_e() {};
    private:
    void make_linear(cfdscheme::scheme& const);
    void update_linear(cfdscheme::scheme& const, turb_k& const, momentum& const, momentum& const, momentum& const);
    void calc_gamma(cfdscheme::scheme& const);
    void calc_wall(cfdscheme::scheme& const);
    void calc_lhs(cfdscheme::scheme& const, momentum& const, momentum& const, momentum& const);
    void calc_rhs(cfdscheme::scheme& const, turb_k& const, momentum& const, momentum& const, momentum& const);
    void calc_bound_lhs(int, int, std::string, int, cfdscheme::scheme& const, momentum& const, momentum& const, momentum& const);
    void calc_bound_rhs(int, int, std::string, int, cfdscheme::scheme& const, momentum& const, momentum& const, momentum& const);
};
class energy : public linear
{
    energy() {};
    private:
    // commons
    void make_linear(cfdscheme::scheme& const);
    void update_linear(cfdscheme::scheme& const, momentum& const, momentum& const, momentum& const);
    void calc_gamma(cfdscheme::scheme& const);
    void calc_wall(cfdscheme::scheme& const);
    void calc_lhs(cfdscheme::scheme& const, momentum& const, momentum& const, momentum& const);
    void calc_rhs(cfdscheme::scheme& const, momentum& const, momentum& const, momentum& const);    
    void calc_bound_lhs(int, int, std::string, int, cfdscheme::scheme& const, std::string, momentum& const, momentum& const, momentum& const);
    void calc_bound_rhs(int, int, std::string, int, cfdscheme::scheme& const, std::string, momentum& const, momentum& const, momentum& const);
};
class s2s : public linear
{
    s2s() {};
    private:
    void make_linear(cfdscheme::scheme& const);
    void update_linear(cfdscheme::scheme& const, energy& const);
    void calc_lhs(cfdscheme::scheme& const);
    void calc_rhs(cfdscheme::scheme& const, energy& const);
};
};
#endif