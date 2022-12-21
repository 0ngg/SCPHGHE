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
    virtual void update_linear();
    virtual void calc_wall();
    protected:
    virtual void make_linear();
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
    pcorrect(cfdscheme::scheme& const scheme_ref, user& user_ref) {this->make_linear(scheme_ref, user_ref);};
    void update_linear(cfdscheme::scheme& const, momentum, momentum, momentum);
    void update_correction(cfdscheme::scheme& const, momentum, momentum, momentum);
    private:
    void make_linear(cfdscheme::scheme& const, user&);
    void calc_lhs(cfdscheme::scheme& const, momentum, momentum, momentum);
    void calc_rhs(cfdscheme::scheme& const, momentum, momentum, momentum);
    void calc_bound_lhs(int, int, std::string, int, cfdscheme::scheme& const, double);
};
class momentum : public linear
{
    public:
    int axis;
    momentum() {};
    momentum(cfdscheme::scheme& const scheme_ref, int axis_in) {this->make_linear(scheme_ref, axis_in);};
    void update_linear(cfdscheme::scheme& const, turb_k, momentum, momentum);
    void calc_wall(cfdscheme::scheme&, momentum, momentum);
    private:
    void make_linear(cfdscheme::scheme& const, int);
    void calc_gamma(cfdscheme::scheme& const);
    void calc_lhs(cfdscheme::scheme& const, momentum, momentum);
    void calc_rhs(cfdscheme::scheme& const, turb_k, momentum, momentum);
    void calc_bound_lhs(int, int, std::string, int, cfdscheme::scheme& const, momentum, momentum);
    void calc_bound_rhs(int, int, std::string, int, cfdscheme::scheme& const, momentum, momentum);
};
class turb_k : public linear
{
    public:
    turb_k() {};
    turb_k(cfdscheme::scheme& const scheme_ref) {this->make_linear(scheme_ref);};
    void update_linear(cfdscheme::scheme& const, turb_e, momentum, momentum, momentum);
    void calc_wall(cfdscheme::scheme& const, turb_e);
    private:
    void make_linear(cfdscheme::scheme& const);
    void calc_gamma(cfdscheme::scheme& const);
    void calc_lhs(cfdscheme::scheme& const, momentum, momentum, momentum);
    void calc_rhs(cfdscheme::scheme& const, turb_e, momentum, momentum, momentum);
    void calc_bound_lhs(int, int, std::string, int, cfdscheme::scheme& const, momentum, momentum, momentum);
    void calc_bound_rhs(int, int, std::string, int, cfdscheme::scheme& const, momentum, momentum, momentum);
};
class turb_e : public linear
{
    public:
    turb_e() {};
    turb_e(cfdscheme::scheme& const scheme_ref) {this->make_linear(scheme_ref);};
    void update_linear(cfdscheme::scheme& const, turb_k, momentum, momentum, momentum);
    void calc_wall(cfdscheme::scheme& const);
    private:
    void make_linear(cfdscheme::scheme& const);
    void calc_gamma(cfdscheme::scheme& const);
    void calc_lhs(cfdscheme::scheme& const, momentum, momentum, momentum);
    void calc_rhs(cfdscheme::scheme& const, turb_k, momentum, momentum, momentum);
    void calc_bound_lhs(int, int, std::string, int, cfdscheme::scheme& const, momentum, momentum, momentum);
    void calc_bound_rhs(int, int, std::string, int, cfdscheme::scheme& const, momentum, momentum, momentum);
};
class energy : public linear
{
    public:
    energy() {};
    energy(cfdscheme::scheme& const scheme_ref, user& user_ref) {this->make_linear(scheme_ref, user_ref);};
    void update_linear(cfdscheme::scheme& const, turb_k, momentum, momentum, momentum, double);
    void calc_wall(cfdscheme::scheme& const, double);
    private:
    // commons
    void make_linear(cfdscheme::scheme& const, user&);
    void calc_gamma(cfdscheme::scheme& const, double);
    void calc_lhs(cfdscheme::scheme& const, turb_k, momentum, momentum, momentum, double);
    void calc_rhs(cfdscheme::scheme& const, momentum, momentum, momentum, double);    
    void calc_bound_lhs(int, int, std::string, int, cfdscheme::scheme& const, std::string, momentum, momentum, momentum, double);
    void calc_bound_rhs(int, int, std::string, int, cfdscheme::scheme& const, std::string, momentum, momentum, momentum, double);
};
class s2s : public linear
{
    public:
    s2s() {};
    s2s(cfdscheme::scheme& const scheme_ref) {this->make_linear(scheme_ref);};
    void update_linear(cfdscheme::scheme& const, energy);
    void update_source_s2s(cfdscheme::scheme& const);
    private:
    void make_linear(cfdscheme::scheme& const);
    void calc_lhs(cfdscheme::scheme& const);
    void calc_rhs(cfdscheme::scheme& const, energy);
};
};
#endif