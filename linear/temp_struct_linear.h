#include"temp_util.h"
#include"temp_util_linear.h"

#ifdef LINEARCFD_H
namespace cfdlinear
{
class linear
{
    public:
    cfdscheme::vinfo value;
    make<make<double>::sp_mat>::map_str gamma;
    make<make<double>::sp_mat>::map_str lhs_cc;
    make<make<double>::sp_mat>::map_str lhs_fc;
    make<make<double>::sp_mat>::map_str rhs_cc;
    make<make<double>::sp_mat>::map_str rhs_fc;
    linear() {};
    private:
    virtual void make_linear();
    virtual void calc_source();
    virtual void calc_lhs_cc();
    virtual void calc_lhs_fc();
    virtual void calc_rhs_cc();
    virtual void calc_rhs_fc();
    virtual void calc_wall();
};
class pcorrect : public linear
{
    public:
    pcorrect() {};
    private:
    void make_linear(cfdscheme::scheme&);
    void calc_source(cfdscheme::scheme&);
    void calc_lhs_cc(cfdscheme::scheme&);
    void calc_lhs_fc(cfdscheme::scheme&);
    void calc_rhs_cc(cfdscheme::scheme&);
    void calc_rhs_fc(cfdscheme::scheme&);
};
class momentum : public linear
{
    public:
    int axis;
    momentum() {};
    private:
    void make_linear(cfdscheme::scheme&);
    void calc_source(cfdscheme::scheme&);
    void calc_lhs_cc(cfdscheme::scheme&);
    void calc_lhs_fc(cfdscheme::scheme&);
    void calc_rhs_cc(cfdscheme::scheme&);
    void calc_rhs_fc(cfdscheme::scheme&);
    void calc_wall(cfdscheme::scheme&);
};
class energy : public linear
{
    energy() {};
    private:
    void make_linear(cfdscheme::scheme&);
    void calc_source(cfdscheme::scheme&);
    void calc_lhs_cc(cfdscheme::scheme&);
    void calc_lhs_fc(cfdscheme::scheme&);
    void calc_rhs_cc(cfdscheme::scheme&);
    void calc_rhs_fc(cfdscheme::scheme&);
    void calc_wall(cfdscheme::scheme&);
};
class turb_k : public linear
{
    turb_k() {};
    private:
    void make_linear(cfdscheme::scheme&);
    void calc_source(cfdscheme::scheme&);
    void calc_lhs_cc(cfdscheme::scheme&);
    void calc_lhs_fc(cfdscheme::scheme&);
    void calc_rhs_cc(cfdscheme::scheme&);
    void calc_rhs_fc(cfdscheme::scheme&);
    void calc_wall(cfdscheme::scheme&);
};
class turb_e : public linear
{
    turb_e() {};
    private:
    void make_linear(cfdscheme::scheme&);
    void calc_source(cfdscheme::scheme&);
    void calc_lhs_cc(cfdscheme::scheme&);
    void calc_lhs_fc(cfdscheme::scheme&);
    void calc_rhs_cc(cfdscheme::scheme&);
    void calc_rhs_fc(cfdscheme::scheme&);
    void calc_wall(cfdscheme::scheme&);
};
class s2s : public linear
{
    s2s() {};
    private:
    void make_linear(cfdscheme::scheme&);
    void calc_lhs_cc(cfdscheme::scheme&);
    void calc_rhs_cc(cfdscheme::scheme&);
};
};
#endif