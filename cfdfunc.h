#ifndef CFDFUNC_H
#define CFDFUNC_H

#include"cfdutil.h"
#include"cfdstruct.h"

namespace cfd
{
// linear
class linear
{
    public:
    linear() {};
    private:
    virtual double a_center();
    virtual double a_face();
    virtual double b_center();
    virtual double no_slip();
    virtual double inlet();
    virtual double outlet();
    virtual double wall();
};
class momentum : public linear
{
    std::string axis;
    make<double>::sp_mat g_src_p;
    make<double>::sp_mat body_src;
    momentum(std::string& axis_in, make<double>::sp_mat& g_src_p_in,
             make<double>::sp_mat& body_src_in): axis(axis_in), g_src_p(g_src_p_in)
    {   
        std::string check("y");
        if(!axis_in.compare(check) == 0)
        {
            this->body_src = body_src_in;
        };
    };
    private:
    double a_center(make<double>::sp_mat&, make<double>::sp_mat&, make<double>::sp_mat&,
                    make<double>::sp_mat&, uint&);
    double a_face(make<double>::sp_mat&, make<double>::sp_mat&, make<uint>::sp_mat&,
                  make<double>::sp_mat&, make<double>::sp_mat&, uint&, uint&);
    double b_center(std::string&, make<double>::sp_mat&, make<double>::sp_mat&,
                    make<double>::sp_mat&, make<double>::sp_mat&, make<uint>::sp_mat&,
                    make<double>::sp_mat&, make<double>::sp_mat&, uint&);
    double no_slip();
    double inlet();
    double outlet();
    double wall();
};
double momentum::a_center(make<double>::sp_mat& g_diff_r, make<double>::sp_mat& g_conv_r,
                          make<double>::sp_mat& rho_v_sf_r, make<double>::sp_mat& prop_f_r,
                          uint& row_r)
{
    double coef_diff = 0.0;
    double coef_conv = 0.0;
    double g_diff_t;
    double g_conv_t;
    double rho_v_sf_t;
    double prop_f_t;
    for(SparseMatrix<double, RowMajor>::InnerIterator it(g_diff_r, row_r); it; it++)
    {
        uint col_t = it.col();
        g_diff_t = g_diff_r.coeffRef(row_r, col_t);
        g_conv_t = g_conv_r.coeffRef(row_r, col_t);
        rho_v_sf_t = rho_v_sf_r.coeffRef(row_r, col_t);
        prop_f_t = prop_f_r.coeffRef(row_r, col_t);
        coef_diff += g_diff_t * (prop_f_t);
        coef_conv += (1 - g_conv_t) * (rho_v_sf_t);
    };
    return (coef_diff + coef_conv);
};
double momentum::a_face(make<double>::sp_mat& g_diff_r, make<double>::sp_mat& g_conv_r,
                        make<uint>::sp_mat& cc_fc_r, make<double>::sp_mat& rho_v_sf_r,
                        make<double>::sp_mat& prop_f_r, uint& row_r, uint& col_r)
{
    uint col_t = cc_fc_r.coeffRef(row_r, col_r);
    double g_diff_t = g_diff_r.coeffRef(row_r, col_t);
    double g_conv_t = g_conv_r.coeffRef(row_r, col_t);
    double rho_v_sf_t = rho_v_sf_r.coeffRef(row_r, col_t);
    double prop_f_t = prop_f_r.coeffRef(row_r, col_t);
    double coef_diff = (-1) * g_diff_t * prop_f_t;
    double coef_conv = g_conv_t * rho_v_sf_t;
    return (coef_diff + coef_conv);
};
double momentum::b_center(std::string& axis_r, make<double>::sp_mat& g_src_diff_r,
                          make<double>::sp_mat& g_src_conv_r, make<double>::sp_mat& g_src_p_r,
                          make<double>::sp_mat& body_src_r, make<uint>::sp_mat& cc_fc_r,
                          make<double>::sp_mat& rho_v_sf_r, make<double>::sp_mat& prop_f_r,
                          uint& row_r)
{
    std::string check("y");
    double src_diff = 0.0;
    double src_conv = 0.0;
    double src_p = 0.0;
    uint col_r;
    uint col_t;
    double rho_v_sf_t;
    double prop_f_t;
    double g_diff_t;
    double g_conv_t;
    double g_p_t;
    double body_t = body_src_r.coeffRef(row_r, row_r);
    for(SparseMatrix<uint, RowMajor>::InnerIterator it(cc_fc_r, row_r); it; it++)
    {
        col_r = it.col();
        col_t = cc_fc_r.coeffRef(row_r, col_r);
        rho_v_sf_t = rho_v_sf_r.coeffRef(row_r, col_t);
        prop_f_t = prop_f_r.coeffRef(row_r, col_t);
        g_diff_t = g_src_diff_r.coeffRef(row_r, col_t);
        g_conv_t = g_src_conv_r.coeffRef(row_r, col_t);
        g_p_t = g_src_p_r.coeffRef(row_r, col_t);
        src_diff += g_diff_t * prop_f_t;
        src_conv += g_conv_t * rho_v_sf_t;
        src_p += g_p_t;
    };
    if(axis_r.compare(check) == 0)
    {
        return (src_diff + src_conv - src_p + body_t);
    }
    else
    {
        return (src_diff + src_conv - src_p);
    };
};
class pcorrect : public linear
{
    private:
    double a_center(make<double>::sp_mat&, uint&);
    double a_face(make<double>::sp_mat&, make<uint>::sp_mat&, uint&, uint&);
    double b_center(make<double>::sp_mat&, make<uint>::sp_mat&, uint&);
    double no_slip();
    double inlet();
    double outlet();
};
double pcorrect::a_center(make<double>::sp_mat& g_r, uint& row_r)
{
    double coef = 0.0;
    double g_t;
    for(SparseMatrix<double, RowMajor>::InnerIterator it(g_r, row_r); it; it++)
    {
        uint col_t = it.col();
        g_t = g_r.coeffRef(row_r, col_t);
        coef += g_t;
    };
    return ((-1) * coef);
};
double pcorrect::a_face(make<double>::sp_mat& g_r, make<uint>::sp_mat& cc_fc_r,
                         uint& row_r, uint& col_r)
{
    uint col_t = cc_fc_r.coeffRef(row_r, col_r);
    double g_t = g_r.coeffRef(row_r, col_t);
    return g_t;
};
double pcorrect::b_center(make<double>::sp_mat& g_src_r, make<uint>::sp_mat& cc_fc_r,
                          uint& row_r)
{
    double src = 0.0;
    uint col_r;
    uint col_t;
    double g_t;
    for(SparseMatrix<uint, RowMajor>::InnerIterator it(cc_fc_r, row_r); it; it++)
    {
        col_r = it.col();
        col_t = cc_fc_r.coeffRef(row_r, col_r);
        g_t = g_src_r.coeffRef(row_r, col_t);
        src += g_t;
    };
    return src;
};
class energy : public linear
{
    private:
    double a_center(make<double>::sp_mat&, make<double>::sp_mat&, make<double>::sp_mat&,
                    make<double>::sp_mat&, uint&);
    double a_face(make<double>::sp_mat&, make<double>::sp_mat&, make<uint>::sp_mat&,
                  make<double>::sp_mat&, make<double>::sp_mat&, uint&, uint&);
    double b_center(make<double>::sp_mat&, make<double>::sp_mat&, make<double>::sp_mat&,
                    make<uint>::sp_mat&, make<double>::sp_mat&, make<double>::sp_mat&,
                    uint&);
    double isotherm();
    double isoflux();
    double conv_amb();
    double conv_conj();
};
double energy::a_center(make<double>::sp_mat& g_diff_r, make<double>::sp_mat& g_conv_r,
                        make<double>::sp_mat& rho_v_sf_r, make<double>::sp_mat& prop_f_r,
                        uint& row_r)
{
    double coef_diff = 0.0;
    double coef_conv = 0.0;
    double g_diff_t;
    double g_conv_t;
    double rho_v_sf_t;
    double prop_f_t;
    for(SparseMatrix<double, RowMajor>::InnerIterator it(g_diff_r, row_r); it; it++)
    {
        uint col_t = it.col();
        g_diff_t = g_diff_r.coeffRef(row_r, col_t);
        g_conv_t = g_conv_r.coeffRef(row_r, col_t);
        rho_v_sf_t = rho_v_sf_r.coeffRef(row_r, col_t);
        prop_f_t = prop_f_r.coeffRef(row_r, col_t);
        coef_diff += g_diff_t * (prop_f_t);
        coef_conv += (1 - g_conv_t) * (rho_v_sf_t);
    };
    return (coef_diff + coef_conv);
};
double energy::a_face(make<double>::sp_mat& g_diff_r, make<double>::sp_mat& g_conv_r,
                      make<uint>::sp_mat& cc_fc_r, make<double>::sp_mat& rho_v_sf_r,
                      make<double>::sp_mat& prop_f_r, uint& row_r, uint& col_r)
{
    uint col_t = cc_fc_r.coeffRef(row_r, col_r);
    double g_diff_t = g_diff_r.coeffRef(row_r, col_t);
    double g_conv_t = g_conv_r.coeffRef(row_r, col_t);
    double rho_v_sf_t = rho_v_sf_r.coeffRef(row_r, col_t);
    double prop_f_t = prop_f_r.coeffRef(row_r, col_t);
    double coef_diff = (-1) * g_diff_t * prop_f_t;
    double coef_conv = g_conv_t * rho_v_sf_t;
    return (coef_diff + coef_conv);
};
double energy::b_center(make<double>::sp_mat& g_src_diff_r, make<double>::sp_mat& g_src_conv_r,
                        make<double>::sp_mat& body_src_r, make<uint>::sp_mat& cc_fc_r,
                        make<double>::sp_mat& rho_v_sf_r, make<double>::sp_mat& prop_f_r,
                        uint& row_r)
{
    double src_diff = 0.0;
    double src_conv = 0.0;
    uint col_r;
    uint col_t;
    double rho_v_sf_t;
    double prop_f_t;
    double g_diff_t;
    double g_conv_t;
    double body_t = body_src_r.coeffRef(row_r, row_r);
    for(SparseMatrix<uint, RowMajor>::InnerIterator it(cc_fc_r, row_r); it; it++)
    {
        col_r = it.col();
        col_t = cc_fc_r.coeffRef(row_r, col_r);
        rho_v_sf_t = rho_v_sf_r.coeffRef(row_r, col_t);
        prop_f_t = prop_f_r.coeffRef(row_r, col_t);
        g_diff_t = g_src_diff_r.coeffRef(row_r, col_t);
        g_conv_t = g_src_conv_r.coeffRef(row_r, col_t);
        src_diff += g_diff_t * prop_f_t;
        src_conv += g_conv_t * rho_v_sf_t;
    };
    return (src_diff + src_conv + body_t);
};
class turb_k : public linear
{
    private:
    double a_center(make<double>::sp_mat&, make<double>::sp_mat&, make<double>::sp_mat&,
                    make<double>::sp_mat&, uint&);
    double a_face(make<double>::sp_mat&, make<double>::sp_mat&, make<uint>::sp_mat&,
                  make<double>::sp_mat&, make<double>::sp_mat&, uint&, uint&);
    double b_center(make<double>::sp_mat&, make<double>::sp_mat&, make<double>::sp_mat&,
                    make<uint>::sp_mat&, make<double>::sp_mat&, make<double>::sp_mat&,
                    uint&);
    double inlet();
    double outlet();
    double wall();
};
double turb_k::a_center(make<double>::sp_mat& g_diff_r, make<double>::sp_mat& g_conv_r,
                        make<double>::sp_mat& rho_v_sf_r, make<double>::sp_mat& prop_f_r,
                        uint& row_r)
{
    double coef_diff = 0.0;
    double coef_conv = 0.0;
    double g_diff_t;
    double g_conv_t;
    double rho_v_sf_t;
    double prop_f_t;
    for(SparseMatrix<double, RowMajor>::InnerIterator it(g_diff_r, row_r); it; it++)
    {
        uint col_t = it.col();
        g_diff_t = g_diff_r.coeffRef(row_r, col_t);
        g_conv_t = g_conv_r.coeffRef(row_r, col_t);
        rho_v_sf_t = rho_v_sf_r.coeffRef(row_r, col_t);
        prop_f_t = prop_f_r.coeffRef(row_r, col_t);
        coef_diff += g_diff_t * (prop_f_t);
        coef_conv += (1 - g_conv_t) * (rho_v_sf_t);
    };
    return (coef_diff + coef_conv);
};
double turb_k::a_face(make<double>::sp_mat& g_diff_r, make<double>::sp_mat& g_conv_r,
                      make<uint>::sp_mat& cc_fc_r, make<double>::sp_mat& rho_v_sf_r,
                      make<double>::sp_mat& prop_f_r, uint& row_r, uint& col_r)
{
    uint col_t = cc_fc_r.coeffRef(row_r, col_r);
    double g_diff_t = g_diff_r.coeffRef(row_r, col_t);
    double g_conv_t = g_conv_r.coeffRef(row_r, col_t);
    double rho_v_sf_t = rho_v_sf_r.coeffRef(row_r, col_t);
    double prop_f_t = prop_f_r.coeffRef(row_r, col_t);
    double coef_diff = (-1) * g_diff_t * prop_f_t;
    double coef_conv = g_conv_t * rho_v_sf_t;
    return (coef_diff + coef_conv);
};
double turb_k::b_center(make<double>::sp_mat& g_src_diff_r, make<double>::sp_mat& g_src_conv_r,
                        make<double>::sp_mat& body_src_r, make<uint>::sp_mat& cc_fc_r,
                        make<double>::sp_mat& rho_v_sf_r, make<double>::sp_mat& prop_f_r,
                        uint& row_r)
{
    double src_diff = 0.0;
    double src_conv = 0.0;
    uint col_r;
    uint col_t;
    double rho_v_sf_t;
    double prop_f_t;
    double g_diff_t;
    double g_conv_t;
    double body_t = body_src_r.coeffRef(row_r, row_r);
    for(SparseMatrix<uint, RowMajor>::InnerIterator it(cc_fc_r, row_r); it; it++)
    {
        col_r = it.col();
        col_t = cc_fc_r.coeffRef(row_r, col_r);
        rho_v_sf_t = rho_v_sf_r.coeffRef(row_r, col_t);
        prop_f_t = prop_f_r.coeffRef(row_r, col_t);
        g_diff_t = g_src_diff_r.coeffRef(row_r, col_t);
        g_conv_t = g_src_conv_r.coeffRef(row_r, col_t);
        src_diff += g_diff_t * prop_f_t;
        src_conv += g_conv_t * rho_v_sf_t;
    };
    return (src_diff + src_conv + body_t);
};
class turb_e : public linear
{
    double a_center(make<double>::sp_mat&, make<double>::sp_mat&, make<double>::sp_mat&,
                    make<double>::sp_mat&, uint&);
    double a_face(make<double>::sp_mat&, make<double>::sp_mat&, make<uint>::sp_mat&,
                  make<double>::sp_mat&, make<double>::sp_mat&, uint&, uint&);
    double b_center(make<double>::sp_mat&, make<double>::sp_mat&, make<double>::sp_mat&,
                    make<uint>::sp_mat&, make<double>::sp_mat&, make<double>::sp_mat&,
                    uint&);
    double inlet();
    double outlet();   
    double wall(); 
};
double turb_e::a_center(make<double>::sp_mat& g_diff_r, make<double>::sp_mat& g_conv_r,
                        make<double>::sp_mat& rho_v_sf_r, make<double>::sp_mat& prop_f_r,
                        uint& row_r)
{
    double coef_diff = 0.0;
    double coef_conv = 0.0;
    double g_diff_t;
    double g_conv_t;
    double rho_v_sf_t;
    double prop_f_t;
    for(SparseMatrix<double, RowMajor>::InnerIterator it(g_diff_r, row_r); it; it++)
    {
        uint col_t = it.col();
        g_diff_t = g_diff_r.coeffRef(row_r, col_t);
        g_conv_t = g_conv_r.coeffRef(row_r, col_t);
        rho_v_sf_t = rho_v_sf_r.coeffRef(row_r, col_t);
        prop_f_t = prop_f_r.coeffRef(row_r, col_t);
        coef_diff += g_diff_t * (prop_f_t);
        coef_conv += (1 - g_conv_t) * (rho_v_sf_t);
    };
    return (coef_diff + coef_conv);
};
double turb_e::a_face(make<double>::sp_mat& g_diff_r, make<double>::sp_mat& g_conv_r,
                      make<uint>::sp_mat& cc_fc_r, make<double>::sp_mat& rho_v_sf_r,
                      make<double>::sp_mat& prop_f_r, uint& row_r, uint& col_r)
{
    uint col_t = cc_fc_r.coeffRef(row_r, col_r);
    double g_diff_t = g_diff_r.coeffRef(row_r, col_t);
    double g_conv_t = g_conv_r.coeffRef(row_r, col_t);
    double rho_v_sf_t = rho_v_sf_r.coeffRef(row_r, col_t);
    double prop_f_t = prop_f_r.coeffRef(row_r, col_t);
    double coef_diff = (-1) * g_diff_t * prop_f_t;
    double coef_conv = g_conv_t * rho_v_sf_t;
    return (coef_diff + coef_conv);
};
double turb_e::b_center(make<double>::sp_mat& g_src_diff_r, make<double>::sp_mat& g_src_conv_r,
                        make<double>::sp_mat& body_src_r, make<uint>::sp_mat& cc_fc_r,
                        make<double>::sp_mat& rho_v_sf_r, make<double>::sp_mat& prop_f_r,
                        uint& row_r)
{
    double src_diff = 0.0;
    double src_conv = 0.0;
    uint col_r;
    uint col_t;
    double rho_v_sf_t;
    double prop_f_t;
    double g_diff_t;
    double g_conv_t;
    double body_t = body_src_r.coeffRef(row_r, row_r);
    for(SparseMatrix<uint, RowMajor>::InnerIterator it(cc_fc_r, row_r); it; it++)
    {
        col_r = it.col();
        col_t = cc_fc_r.coeffRef(row_r, col_r);
        rho_v_sf_t = rho_v_sf_r.coeffRef(row_r, col_t);
        prop_f_t = prop_f_r.coeffRef(row_r, col_t);
        g_diff_t = g_src_diff_r.coeffRef(row_r, col_t);
        g_conv_t = g_src_conv_r.coeffRef(row_r, col_t);
        src_diff += g_diff_t * prop_f_t;
        src_conv += g_conv_t * rho_v_sf_t;
    };
    return (src_diff + src_conv + body_t);
};
class s2s : public linear
{
    double a_center();
    double a_face();
    double b_center();
};
};
#endif