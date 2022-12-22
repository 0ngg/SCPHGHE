#ifndef SOLVERCFD_H
#define SOLVERCFD_H

#include"temp_util.h"
#include"temp_struct_linear.h"
#include"temp_struct_export.h"

namespace cfdsolver
{
// commons
template <class Obj, typename Type>
struct interpolate
{
    Type linear_face_itr(int C_id, int F_id, Obj& const eq__, cfdscheme::scheme& const scheme_ref, std::string what, std::string domain__)
    {
        switch(what)
        {
            case "value":
            double gc__ = scheme_ref.mesh.constants["gc"][domain__].coeffRef(C_id, F_id);
            double Cval__ = eq__.value.cvalue[domain__][C_id];
            double Fval__ = eq__.value.cvalue[domain__][F_id];
            return gc__ * Cval__ + (1 - gc__) * Fval__;
            case "grad":
            double gc__ = scheme_ref.mesh.constants["gc"][domain__].coeffRef(C_id, F_id);
            coor Cgrad__ = eq__.value.cgrad[domain__][C_id];
            coor Fgrad__ = eq__.value.cgrad[domain__][F_id];
            return gc__ * Cgrad__ + (1 - gc__) * Fgrad__;
            case "gamma":
            double gc__ = scheme_ref.mesh.constants["gc"][domain__]coeffRef(C_id, F_id);
            double Cgamma__ = eq__.gamma[domain__][C_id];
            double Fgamma__ = eq__.gamma[domain__][F_id];
            return (gc__ * Fval__ + (1 - gc__) * Cval__) / (Cval__ * Fval__);
        };
    };
    Type quick_face_itr(int C_id, int F_id, Obj& const eq__, cfdscheme::scheme& const scheme_ref, std::string what, std::string domain__)
    {
        switch(what)
        {
            case "value":
            int f_id = eq__.mesh.cc_fc[domain__].coeffRef(C_id, F_id);
            coor cdCf__ = scheme_ref.mesh.geom["dCf"][domain__].axes_to_coor(C_id, f_id);
            coor Cgrad__ = eq__.value.cgrad[C_id];
            coor Fgrad__ = eq__.value.cgrad[F_id];
            coor fgrad__ = eq__.value.fgrad[f_id];
            double Cval__ = eq__.value.cvalue[C_id];
            return Cval__ + fgrad_itr__.dot(cdCf__) / 2;
            case "grad":
            int f_id = eq__.mesh.cc_fc.coeffRef(C_id, F_id);
            double gc__ = scheme_ref.mesh.constants.gc.coeffRef(C_id, F_id);
            coor ceCF__ = scheme_ref.mesh.geom.eCF.axes_to_coor(C_id, F_id);
            double vdCf__ = scheme_ref.mesh.geom.dCf.axes_to_value(C_id,  f_id);
            coor Cgrad__ = eq__.value.cgrad[C_id];
            coor Fgrad__ = eq__.value.cgrad[F_id];
            double Cval__ = eq__.value.cvalue[C_id];
            double Fval__ = eq__.value.cvalue[F_id];
            coor fgrad_itr_ = linear_face_itr(C_id, F_id, eq__, scheme_, "grad");
            return fgrad_itr__ + ((Fval__ - Cval__) / vdCF__) * ceCF__ - (fgrad__.dot(ceCF__) * ceCF__);
        };
    };
    void least_square_itr(Obj& const eq__, cfdscheme::scheme& const scheme_ref, std::string domain__)
    {
        // calculate cgrad with optimization [A][B] = [C]
        make<int>::sp_mat& lhs_cc_ = eq__.lhs_cc_[domain__];
        make<double>::sp_mat& dCF_ = scheme_ref.mesh.geom["dCF"][domain__];
        double Cval__; double Fval__;
        coor cdCF__; double wk__;
        int F_id;
        for(int i = 0; i < cc_.outerSize(); i++)
        {
            Eigen::Matrix3d::Zero ls_lhs__;
            Eigen::Vector3d::Zero ls_rhs__;
            for(make<int>::sp_mat::InnerIterator it(lhs_cc_, i); it; it++)
            {
                Cval__ = eq__.value.cvalue[domain__][it.row()];
                Fval__ = eq__.value.cvalue[domain__][it.col()];
                cdCF__ = dCF_.axes_to_coor(it.row(), it.col());
                wk__ = 1 / dCF_.axes_to_val(it.row(), it.col());
                for(int j = 0; j < 2; j++)
                {
                    for(int k = 0; k < 2; k++)
                    {
                        ls_lhs__(i, j) += wk__ * cdCF__(i) * cdCF__(j);
                    };
                    ls_rhs__(i) += wk__ * cdCF__(i) * (Fval__ - Cval__);
                };
            };
            Eigen::Vector3d x = ls_lhs__.template bdcSvd<Eigen::ComputeThinU | Eigen::ComputeThinV>().solve(ls_rhs__);
            eq__.value.cgrad[domain__][it.row()] = x;
        };
    };
};
double calc_fluid_prop(std::string what, double P, double T, double W)
{
    // rho, miu, cp, eps, alpha
    // W = absolute humidity
    if(what.compare("rho") == 0)
    {
        std::string rho_("Vha"); std::string T_("T");
        std::string P_("P"); std::string W_("W");
        return 1 / HumidAir::HAPropsSI(rho_, P_, P, T_, T, W_, W);
    }
    else if(what.compare("miu") == 0)
    {
        std::string mu_("mu"); std::string T_("T");
        std::string P_("P"); std::string W_("W");
        return HumidAir::HAPropsSI(mu_, P_, P, T_, T, W_, W);
    }
    else if(what.compare("cp") == 0)
    {
        std::string cp_("cp_ha"); std::string T_("T");
        std::string P_("P"); std::string W_("W");
        return HumidAir::HAPropsSI(cp_, P_, P, T_, T, W_, W);
    }
    else if(what.compare("k") == 0)
    {
        std::string k_("k"); std::string T_("T");
        std::string P_("P"); std::string W_("W");
        return HumidAir::HAPropsSI(k_, P_, P, T_, T, W_, W);
    }
    else if(what.compare("alpha") == 0)
    {
        // k / (rho * cp)
        std::string rho_("Vha"); std::string k_("k"); std::string cp_("cp_ha");
        std::string T_("T"); std::string P_("P"); std::string W_("W");
        double k__ = HumidAir::HAPropsSI(k_, P_, P, T_, T, W_, W);
        double rho__ = 1 / HumidAir::HAPropsSI(rho_, P_, P, T_, T, W_, W);
        double cp__ = HumidAir::HAPropsSI(cp_, P_, P, T_, T, W_, W);
        return k__ / (rho__ * cp__);
    }
    else
    {
        return 0.0;
    };
};
// scheme object update functions outside cfdscheme namespace
void update_wall(cfdscheme::scheme& scheme_ref, cfdlinear::momentum& const u_ref, cfdlinear::momentum& const v_ref, cfdlinear::momentum& const w_ref,
                 cfdlinear::turb_k& const k_ref, cfdlinear::turb_e& const e_ref, cfdlinear::energy& const energy_ref)
{
    for(std::pair<int, double> entry : scheme_ref.wall.wall_dist)
    {
        // utau
        coor vc__{u_ref.value.cvalue["fluid"][entry.first], v_ref.value.cvalue["fluid"][entry.first], w_ref.value.cvalue["fluid"][entry.first]};
        double v__ = sqrt_sum(vc__);
        double dist__ = scheme_ref.wall.wall_dist[entry.first];
        double viu__ = scheme_ref.prop.miu["cell"][entry.first] / scheme_ref.prop.rho["cell"][entry.first];
        double u_star__ = pow(0.09, 0.25) * pow(k_ref.value.cvalue["fluid"][entry.first], 0.5);
        double utau_check__ = v__ / (5.25 + std::log(dist__) / 0.41);
        scheme_ref.wall.utau[entry.first] = std::max(utau_check__, 11.06);
        // ts
        double ts_vec1__ = k_ref.value.cvalue["fluid"][entry.first] / e_ref.value.cvalue["fluid"][entry.first];
        coor u_gradC__ = u_ref.value.cgrad["fluid"][entry.first]; coor v_gradC__ = v_ref.value.cgrad["fluid"][entry.first]; coor w_gradC__ = w_ref.value.cgrad["fluid"][entry.first];
        double St__ = pow(u_gradC__(0), 2) + pow(v_gradC__(1), 2) + pow(w_gradC__(2), 2) +
                    (pow(u_gradC__(1) + v_gradC__(0), 2) / 2) + (pow(u_gradC__(2) + w_gradC__(0), 2) / 2) + (pow(v_gradC__(2) + w_gradC__(1), 2) / 2);
        St__ = pow(St__, 0.5);
        double ts_vec2__ = 0.6 / (pow(6, 0.5) * 0.09 * St__);
        scheme_ref.wall.miut[entry.first] = std::min(ts_vec1__, ts_vec2__);
        // miut
        double rho_C__ = scheme_ref.prop.rho["cell"][entry.first];
        double miut__ = rho_C__ * 0.09 * k_ref.value.cvalue["fluid"][entry.first] * scheme_ref.wall.ts[entry.first];
        scheme_ref.wall.miut[entry.first] = miut__;
    };
};
void update_fluid_prop(cfdscheme::scheme& scheme_ref, cfdlinear::energy& const energy_ref, double AH)
{
    double pval__; double Tval__;
    for(std::pair<std::string, make<double>::map_int> entry_str : scheme_ref.prop.rho)
    {
        for(std::pair<int, double> entry_int : entry_str.second)
        {
            if(entry_str.first.compare("cell") == 0)
            {
                pval__ = scheme_ref.pressure.cvalue["fluid"][entry_int.first];
                Tval__ = energy_ref.value.cvalue["fluid"][entry_int.first];
            }
            else
            {
                pval__ = scheme_ref.pressure.fvalue["fluid"][entry_int.first];
                Tval__ = energy_ref.value.fvalue["fluid"][entry_int.first];
            };
            scheme_ref.prop.rho[entry_str.first][entry_int.first] = calc_fluid_prop("rho", pval__, Tval__, AH);
            scheme_ref.prop.miu[entry_str.first][entry_int.first] = calc_fluid_prop("miu", pval__, Tval__, AH);
            scheme_ref.prop.cp[entry_str.first][entry_int.first] = calc_fluid_prop("cp", pval__, Tval__, AH);
        };
    };
};
void update_body_term(cfdscheme::scheme& scheme_ref, cfdlinear::momentum& const u_ref, cfdlinear::momentum& const v_ref, cfdlinear::momentum& const w_ref)
{
    // rho_v_sf
    make<double>::sp_mat& rho_v_sf_ = scheme_ref.rho_v_sf;
    make<double>::map_int& rho_f_ = scheme_ref.prop.rho["face"];
    make<double>::map_int& u_fvalue_ = u_ref.value.fvalue["fluid"];
    make<double>::map_int& v_fvalue_ = v_ref.value.fvalue["fluid"];
    make<double>::map_int& w_fvalue_ = w_ref.value.fvalue["fluid"];
    axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
    coor Sf__;
    for(int i = 0; i < rho_v_sf_.outerSize(); i++)
    {
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(rho_v_sf_, i); it; it++)
        {
            coor v__(u_fvalue_[it.col()], v_fvalue_[it.col()], w_fvalue_[it.col()]);
            Sf__ = Sf_.axes_to_coor(it.row(), it.col());
            rho_v_sf_.coeffRef(it.row(), it.col()) = rho_f_[it.col()] * (v__.dot(Sf__));
        };
    };
    // phi_v
    make<double>::map_int& phi_v_ = scheme_ref.phi_v;
    make<coor>::map_int& u_fgrad_ = u_ref.value.fgrad["fluid"];
    make<coor>::map_int& v_fgrad_ = v_ref.value.fgrad["fluid"];
    make<coor>::map_int& w_fgrad_ = w_ref.value.fgrad["fluid"];
    coor u_fgrad__; coor v_fgrad__; coor w_fgrad__;
    double phi_v__; 
    for(std::pair<int, double> entry : phi_v_)
    {
        u_fgrad__ = u_fgrad_[entry.first]; v_fgrad__ = v_fgrad_[entry.first]; w_fgrad__ = w_fgrad_[entry.first]; 
        phi_v__ = 2 * (pow(u_fgrad__[0], 2) + pow(v_fgrad__[1], 2) + pow(w_fgrad__[2], 2))
                  + pow(u_fgrad__[1] + v_fgrad__[0], 2) + pow(u_fgrad__[2] + w_fgrad__[0], 2)
                  + pow(v_fgrad__[2] + w_fgrad__[1], 2);
        phi_v_[entry.first] = phi_v__;
    };
};
// scphghe iteration util
template <class V>
double check_convergence(V solv__)
{
    make<double>::sp_mat lhs__ = solv__.lhs;
    Eigen::VectorXf rhs__ = solv__.rhs.column(0);
    EIgen::VectorXf::Zeros cvalue__(solv__.rhs.rows());
    for(std::pair<std::string, make<double>::map_int> entry1 : solv__.eq.value.cvalue)
    {
        for(std::pair<int, double> entry2 : entry1.second)
        {
            cvalue__(entry2.first) = entry2.second;
        };
    };
    Eigen::VectorXf x = lhs__ * cvalue__;
    x = rhs__ - x;
    double res = 0.0;
    for(auto i : x)
    {
        res += pow(i, 2);
    };
    res = pow(res / rhs__.size(), 0.5);
    return res;
};
};
#endif