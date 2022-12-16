#include"temp_struct_linear.h"
#include"temp_func_solver.h"

#ifdef LINEARCFD_H
namespace cfdlinear
{
// make_linear
void pcorrect::make_linear(cfdscheme::scheme& const scheme_ref)
{
    append_template(make<std::string>::vec{"fluid"}, this, scheme_ref);
};
void momentum::make_linear(cfdscheme::scheme& const scheme_ref, int axis_in)
{
    this->axis = axis_in;
    append_template(make<std::string>::vec{"fluid"}, this, scheme_ref);
};
void turb_k::make_linear(cfdscheme::scheme& const scheme_ref)
{
    append_template(make<std::string>::vec{"fluid"}, this, scheme_ref);
};
void turb_e::make_linear(cfdscheme::scheme& const scheme_ref)
{
    append_template(make<std::string>::vec{"fluid"}, this, scheme_ref);
};
void energy::make_linear(cfdscheme::scheme& const scheme_ref)
{
    make<std::string>::vec which;
    for(std::pair<std::string, make<double>::sp_mat> entry : scheme_ref.mesh.cc)
    {
        if(entry.first.compare("s2s") != 0)
        {
            which.push_back(entry.first);
        };
    };
    append_template(which, this, scheme_ref);
};
void s2s::make_linear(cfdscheme::scheme& const scheme_ref)
{
    append_template(make<std::string>::vec{"s2s"}, this, scheme_ref);
};
// update_linear
void pcorrect::update_linear(cfdscheme::scheme& const scheme_ref, momentum& const u_ref, momentum& const v_ref, momentum& const w_ref)
{
    this->calc_lhs(scheme_ref, u_ref, v_ref, w_ref);
    this->calc_rhs(scheme_ref, u_ref, v_ref, w_ref);
    //solve
    this->update_correction(scheme_ref, u_ref, v_ref, w_ref);
};
void momentum::update_linear(cfdscheme::scheme& const scheme_ref, turb_k& const k_ref, momentum& const v1_ref, momentum& const v2_ref)
{
    this->calc_gamma(scheme_ref);
    this->calc_lhs(scheme_ref, v1_ref, v2_ref);
    this->calc_rhs(scheme_ref, k_ref, v1_ref, v2_ref);
    //solve
    this->calc_wall(scheme_ref, v1_ref, v2_ref);
};
void turb_k::update_linear(cfdscheme::scheme& const scheme_ref, turb_e& const e_ref, momentum& const u_ref,
                           momentum& const v_ref, momentum& const w_ref)
{
    this->calc_gamma(scheme_ref);
    this->calc_lhs(scheme_ref, u_ref, v_ref, w_ref);
    this->calc_rhs(scheme_ref, e_ref, u_ref, v_ref, w_ref);
    //solve
    this->calc_wall(scheme_ref, e_ref);
};
void turb_e::update_linear(cfdscheme::scheme& const scheme_ref, turb_k& const k_ref, momentum& const u_ref, momentum& const v_ref, momentum& const w_ref)
{
    this->calc_gamma(scheme_ref);
    this->calc_lhs(scheme_ref, u_ref, v_ref, w_ref);
    this->calc_rhs(scheme_ref, k_ref, u_ref, v_ref, w_ref);
    //solve
    this->calc_wall(scheme_ref);
};
void energy::update_linear(cfdscheme::scheme& const scheme_ref, momentum& const u_ref, momentum& const v_ref, momentum& const w_ref)
{
    this->calc_gamma(scheme_ref);
    this->calc_lhs(scheme_ref, u_ref, v_ref, w_ref);
    this->calc_rhs(scheme_ref, u_ref, v_ref, w_ref);
    //solve
    this->calc_wall(scheme_ref);
};
void s2s::update_linear(cfdscheme::scheme& const scheme_ref, energy& const energy_ref)
{
    this->calc_lhs(scheme_ref);
    this->calc_rhs(scheme_ref, energy_ref);
    //solve
    //update source
    for(std::pair<int, double> entry : this->value.cvalue["s2s"])
    {
        scheme_ref.source.value["s2s"][entry.first] = entry.second;
    };
};
// update_correcttion
void pcorrect::update_correction(cfdscheme::scheme& const scheme_ref, momentum& const u_ref, momentum& const v_ref, momentum& w_ref)
{
    make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref.prop.rho["face"];
    make<double>::map_int& rho_c_ = scheme_ref.prop.rho["cell"];
    make<double>::map_int& vol_ = scheme_ref.mesh.size["volume"];
    axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
    axes& dCF_ = scheme_ref.mesh.geom["dCF"]["fluid"];
    make<double>::sp_mat& u_aC_ = u_ref.lhs_cc["fluid"];
    make<double>::sp_mat& v_aC_ = v_ref.lhs_cc["fluid"];
    make<double>::sp_mat& w_aC_ = w_ref.lhs_cc["fluid"];
    make<coor>::map_int& pcorrect_cgrad_ = this->value.cgrad["fluid"];
    make<coor>::map_int& pcorrect_fgrad_ = this->value.fgrad["fluid"];
    double g_pcorrect__; coor pcorrect_cgrad__; coor pcorrect_fgrad__;
    coor Sf__; coor Ef__; coor Tf__; coor dCF__;
    double u_Df__; double v_Df__; double w_Df__;
    double u_DC__; double v_DC__; double w_DC__;
    int row;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
        {
            u_Df__ = (vol_[it.row()]/u_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/u_aC_.coeffRef(it.col(), it.col())) / 2;
            v_Df__ = (vol_[it.row()]/v_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/v_aC_.coeffRef(it.col(), it.col())) / 2;
            w_Df__ = (vol_[it.row()]/w_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/w_aC_.coeffRef(it.col(), it.col())) / 2;
            Sf__ = Sf_.axes_to_coor(it.row(), it.value());
            dCF__ = dCF_.axes_to_coor(it.row(), it.col());
            pcorrect_fgrad__ = pcorrect_fgrad_[it.value()];
            g_pcorrect__ = (-1) * rho_f_[it.value()] * pcorrect_fgrad__.dot(Sf__);
            u_ref.value.fvalue["fluid"][it.value()] += u_Df__ * g_pcorrect__;
            v_ref.value.fvalue["fluid"][it.value()] += v_Df__ * g_pcorrect__;
            w_ref.value.fvalue["fluid"][it.value()] += w_Df__ * g_pcorrect__;
            row = it.row();
        };
        pcorrect_cgrad__ = pcorrect_cgrad_[row];
        u_DC__ = vol_[row]/u_aC_.coeffRef(row, row);
        v_DC__ = vol_[row]/v_aC_.coeffRef(row, row);
        w_DC__ = vol_[row]/v_aC_.coeffRef(row, row);
        u_ref.value.cvalue["fluid"][row] += (-1) * rho_c_[row] * u_DC__ * pcorrect_cgrad__(0);
        v_ref.value.cvalue["fluid"][row] += (-1) * rho_c_[row] * v_DC__ * pcorrect_cgrad__(1);
        w_ref.value.cvalue["fluid"][row] += (-1) * rho_c_[row] * w_DC__ * pcorrect_cgrad__(2);
        scheme_ref.pressure.cvalue["fluid"][row] += this->value.cvalue["fluid"][row];
    };
};
//coef calc
//pcorrect
void pcorrect::calc_lhs(cfdscheme::scheme& const scheme_ref, momentum& const u_ref, momentum& const v_ref,
                        momentum& const w_ref)
{
    make<double>::sp_mat& fc_ = this->lhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->lhs_cc["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref.prop.rho["face"];
    make<double>::map_int& vol_ = scheme_ref.mesh.size["volume"];
    axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
    axes& dCF_ = scheme_ref.mesh.geom["dCF"]["fluid"];
    make<double>::sp_mat& u_aC_ = u_ref.lhs_cc["fluid"];
    make<double>::sp_mat& v_aC_ = v_ref.lhs_cc["fluid"];
    make<double>::sp_mat& w_aC_ = w_ref.lhs_cc["fluid"];
    coor Sf__; coor Ef__; coor Tf__; coor dCF__;
    double u_Df__; double v_Df__; double w_Df__;
    double Dauf__; double DauC__ = 0.0;
    int row; double aC;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        aC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
        {
            u_Df__ = (vol_[it.row()]/u_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/u_aC_.coeffRef(it.col(), it.col())) / 2;
            v_Df__ = (vol_[it.row()]/v_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/v_aC_.coeffRef(it.col(), it.col())) / 2;
            w_Df__ = (vol_[it.row()]/w_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/w_aC_.coeffRef(it.col(), it.col())) / 2;
            Sf__ = Sf_.axes_to_coor(it.row(), it.value());
            dCF__ = dCF_.axes_to_coor(it.row(), it.col());
            Dauf__ = (pow(u_Df__*Sf__(0) ,2) + pow(v_Df__*Sf__(1),2) + pow(w_Df__*Sf__(2),2)) /
                     (dCF__(0)*u_Df__*Sf__(0) + dCF__(1)*v_Df__*Sf__(1) + dCF__(2)*w_Df__*Sf__(2));
            DauC__ += (-1) * Dauf__;
            fc_.coeffRef(it.row(), it.value()) = rho_f_[it.value()]*Dauf__;
            cc_.coeffRef(it.row(), it.col()) = (-1)*rho_f_[it.value()]*Dauf__;
            row = it.row();
        };
        // boundary
        make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int::iterator it_bound = scheme_ref.mesh.bid["fluid"].find(row);
        make<std::pair<std::string, int>>::vec fc_bound;
        if(it_bound != scheme_ref.mesh.bid["fluid"].end())
        {
            for(std::pair<int, make<std::pair<std::string, int>>::vec> entry : scheme_ref.mesh.bid["fluid"][row])
            {
                fc_bound = scheme_ref.mesh.bid["fluid"][row][entry.first];
                for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                {
                    this->calc_bound_lhs(row, entry.first, j->first, j->second, scheme_ref, DauC__);
                };
            };
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; it++)
        {
            row = it.row();
            aC += fc_.coeffRef(it.row(), it.col());
        };
        cc_.coeffRef(row, row) = aC;
    };
};
void pcorrect::calc_rhs(cfdscheme::scheme& const scheme_ref, momentum& const u_ref, momentum& const v_ref,
                        momentum& const w_ref)
{
    make<double>::sp_mat& fc_ = this->rhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->rhs_cc["fluid"];
    make<double>::sp_mat& gc_ = scheme_ref.mesh.constants["gc"]["fluid"];
    cfdscheme::vinfo& pressure_ = scheme_ref.pressure;
    make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref.prop.rho["face"];
    make<double>::map_int& vol_ = scheme_ref.mesh.size["volume"];
    axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
    axes& Ef_ = scheme_ref.mesh.geom["Ef"]["fluid"];
    axes& Tf_ = scheme_ref.mesh.geom["Tf"]["fluid"];
    axes& dCF_ = scheme_ref.mesh.geom["dCF"]["fluid"];
    make<double>::sp_mat& u_aC_ = u_ref.lhs_cc["fluid"];
    make<double>::sp_mat& v_aC_ = v_ref.lhs_cc["fluid"];
    make<double>::sp_mat& w_aC_ = w_ref.lhs_cc["fluid"];
    coor mom_f__; coor mom_C__; coor mom_F__;
    double gc__;
    coor prev_pgradf__; coor prev_pgradf_itr__; coor prev_pgradf_min__;
    coor Sf__; coor Ef__; coor Tf__; coor dCF__;
    double u_Df__; double v_Df__; double w_Df__;
    double Dauf__;
    int row; double bC;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        bC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
        {
            gc__ = gc_.coeffRef(it.row(), it.value());
            prev_pgradf__ = pressure_.prev_fgrad["fluid"][it.value()];
            prev_pgradf_itr__ = gc__ * pressure_.prev_cgrad["fluid"][it.row()] + (1 - gc__) * pressure_.prev_cgrad["fluid"][it.col()];
            prev_pgradf_min__ = prev_pgradf__ - prev_pgradf_itr__;
            u_Df__ = (vol_[it.row()]/u_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/u_aC_.coeffRef(it.col(), it.col())) / 2;
            v_Df__ = (vol_[it.row()]/v_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/v_aC_.coeffRef(it.col(), it.col())) / 2;
            w_Df__ = (vol_[it.row()]/w_aC_.coeffRef(it.row(), it.row()) + vol_[it.col()]/w_aC_.coeffRef(it.col(), it.col())) / 2;
            mom_C__(0) = u_ref.value.cvalue["fluid"][it.row()]; mom_C__(1) = v_ref.value.cvalue["fluid"][it.row()]; mom_C__(2) = w_ref.value.cvalue["fluid"][it.row()];
            mom_F__(0) = u_ref.value.cvalue["fluid"][it.col()]; mom_F__(1) = v_ref.value.cvalue["fluid"][it.col()]; mom_F__(2) = w_ref.value.cvalue["fluid"][it.col()];
            mom_f__ = gc__* (mom_C__) + (1 - gc__) * (mom_F__);
            Sf__ = Sf_.axes_to_coor(it.row(), it.value());
            Ef__ = Ef_.axes_to_coor(it.row(), it.value());
            Tf__ = Tf_.axes_to_coor(it.row(), it.value());
            dCF__ = dCF_.axes_to_coor(it.row(), it.col());
            Dauf__ = (pow(u_Df__*Sf__(0) ,2) + pow(v_Df__*Sf__(1),2) + pow(w_Df__*Sf__(2),2)) /
                     (dCF__(0)*u_Df__*Sf__(0) + dCF__(1)*v_Df__*Sf__(1) + dCF__(2)*w_Df__*Sf__(2));
            fc_.coeffRef(it.row(), it.value()) = (rho_f_[it.value()] * (mom_f__.dot(Ef__) + mom_f__.dot(Tf__)) -
                                                 Dauf__ * (prev_pgradf_min__.dot(Ef__) + prev_pgradf_min__.dot(Tf__)));
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; it++)
        {
            row = it.row();
            bC += fc_.coeffRef(it.row(), it.col());
        };
        cc_.coeffRef(row, 0) = bC;
    };
};
void pcorrect::calc_bound_lhs(int row, int col, std::string check, int unique, cfdscheme::scheme& const scheme_ref, double DauC__)
{
    if(check.compare("ns") == 0)
    {
        cfdscheme::vinfo& pressure_ = scheme_ref.pressure;
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& Tf_ = scheme_ref.mesh.geom["Tf"]["fluid"];
        double pc__ = pressure_.cvalue["fluid"][row];
        coor prev_pgradc__ = pressure_.prev_cgrad["fluid"][row];
        coor prev_pgradf__ = pressure_.prev_fgrad["fluid"][col];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor Tf__ = Tf_.axes_to_coor(row, col);
        pressure_.fvalue["fluid"][col] = pc__ + (prev_pgradc__.dot(Sf__) - prev_pgradf__.dot(Tf__)) / DauC__;
    }
    else if(check.compare("in") == 0)
    {
        this->lhs_fc["fluid"].coeffRef(row, col) = scheme_ref.prop.rho["face"][col] * DauC__; 
    }
    else if(check.compare("out") == 0)
    {
        this->lhs_fc["fluid"].coeffRef(row, col) = scheme_ref.prop.rho["face"][col] * DauC__; 
    }
    else
    {
        return;
    };
};
//momentum
void momentum::calc_gamma(cfdscheme::scheme& const scheme_ref)
{
    make<double>::map_int& gamma_c_ = this->gamma["fluid"]["cell"];
    make<double>::map_int& gamma_f_ = this->gamma["fluid"]["face"];
    make<double>::sp_mat& gc_ = scheme_ref.mesh.constants["gc"]["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc["fluid"];
    make<double>::map_int& miu_c_ = scheme_ref.prop.miu["cell"];
    make<double>::map_int& miut_ = scheme_ref.wall.miut;
    for(int i = 0; i < sizeof(gamma_c_); i++)
    {
        gamma_c_[i] = miu_c_[i] + miut_[i];
    };
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
        {
            gamma_f_[it.value()] = gamma_c_[it.row()] * gamma_c_[it.col()] / (gamma_c_[it.col()] *
                                   gc_.coeffRef(it.row(), it.col()) + gamma_c_[it.row()] * (1 -
                                   gc_.coeffRef(it.row(), it.col())));
        };
    };
};
void momentum::calc_wall(cfdscheme::scheme& const scheme_ref, momentum& const v1_ref, momentum& const v2_ref)
{
    make<int>::vec& wall_cell_ = scheme_ref.mesh.cid["misc"]["wall"];
    make<double>::map_int& valuec_ = this->value.cvalue["fluid"];
    make<coor>::map_int& gradc_ = this->value.cgrad["fluid"];
    coor gradc__; double valuec__;
    coor v__; coor v_norm__; coor v_cross__;
    coor wall_perp__; coor wall_parallel__;
    double theta__;
    double d_plus__; double d_perp__; double u_tau__; double viu__;
    for(auto i = wall_cell_.begin(); i != wall_cell_.end(); i++)
    {
        v__(this->axis) = valuec_[*i];
        v__(v1_ref.axis) = v1_ref.value.cvalue["fluid"][*i];
        v__(v2_ref.axis) = v2_ref.value.cvalue["fluid"][*i];
        v_norm__ = v__ / v__.norm();
        wall_perp__ = scheme_ref.wall.wall_parallel[*i];
        v_cross__ = v_norm__.cross(wall_perp__);
        wall_parallel__ = wall_perp__.cross(v_cross__);
        wall_parallel__ = wall_parallel__ / wall_parallel__.norm();
        theta__ = std::acos(v__.dot(wall_parallel__) / v__.norm());
        if(theta__ > 90)
        {
            wall_parallel__ = (-1) * wall_parallel__;
        };
        gradc__ = gradc_[*i];
        valuec__ = valuec_[*i];
        u_tau__ = scheme_ref.wall.utau[*i];
        viu__ = scheme_ref.prop.miu["cell"][*i] / scheme_ref.prop.rho["cell"][*i];
        d_perp__ = pow(gradc__.squaredNorm() + 2*valuec__, 0.5) - gradc__.norm();
        d_plus__ = d_perp__ * u_tau__ / viu__;
        valuec_[*i] = ((std::log(d_plus__) / 0.41) + 5.25) * wall_parallel__(this->axis);
    };
};
void momentum::calc_lhs(cfdscheme::scheme& const scheme_ref, momentum& const v1_ref, momentum& const v2_ref)
{
    make<double>::sp_mat& fc_ = this->lhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->lhs_cc["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc["fluid"];
    make<double>::sp_mat& fcg_conv_aC_ = scheme_ref.mesh.constants["g_conv_aC"]["fluid"]; 
    make<double>::sp_mat& fcg_conv_aF_ = scheme_ref.mesh.constants["g_conv_aF"]["fluid"];
    make<double>::sp_mat& fcg_diff_aC_ = scheme_ref.mesh.constants["g_diff_aC"]["fluid"];
    make<double>::sp_mat& fcg_diff_aF_ = scheme_ref.mesh.constants["g_diff_aF"]["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref.prop.rho["face"];
    make<double>::sp_mat& rho_v_sf_ = scheme_ref.rho_v_sf;
    make<double>::map_int& gamma_f_ = this->gamma["fluid"]["face"];
    axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
    axes& Ef_ = scheme_ref.mesh.geom["Ef"]["fluid"];
    axes& Tf_ = scheme_ref.mesh.geom["Tf"]["fluid"];
    axes& eCF_ = scheme_ref.mesh.geom["eCF"]["fluid"];
    axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
    axes& dCF_ = scheme_ref.mesh.geom["dCF"]["fluid"];
    coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__;
    double dCFval__;
    int row; double aC;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        aC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
        {
            fc_.coeffRef(it.row(), it.value()) = fcg_conv_aC_.coeffRef(it.row(), it.value()) *
                                                 rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                                 fcg_diff_aC_.coeffRef(it.row(), it.value()) *
                                                 gamma_f_[it.value()];
            cc_.coeffRef(it.row(), it.col()) = fcg_conv_aF_.coeffRef(it.row(), it.value()) *
                                               rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                               fcg_diff_aF_.coeffRef(it.row(), it.value()) *
                                               gamma_f_[it.value()];
            row = it.row();
        };
        // boundary
        make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int::iterator it_bound = scheme_ref.mesh.bid["fluid"].find(row);
        make<std::pair<std::string, int>>::vec fc_bound;
        if(it_bound != scheme_ref.mesh.bid["fluid"].end())
        {
            for(std::pair<int, make<std::pair<std::string, int>>::vec> entry : scheme_ref.mesh.bid["fluid"][row])
            {
                fc_bound = scheme_ref.mesh.bid["fluid"][row][entry.first];
                for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                {
                    this->calc_bound_lhs(row, entry.first, j->first, j->second, scheme_ref, v1_ref, v2_ref);
                };
            };
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; it++)
        {
            row = it.row();
            aC += fc_.coeffRef(it.row(), it.col());
        };
        cc_.coeffRef(row, row) = aC;
    };
};
void momentum::calc_rhs(cfdscheme::scheme& const scheme_ref, turb_k& const k_ref, momentum& const v1_ref,
                        momentum& const v2_ref)
{
    make<double>::sp_mat& fc_ = this->rhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->rhs_cc["fluid"];
    cfdscheme::vinfo& pressure_ = scheme_ref.pressure;
    make<double>::sp_mat& gc_ = scheme_ref.mesh.constants["gc"]["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc["fluid"];
    make<double>::sp_mat& fcg_conv_aC_ = scheme_ref.mesh.constants["g_conv_aC"]["fluid"]; 
    make<double>::sp_mat& fcg_conv_aF_ = scheme_ref.mesh.constants["g_conv_aF"]["fluid"];
    make<double>::sp_mat& fcg_diff_aC_ = scheme_ref.mesh.constants["g_diff_aC"]["fluid"];
    make<double>::sp_mat& fcg_diff_aF_ = scheme_ref.mesh.constants["g_diff_aF"]["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref.prop.rho["face"];
    make<double>::map_int& rho_c_ = scheme_ref.prop.rho["cell"];
    make<double>::sp_mat& rho_v_sf_ = scheme_ref.rho_v_sf;
    make<double>::map_int& gamma_f_ = this->gamma["face"]["fluid"];
    make<double>::map_int& vol_ = scheme_ref.mesh.size["volume"];
    double gc__;
    axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
    axes& Ef_ = scheme_ref.mesh.geom["Ef"]["fluid"];
    axes& Tf_ = scheme_ref.mesh.geom["Tf"]["fluid"];
    axes& eCF_ = scheme_ref.mesh.geom["eCF"]["fluid"];
    axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
    axes& dCf_ = scheme_ref.mesh.geom["dCf"]["fluid"];
    coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__; coor dCf__;
    double Pf__; double kf__;
    coor gradf_itr__; coor gradc__;
    int row; double bC; coor bf_conv; coor bf_diff; double bf_body;
    make<std::pair<std::string, int>>::vec fc_bound;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        bC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
        {
            gc__ = gc_.coeffRef(it.row(), it.value());
            Pf__ = pressure_.fvalue["fluid"][it.value()];
            kf__ = k_ref.value.fvalue["fluid"][it.value()];
            gradc__ = this->value.cgrad["fluid"][it.row()];
            gradf_itr__ = gc__ * gradc__ + (1 - gc__) * this->value.cgrad["fluid"][it.col()];
            Ef__ = Ef_.axes_to_coor(it.row(), it.value());
            Tf__ = Tf_.axes_to_coor(it.row(), it.value());
            eCF__ = eCF_.axes_to_coor(it.row(), it.col());
            dCf__ = dCf_.axes_to_coor(it.row(), it.value());
            bf_conv = (gradf_itr__.dot(eCF__) * eCF__ - (gradc__ + gradf_itr__)) / 2;
            bf_diff = gradf_itr__ - (gradf_itr__.dot(eCF__) * eCF__);
            bf_body = Pf__ + (2 * rho_f_[it.value()] * kf__ / 3) * Sf__(this->axis);
            fc_.coeffRef(it.row(), it.value()) = (bf_conv.dot(dCf__) * rho_v_sf_.coeffRef(it.row(), it.value())) +
                                                 (bf_diff.dot(Ef__) * gamma_f_[it.value()]) + 
                                                 (bf_diff.dot(Tf__) * gamma_f_[it.value()]) + bf_body;
            row = it.row();
        };
        // boundary
        make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int::iterator it_bound = scheme_ref.mesh.bid["fluid"].find(row);
        make<std::pair<std::string, int>>::vec fc_bound;
        if(it_bound != scheme_ref.mesh.bid["fluid"].end())
        {
            for(std::pair<int, make<std::pair<std::string, int>>::vec> entry : scheme_ref.mesh.bid["fluid"][row])
            {
                fc_bound = scheme_ref.mesh.bid["fluid"][row][entry.first];
                for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                {
                    this->calc_bound_lhs(row, entry.first, j->first, j->second, scheme_ref, v1_ref, v2_ref);
                };
            };
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; it++)
        {
            row = it.row();
            bC += fc_.coeffRef(it.row(), it.col());
            if(this->axis == 1)
            {
                bC += rho_c_[it.row()] * 9.81 * vol_[it.row()];
            };
        };
        cc_.coeffRef(row, row) = bC;
    };
};
void momentum::calc_bound_lhs(int row, int col, std::string check, int unique, cfdscheme::scheme& const const scheme_ref,
                              momentum& const v1_ref, momentum& const v2_ref)
{
    if(check.compare("ns") == 0)
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        double Sfval__ = scheme_ref.mesh.size["area"][col];
        double miu_f__ = scheme_ref.prop.miu["face"][col];
        double d_perp__ = dCf__.dot(Sf__) / Sfval__;
        this->lhs_fc["fluid"].coeffRef(row, col) = (miu_f__ * Sfval__ * (1 - pow(eCf__(this->axis), 2)))/ d_perp__;
    }
    else if(check.compare("out") == 0)
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = this->value.cgrad["fluid"][row];
        v1_gradc__ = v1_ref.value.cgrad["fluid"][row];
        v2_gradc__ = v2_ref.value.cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = this->value.cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v1_ref.value.cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = v2_ref.value.cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(this->axis) = v0_f__; vf__(v1_ref.axis) = v1_f__; vf__(v2_ref.axis) = v2_f__; 
        this->lhs_fc["fluid"].coeffRef(row, col) = scheme_ref.prop.rho["face"][col] * (vf__.dot(Sf__)); 
    }
    else
    {
        return;
    };
};
void momentum::calc_bound_rhs(int row, int col, std::string check, int unique, cfdscheme::scheme& const scheme_ref,
                              momentum& const v1_ref, momentum& const v2_ref)
{
    if(check.compare("ns") == 0)
    {
        cfdscheme::vinfo& pressure_ = scheme_ref.pressure;
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        double Sfval__ = scheme_ref.mesh.size["area"][col];
        double miu_f__ = scheme_ref.prop.miu["face"][col];
        double d_perp__ = dCf__.dot(Sf__) / Sfval__;
        this->rhs_fc["fluid"].coeffRef(row, col) = (miu_f__ * Sfval__ / d_perp__) * (
                                          this->value.fvalue["fluid"][col] * (1 - pow(eCf__(this->axis), 2)) +
                                          (v1_ref.value.cvalue["fluid"][row] - v1_ref.value.fvalue["fluid"][col]) * eCf__(v1_ref.axis) * eCf__(this->axis) +
                                          (v2_ref.value.cvalue["fluid"][row] - v2_ref.value.fvalue["fluid"][col]) * eCf__(v2_ref.axis) * eCf__(this->axis)
                                          ) - (pressure_.fvalue["fluid"][col] * Sf__(this->axis));
    }
    else if(check.compare("in") == 0)
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = (-1) * eCf_.axes_to_coor(row, col);
        coor rho_v_sf__ = scheme_ref.source.value[check][unique] * eCf__;
        this->rhs_fc["fluid"].coeffRef(row, col) = scheme_ref.prop.rho["face"][col] * rho_v_sf__.dot(Sf__); 
    }
    else if(check.compare("out") == 0)
    {
        cfdscheme::vinfo& pressure_ = scheme_ref.pressure;
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = this->value.cgrad["fluid"][row];
        v1_gradc__ = v1_ref.value.cgrad["fluid"][row];
        v2_gradc__ = v2_ref.value.cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = this->value.cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v1_ref.value.cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = v2_ref.value.cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(this->axis) = v0_f__; vf__(v1_ref.axis) = v1_f__; vf__(v2_ref.axis) = v2_f__; 
        this->rhs_fc["fluid"].coeffRef(row, col) = (-1) * scheme_ref.prop.rho["face"][col] * (vf__.dot(Sf__)) * (v0_gradc__.dot(dCf__)) -
                                          (pressure_.fvalue["fluid"][col] * Sf__(this->axis));
    }
    else
    {
        return;
    };
};
//turb_k
void turb_k::calc_gamma(cfdscheme::scheme& const scheme_ref)
{
    make<double>::map_int& gamma_c_ = this->gamma["cell"]["fluid"];
    make<double>::map_int& gamma_f_ = this->gamma["face"]["fluid"];
    make<double>::sp_mat& gc_ = scheme_ref.mesh.constants["gc"]["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc["fluid"];
    make<double>::map_int& miu_c_ = scheme_ref.prop.miu["cell"];
    make<double>::map_int& miut_ = scheme_ref.wall.miut;
    for(int i = 0; i < sizeof(gamma_c_); i++)
    {
        gamma_c_[i] = miu_c_[i] + miut_[i];
    };
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
        {
            gamma_f_[it.value()] = gamma_c_[it.row()] * gamma_c_[it.col()] / (gamma_c_[it.col()] *
                                   gc_.coeffRef(it.row(), it.col()) + gamma_c_[it.row()] * (1 -
                                   gc_.coeffRef(it.row(), it.col())));
        };
    };
};
void turb_k::calc_wall(cfdscheme::scheme& const scheme_ref, turb_e& const e_ref)
{
    make<int>::vec& wall_cell_ = scheme_ref.mesh.cid["misc"]["wall"];
    make<double>::map_int& valuec_ = this->value.cvalue["fluid"];
    double Re_t__; double cmiu__;
    for(auto i = wall_cell_.begin(); i != wall_cell_.end(); i++)
    {
        Re_t__ = scheme_ref.prop.rho["cell"][*i] * pow(this->value.cvalue["fluid"][*i], 2) / (scheme_ref.prop.miu["cell"][*i] * e_ref.value.cvalue["fluid"][*i]);
        cmiu__ = 0.09 * std::exp(-3.4 / pow(1 + (Re_t__/50), 2));
        valuec_[*i] = 1 / pow(cmiu__, 0.5);
    };
};
void turb_k::calc_lhs(cfdscheme::scheme& const scheme_ref, momentum& const u_ref, momentum& const v_ref, momentum& const w_ref)
{
    make<double>::sp_mat& fc_ = this->lhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->lhs_cc["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc["fluid"];
    make<double>::sp_mat& fcg_conv_aC_ = scheme_ref.mesh.constants["g_conv_aC"]["fluid"]; 
    make<double>::sp_mat& fcg_conv_aF_ = scheme_ref.mesh.constants["g_conv_aF"]["fluid"];
    make<double>::sp_mat& fcg_diff_aC_ = scheme_ref.mesh.constants["g_diff_aC"]["fluid"];
    make<double>::sp_mat& fcg_diff_aF_ = scheme_ref.mesh.constants["g_diff_aF"]["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref.prop.rho["face"];
    make<double>::sp_mat& rho_v_sf_ = scheme_ref.rho_v_sf;
    make<double>::map_int& gamma_f_ = this->gamma["face"]["fluid"];
    axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
    axes& Ef_ = scheme_ref.mesh.geom["Ef"]["fluid"];
    axes& Tf_ = scheme_ref.mesh.geom["Tf"]["fluid"];
    axes& eCF_ = scheme_ref.mesh.geom["eCF"]["fluid"];
    axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
    axes& dCF_ = scheme_ref.mesh.geom["dCF"]["fluid"];
    coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__;
    double dCFval__;
    int row; double aC;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        aC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
        {
            fc_.coeffRef(it.row(), it.value()) = fcg_conv_aC_.coeffRef(it.row(), it.value()) *
                                                 rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                                 fcg_diff_aC_.coeffRef(it.row(), it.value()) *
                                                 gamma_f_[it.value()];
            cc_.coeffRef(it.row(), it.col()) = fcg_conv_aF_.coeffRef(it.row(), it.value()) *
                                               rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                               fcg_diff_aF_.coeffRef(it.row(), it.value()) *
                                               gamma_f_[it.value()];
            row = it.row();
        };
        // boundary
        make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int::iterator it_bound = scheme_ref.mesh.bid["fluid"].find(row);
        make<std::pair<std::string, int>>::vec fc_bound;
        if(it_bound != scheme_ref.mesh.bid["fluid"].end())
        {
            for(std::pair<int, make<std::pair<std::string, int>>::vec> entry : scheme_ref.mesh.bid["fluid"][row])
            {
                fc_bound = scheme_ref.mesh.bid["fluid"][row][entry.first];
                for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                {
                    this->calc_bound_lhs(row, entry.first, j->first, j->second, scheme_ref, u_ref, v_ref, w_ref);
                };
            };
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; it++)
        {
            row = it.row();
            aC += fc_.coeffRef(it.row(), it.col());
        };
        cc_.coeffRef(row, row) = aC;
    };
};
void turb_k::calc_rhs(cfdscheme::scheme& const scheme_ref, turb_e& const e_ref, momentum& const u_ref,
                      momentum& const v_ref, momentum& const w_ref)
{
    make<double>::sp_mat& fc_ = this->lhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->lhs_cc["fluid"];
    cfdscheme::vinfo& pressure_ = scheme_ref.pressure;
    make<double>::sp_mat& gc_ = scheme_ref.mesh.constants["gc"]["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc["fluid"];
    make<double>::sp_mat& fcg_conv_aC_ = scheme_ref.mesh.constants["g_conv_aC"]["fluid"]; 
    make<double>::sp_mat& fcg_conv_aF_ = scheme_ref.mesh.constants["g_conv_aF"]["fluid"];
    make<double>::sp_mat& fcg_diff_aC_ = scheme_ref.mesh.constants["g_diff_aC"]["fluid"];
    make<double>::sp_mat& fcg_diff_aF_ = scheme_ref.mesh.constants["g_diff_aF"]["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref.prop.rho["face"];
    make<double>::map_int& rho_c_ = scheme_ref.prop.rho["cell"];
    make<double>::sp_mat& rho_v_sf_ = scheme_ref.rho_v_sf;
    make<double>::map_int& gamma_f_ = this->gamma["face"]["fluid"];
    make<double>::map_int& vol_ = scheme_ref.mesh.size["volume"];
    double gc__;
    axes& Sf_ = scheme_ref.mesh.geom["Sf"]["face"];
    axes& Ef_ = scheme_ref.mesh.geom["Ef"]["face"];
    axes& Tf_ = scheme_ref.mesh.geom["Tf"]["face"];
    axes& eCF_ = scheme_ref.mesh.geom["eCF"]["face"];
    axes& eCf_ = scheme_ref.mesh.geom["eCf"]["face"];
    axes& dCf_ = scheme_ref.mesh.geom["dCf"]["face"];
    coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__; coor dCf__;
    double Pf__; double kf__;
    coor gradf_itr__; coor gradc__;
    int row; double bC; coor bf_conv; coor bf_diff;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        bC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
        {
            gc__ = gc_.coeffRef(it.row(), it.value());
            Pf__ = pressure_.fvalue["fluid"][it.value()];
            gradc__ = this->value.cgrad["fluid"][it.row()];
            gradf_itr__ = gc__ * gradc__ + (1 - gc__) * this->value.cgrad["fluid"][it.col()];
            Ef__ = Ef_.axes_to_coor(it.row(), it.value());
            Tf__ = Tf_.axes_to_coor(it.row(), it.value());
            eCF__ = eCF_.axes_to_coor(it.row(), it.col());
            dCf__ = dCf_.axes_to_coor(it.row(), it.value());
            bf_conv = (gradf_itr__.dot(eCF__) * eCF__ - (gradc__ + gradf_itr__)) / 2;
            bf_diff = gradf_itr__ - (gradf_itr__.dot(eCF__) * eCF__);
            fc_.coeffRef(it.row(), it.value()) = (bf_conv.dot(dCf__) * rho_v_sf_.coeffRef(it.row(), it.value())) +
                                                 (bf_diff.dot(Ef__) * gamma["face"]["fluid"][it.value()]) + 
                                                 (bf_diff.dot(Tf__) * gamma["face"]["fluid"][it.value()]);
        row = it.row();
        };
        // boundary
        make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int::iterator it_bound = scheme_ref.mesh.bid["fluid"].find(row);
        make<std::pair<std::string, int>>::vec fc_bound;
        if(it_bound != scheme_ref.mesh.bid["fluid"].end())
        {
            for(std::pair<int, make<std::pair<std::string, int>>::vec> entry : scheme_ref.mesh.bid["fluid"][row])
            {
                fc_bound = scheme_ref.mesh.bid["fluid"][row][entry.first];
                for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                {
                    this->calc_bound_lhs(row, entry.first, j->first, j->second, scheme_ref, u_ref, v_ref, w_ref);
                };
            };
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; it++)
        {
            row = it.row();
            bC += fc_.coeffRef(it.row(), it.col());
        };
        double miut__ = scheme_ref.wall.miut[row];
        double phi_v__ = scheme_ref.phi_v[row];
        double rho_c__ = rho_c_[row];
        double e_c__ = e_ref.value.cvalue["fluid"][row];
        double vol__ = vol_[row];
        cc_.coeffRef(row, row) = bC + (((miut__ * phi_v__) - (rho_c__ * e_c__)) * vol__);
    };
};
void turb_k::calc_bound_lhs(int row, int col, std::string check, int unique, cfdscheme::scheme& const scheme_ref,
                            momentum& const u_ref, momentum& const v_ref, momentum& const w_ref)
{
    if(check.compare("in") == 0)
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        this->lhs_fc["fluid"].coeffRef(row, col) = this->gamma["face"]["fluid"][col] * Sf__.norm() / dCf__.norm(); 
    }
    else if(check.compare("out") == 0)
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        // momentum
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = u_ref.value.cgrad["fluid"][row];
        v1_gradc__ = v_ref.value.cgrad["fluid"][row];
        v2_gradc__ = w_ref.value.cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = u_ref.value.cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v_ref.value.cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = w_ref.value.cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(u_ref.axis) = v0_f__; vf__(v_ref.axis) = v1_f__; vf__(w_ref.axis) = v2_f__; 
        this->lhs_fc["fluid"].coeffRef(row, col) = scheme_ref.prop.rho["face"][col] * (vf__.dot(Sf__)); 
    }
    else
    {
        return;
    };
}; 
void turb_k::calc_bound_rhs(int row, int col, std::string check, int unique, cfdscheme::scheme& const scheme_ref,
                            momentum& const u_ref, momentum& const v_ref, momentum& const w_ref)
{
    if(check.compare("in") == 0)
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        coor vf__ = (-1) * scheme_ref.source.value[check][unique] * Sf__;
        double inlet_k = 1/2 * 0.01 * (vf__.dot(vf__));
        this->rhs_fc["fluid"].coeffRef(row, col) = (-1) * this->gamma["face"]["fluid"][col] * Sf__.norm() * inlet_k / dCf__.norm();
    }
    else if(check.compare("out") == 0)
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        // momentum
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = u_ref.value.cgrad["fluid"][row];
        v1_gradc__ = v_ref.value.cgrad["fluid"][row];
        v2_gradc__ = w_ref.value.cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = u_ref.value.cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v_ref.value.cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = w_ref.value.cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(u_ref.axis) = v0_f__; vf__(v_ref.axis) = v1_f__; vf__(w_ref.axis) = v2_f__;
        // turb_k
        coor k_gradc__; coor k_gradf__;
        k_gradc__ = this->value.cgrad["fluid"][row];
        k_gradf__ = k_gradc__ - (k_gradc__.dot(eCf__) * eCf__);
        this->rhs_fc["fluid"].coeffRef(row, col) = (-1) * scheme_ref.prop.rho["face"][col] * (vf__.dot(Sf__)) * (k_gradf__.dot(dCf__)); 
    }
    else
    {
        return;
    };
};
//turb_e
void turb_e::calc_gamma(cfdscheme::scheme& const scheme_ref)
{
    make<double>::map_int& gamma_c_ = this->gamma["cell"]["fluid"];
    make<double>::map_int& gamma_f_ = this->gamma["face"]["fluid"];
    make<double>::sp_mat& gc_ = scheme_ref.mesh.constants["gc"]["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc["fluid"];
    make<double>::map_int& miu_c_ = scheme_ref.prop.miu["cell"];
    make<double>::map_int& miut_ = scheme_ref.wall.miut;
    for(int i = 0; i < sizeof(gamma_c_); i++)
    {
        gamma_c_[i] = miu_c_[i] + miut_[i] / 1.3;
    };
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
        {
            gamma_f_[it.value()] = gamma_c_[it.row()] * gamma_c_[it.col()] / (gamma_c_[it.col()] *
                                   gc_.coeffRef(it.row(), it.col()) + gamma_c_[it.row()] * (1 -
                                   gc_.coeffRef(it.row(), it.col())));
        };
    };
};
void turb_e::calc_wall(cfdscheme::scheme& const scheme_ref)
{
    make<int>::vec& wall_cell_ = scheme_ref.mesh.cid["misc"]["wall"];
    make<double>::map_int& valuec_ = this->value.cvalue["fluid"];
    make<coor>::map_int& gradc_ = this->value.cgrad["fluid"];
    coor gradc__; double valuec__;
    double d_perp__; double u_tau__; double viu__;
    for(auto i = wall_cell_.begin(); i != wall_cell_.end(); i++)
    {
        gradc__ = gradc_[*i];
        valuec__ = valuec_[*i];
        u_tau__ = scheme_ref.wall.utau[*i];
        viu__ = scheme_ref.prop.miu["cell"][*i] / scheme_ref.prop.rho["cell"][*i];
        d_perp__ = pow(gradc__.squaredNorm() + 2*valuec__, 0.5) - gradc__.norm();
        valuec_[*i] = viu__ / (u_tau__ * 0.41 * d_perp__);
    };
};
void turb_e::calc_lhs(cfdscheme::scheme& const scheme_ref, momentum& const u_ref, momentum& const v_ref, momentum& const w_ref)
{
    make<double>::sp_mat& fc_ = this->lhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->lhs_cc["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc["fluid"];
    make<double>::sp_mat& fcg_conv_aC_ = scheme_ref.mesh.constants["g_conv_aC"]["fluid"]; 
    make<double>::sp_mat& fcg_conv_aF_ = scheme_ref.mesh.constants["g_conv_aF"]["fluid"];
    make<double>::sp_mat& fcg_diff_aC_ = scheme_ref.mesh.constants["g_diff_aC"]["fluid"];
    make<double>::sp_mat& fcg_diff_aF_ = scheme_ref.mesh.constants["g_diff_aF"]["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref.prop.rho["face"];
    make<double>::sp_mat& rho_v_sf_ = scheme_ref.rho_v_sf;
    make<double>::map_int& gamma_f_ = this->gamma["face"]["fluid"];
    axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
    axes& Ef_ = scheme_ref.mesh.geom["Ef"]["fluid"];
    axes& Tf_ = scheme_ref.mesh.geom["Tf"]["fluid"];
    axes& eCF_ = scheme_ref.mesh.geom["eCF"]["fluid"];
    axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
    axes& dCF_ = scheme_ref.mesh.geom["dCF"]["fluid"];
    coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__;
    double dCFval__;
    int row; double aC;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        aC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
        {
            fc_.coeffRef(it.row(), it.value()) = fcg_conv_aC_.coeffRef(it.row(), it.value()) *
                                                 rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                                 fcg_diff_aC_.coeffRef(it.row(), it.value()) *
                                                 gamma_f_[it.value()];
            cc_.coeffRef(it.row(), it.col()) = fcg_conv_aF_.coeffRef(it.row(), it.value()) *
                                               rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                               fcg_diff_aF_.coeffRef(it.row(), it.value()) *
                                               gamma_f_[it.value()];
            row = it.row();
        };
        // boundary
        make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int::iterator it_bound = scheme_ref.mesh.bid["fluid"].find(row);
        make<std::pair<std::string, int>>::vec fc_bound;
        if(it_bound != scheme_ref.mesh.bid["fluid"].end())
        {
            for(std::pair<int, make<std::pair<std::string, int>>::vec> entry : scheme_ref.mesh.bid["fluid"][row])
            {
                fc_bound = scheme_ref.mesh.bid["fluid"][row][entry.first];
                for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                {
                    this->calc_bound_lhs(row, entry.first, j->first, j->second, scheme_ref, u_ref, v_ref, w_ref);
                };
            };
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; it++)
        {
            row = it.row();
            aC += fc_.coeffRef(it.row(), it.col());
        };
        cc_.coeffRef(row, row) = aC;
    };
};
void turb_e::calc_rhs(cfdscheme::scheme& const scheme_ref, turb_k& const k_ref, momentum& const u_ref, momentum& const v_ref, momentum& const w_ref)
{
    make<double>::sp_mat& fc_ = this->lhs_fc["fluid"];
    make<double>::sp_mat& cc_ = this->lhs_cc["fluid"];
    cfdscheme::vinfo& pressure_ = scheme_ref.pressure;
    make<double>::sp_mat& gc_ = scheme_ref.mesh.constants["gc"]["fluid"];
    make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc["fluid"];
    make<double>::sp_mat& fcg_conv_aC_ = scheme_ref.mesh.constants["g_conv_aC"]["fluid"]; 
    make<double>::sp_mat& fcg_conv_aF_ = scheme_ref.mesh.constants["g_conv_aF"]["fluid"];
    make<double>::sp_mat& fcg_diff_aC_ = scheme_ref.mesh.constants["g_diff_aC"]["fluid"];
    make<double>::sp_mat& fcg_diff_aF_ = scheme_ref.mesh.constants["g_diff_aC"]["fluid"];
    make<double>::map_int& rho_f_ = scheme_ref.prop.rho["face"];
    make<double>::map_int& rho_c_ = scheme_ref.prop.rho["cell"];
    make<double>::sp_mat& rho_v_sf_ = scheme_ref.rho_v_sf;
    make<double>::map_int& gamma_f_ = this->gamma["face"]["fluid"];
    make<double>::map_int& vol_ = scheme_ref.mesh.size["volume"];
    double gc__;
    axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
    axes& Ef_ = scheme_ref.mesh.geom["Ef"]["fluid"];
    axes& Tf_ = scheme_ref.mesh.geom["Tf"]["fluid"];
    axes& eCF_ = scheme_ref.mesh.geom["eCF"]["fluid"];
    axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
    axes& dCf_ = scheme_ref.mesh.geom["dCf"]["fluid"];
    coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__; coor dCf__;
    double Pf__; double kf__;
    coor gradf_itr__; coor gradc__;
    int row; double bC; coor bf_conv; coor bf_diff;
    double ce_1__ = 1.44; double ce_2__; double Re_t__;
    for(int i = 0; i < cc_fc_.outerSize(); i++)
    {
        bC = 0.0;
        for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
        {
            gc__ = gc_.coeffRef(it.row(), it.value());
            Pf__ = pressure_.fvalue["fluid"][it.value()];
            gradc__ = this->value.cgrad["fluid"][it.row()];
            gradf_itr__ = gc__ * gradc__ + (1 - gc__) * this->value.cgrad["fluid"][it.col()];
            Ef__ = Ef_.axes_to_coor(it.row(), it.value());
            Tf__ = Tf_.axes_to_coor(it.row(), it.value());
            eCF__ = eCF_.axes_to_coor(it.row(), it.col());
            dCf__ = dCf_.axes_to_coor(it.row(), it.value());
            bf_conv = (gradf_itr__.dot(eCF__) * eCF__ - (gradc__ + gradf_itr__)) / 2;
            bf_diff = gradf_itr__ - (gradf_itr__.dot(eCF__) * eCF__);
            fc_.coeffRef(it.row(), it.value()) = (bf_conv.dot(dCf__) * rho_v_sf_.coeffRef(it.row(), it.value())) +
                                                 (bf_diff.dot(Ef__) * gamma_f_[it.value()]) + 
                                                 (bf_diff.dot(Tf__) * gamma_f_[it.value()]);
            row = it.row();
        };
        // boundary
        make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int::iterator it_bound = scheme_ref.mesh.bid["fluid"].find(row);
        make<std::pair<std::string, int>>::vec fc_bound;
        if(it_bound != scheme_ref.mesh.bid["fluid"].end())
        {
            for(std::pair<int, make<std::pair<std::string, int>>::vec> entry : scheme_ref.mesh.bid["fluid"][row])
            {
                fc_bound = scheme_ref.mesh.bid["fluid"][row][entry.first];
                for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                {
                    this->calc_bound_lhs(row, entry.first, j->first, j->second, scheme_ref, u_ref, v_ref, w_ref);
                };
            };
        };
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; it++)
        {
            row = it.row();
            bC += fc_.coeffRef(it.row(), it.col());
        };
        double ts__ = scheme_ref.wall.ts[row];
        double miut__ = scheme_ref.wall.miut[row];
        double phi_v__ = scheme_ref.phi_v[row];
        double rho_c__ = rho_c_[row];
        double e_c__ = this->value.cvalue["fluid"][row];
        double vol__ = vol_[row];
        Re_t__ = scheme_ref.prop.rho["cell"][row] * pow(k_ref.value.cvalue["fluid"][row], 2) / (scheme_ref.prop.miu["cell"][row] * this->value.cvalue["fluid"][row]);
        ce_2__ = 1.92*(1 - 0.3*std::exp((-1) * pow(Re_t__, 2)));
        cc_.coeffRef(row, row) = bC + ((ce_1__ * miut__ * phi_v__ / ts__) - (ce_2__ * scheme_ref.prop.rho["cell"][row] *
                                 this->value.cvalue["fluid"][row] / ts__))*vol__;
    };
};
void turb_e::calc_bound_lhs(int row, int col, std::string check, int unique, cfdscheme::scheme& const scheme_ref,
                            momentum& const u_ref, momentum& const v_ref, momentum& const w_ref)
{
    if(check.compare("in") == 0)
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& dCf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        this->lhs_fc["fluid"].coeffRef(row, col) = this->gamma["face"]["fluid"][col] * Sf__.norm() / dCf__.norm(); 
    }
    else if(check.compare("out") == 0)
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        // momentum
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = u_ref.value.cgrad["fluid"][row];
        v1_gradc__ = v_ref.value.cgrad["fluid"][row];
        v2_gradc__ = w_ref.value.cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = u_ref.value.cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v_ref.value.cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = w_ref.value.cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(u_ref.axis) = v0_f__; vf__(v_ref.axis) = v1_f__; vf__(w_ref.axis) = v2_f__; 
        this->lhs_fc["fluid"].coeffRef(row, col) = scheme_ref.prop.rho["face"][col] * (vf__.dot(Sf__)); 
    }
    else
    {
        return;
    };
}; 
void turb_e::calc_bound_rhs(int row, int col, std::string check, int unique, cfdscheme::scheme& const scheme_ref,
                            momentum& const u_ref, momentum& const v_ref, momentum& const w_ref)
{
    if(check.compare("in") == 0)
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        coor vf__ = (-1) * scheme_ref.source.value[check][unique] * Sf__;
        make<double>::map_int& rho_f_ = scheme_ref.prop.rho["face"];
        make<double>::map_int& miu_f_ = scheme_ref.prop.miu["face"];
        make<double>::map_int& miut_ = scheme_ref.wall.miut;
        double rho_f__ = rho_f_[col]; double miut__ = miut_[row];
        double inlet_k = 1/2 * 0.01 * (vf__.dot(vf__));
        double inlet_e = 0.09 * rho_f__ * pow(inlet_k, 2) / miut__;
        this->rhs_fc["fluid"].coeffRef(row, col) = (-1) * this->gamma["face"]["fluid"][col] * Sf__.norm() * inlet_e / dCf__.norm(); 
    }
    else if(check.compare("out") == 0)
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        // momentum
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = u_ref.value.cgrad["fluid"][row];
        v1_gradc__ = v_ref.value.cgrad["fluid"][row];
        v2_gradc__ = w_ref.value.cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = u_ref.value.cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v_ref.value.cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = w_ref.value.cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(u_ref.axis) = v0_f__; vf__(v_ref.axis) = v1_f__; vf__(w_ref.axis) = v2_f__;
        // turb_e
        coor e_gradc__; coor e_gradf__;
        e_gradc__ = this->value.cgrad["fluid"][row];
        e_gradf__ = e_gradc__ - (e_gradc__.dot(eCf__) * eCf__);
        this->rhs_fc["fluid"].coeffRef(row, col) = (-1) * scheme_ref.prop.rho["face"][col] * (vf__.dot(Sf__)) * (e_gradf__.dot(dCf__)); 
    }
    else
    {
        return;
    };
};
//energy
void energy::calc_gamma(cfdscheme::scheme& const scheme_ref)
{
    for(std::pair<std::string, make<double>::sp_mat> entry : this->lhs_cc)
    {
        make<double>::map_int& gamma_c_ = this->gamma["cell"][entry.first];
        make<double>::map_int& gamma_f_ = this->gamma["face"][entry.first];
        make<double>::sp_mat& gc_ = scheme_ref.mesh.constants["gc"][entry.first];
        make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc[entry.first];
        make<double>::map_int& k_ = scheme_ref.prop.k["face"];
        make<double>::map_int& cp_c_ = scheme_ref.prop.cp["cell"];
        make<double>::map_int& alpha_ = scheme_ref.prop.alpha["cell"];
        make<double>::map_int& rho_c_ = scheme_ref.prop.rho["cell"];
        make<double>::map_int& miu_c_ = scheme_ref.prop.miu["cell"];
        make<double>::map_int& miut_ = scheme_ref.wall.miut;
        double Pr__;
        if(entry.first.compare("fluid") == 0)
        {
            for(int i = 0; i < cc_fc_.outerSize(); i++)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
                {
                    Pr__ = miu_c_[it.row()] / (rho_c_[it.row()] * alpha_[it.row()]);
                    gamma_c_[it.row()] = cp_c_[it.row()] * (miu_c_[it.row()] / Pr__ + miut_[it.row()] / Pr__ );
                    break;
                };
            };
            for(int i = 0; i < cc_fc_.outerSize(); i++)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
                {
                    gamma_f_[it.value()] = gamma_c_[it.row()] * gamma_c_[it.col()] / (gamma_c_[it.col()] *
                                        gc_.coeffRef(it.row(), it.col()) + gamma_c_[it.row()] * (1 -
                                        gc_.coeffRef(it.row(), it.col())));
                };
            };
        }
        else if(entry.first.compare("solid") == 0)
        {
            for(int i = 0; i < cc_fc_.outerSize(); i++)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
                {
                    gamma_c_[it.row()] = k_[it.row()];
                };
            };
            for(int i = 0; i < cc_fc_.outerSize(); i++)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
                {
                    gamma_f_[it.value()] = gamma_c_[it.row()] * gamma_c_[it.col()] / (gamma_c_[it.col()] *
                                        gc_.coeffRef(it.row(), it.col()) + gamma_c_[it.row()] * (1 -
                                        gc_.coeffRef(it.row(), it.col())));
                };
            };
        }
        else if(entry.first.compare("conj") == 0)
        {
            for(int i = 0; i < cc_fc_.outerSize(); i++)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
                {
                    gamma_f_[it.value()] = gamma_c_[it.row()] * gamma_c_[it.col()] / (gamma_c_[it.col()] *
                                        gc_.coeffRef(it.row(), it.col()) + gamma_c_[it.row()] * (1 -
                                        gc_.coeffRef(it.row(), it.col())));
                };
            };
        };
    };
};
void energy::calc_wall(cfdscheme::scheme& const scheme_ref)
{
    make<int>::vec& wall_cell_ = scheme_ref.mesh.cid["misc"]["wall"];
    make<double>::map_int& valuec_ = this->value.cvalue["fluid"];
    make<coor>::map_int& gradc_ = this->value.cgrad["fluid"];
    coor gradc__; double valuec__;
    double d_plus__; double d_perp__; double u_tau__; double viu__; double Pr__; double beta__;
    for(auto i = wall_cell_.begin(); i != wall_cell_.end(); i++)
    {
        gradc__ = gradc_[*i];
        valuec__ = valuec_[*i];
        u_tau__ = scheme_ref.wall.utau[*i];
        viu__ = scheme_ref.prop.miu["cell"][*i] / scheme_ref.prop.rho["cell"][*i];
        d_perp__ = pow(gradc__.squaredNorm() + 2*valuec__, 0.5) - gradc__.norm();
        d_plus__ = d_perp__ * u_tau__ / viu__;
        Pr__ = viu__ / scheme_ref.prop.alpha["cell"][*i];
        beta__ = pow(3.85 * pow(Pr__, 1/3), 2)  + 2.12 * std::log(Pr__);
        valuec_[*i] = 2.12 * std::log(d_plus__) + beta__;
    };
};
void energy::calc_lhs(cfdscheme::scheme& const scheme_ref, momentum& const u_ref, momentum& const v_ref, momentum& const w_ref)
{
    for(std::pair<std::string, make<double>::sp_mat> entry : this->lhs_cc)
    {
        make<double>::sp_mat& fc_ = this->lhs_fc[entry.first];
        make<double>::sp_mat& cc_ = this->lhs_cc[entry.first];
        make<double>::sp_mat& gc_ = scheme_ref.mesh.constants["gc"][entry.first];
        make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc[entry.first];
        make<double>::sp_mat& fcg_conv_aC_ = scheme_ref.mesh.constants["g_conv_aC"][entry.first];
        make<double>::sp_mat& fcg_conv_aF_ = scheme_ref.mesh.constants["g_conv_aF"][entry.first];
        make<double>::sp_mat& fcg_diff_aC_ = scheme_ref.mesh.constants["g_diff_aC"][entry.first];
        make<double>::sp_mat& fcg_diff_aF_ = scheme_ref.mesh.constants["g_diff_aF"][entry.first];
        make<double>::map_int& rho_f_ = scheme_ref.prop.rho["face"];
        make<double>::sp_mat& rho_v_sf_ = scheme_ref.rho_v_sf;
        make<double>::map_int& gamma_f_ = this->gamma["face"][entry.first];
        axes& Sf_ = scheme_ref.mesh.geom["Sf"][entry.first];
        axes& Ef_ = scheme_ref.mesh.geom["Ef"][entry.first];
        axes& Tf_ = scheme_ref.mesh.geom["Tf"][entry.first];
        axes& eCF_ = scheme_ref.mesh.geom["eCF"][entry.first];
        axes& eCf_ = scheme_ref.mesh.geom["eCf"][entry.first];
        axes& dCF_ = scheme_ref.mesh.geom["dCF"][entry.first];
        coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__; coor dCf__;
        double Pf__; double kf__;
        coor gradf_itr__; coor gradc__;
        int row; double bC; coor bf_conv; coor bf_diff;
        double dCFval__; double cp_f__;
        int row; double aC;
        for(int i = 0; i < cc_.outerSize(); i++)
        {
            aC = 0.0;
            if(entry.first.compare("fluid") == 0)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
                {
                    cp_f__ = scheme_ref.prop.cp["face"][it.value()];
                    fc_.coeffRef(it.row(), it.value()) = cp_f__ * fcg_conv_aC_.coeffRef(it.row(), it.value()) *
                                                        rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                                        fcg_diff_aC_.coeffRef(it.row(), it.value()) *
                                                        gamma_f_[it.value()];
                    cc_.coeffRef(it.row(), it.col()) = cp_f__ * fcg_conv_aF_.coeffRef(it.row(), it.value()) *
                                                    rho_v_sf_.coeffRef(it.row(), it.value()) + 
                                                    fcg_diff_aF_.coeffRef(it.row(), it.value()) *
                                                    gamma_f_[it.value()];
                    row = it.row();
                };
                // boundary
                make<make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>::map_str::iterator it_domain = scheme_ref.mesh.bid.find("fluid");
                if(it_domain != scheme_ref.mesh.bid.end())
                {
                    make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int::iterator it_bound = scheme_ref.mesh.bid["fluid"].find(row);
                    make<std::pair<std::string, int>>::vec fc_bound;
                    if(it_bound != scheme_ref.mesh.bid[entry.first].end())
                    {
                        for(std::pair<int, make<std::pair<std::string, int>>::vec> entry_bound : scheme_ref.mesh.bid["fluid"][row])
                        {
                            fc_bound = scheme_ref.mesh.bid["fluid"][row][entry_bound.first];
                            for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                            {
                                this->calc_bound_lhs(row, entry_bound.first, j->first, j->second, scheme_ref, "fluid", u_ref, v_ref, w_ref);
                            };
                        };
                    };
                };
                for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; it++)
                {
                    row = it.row();
                    aC += fc_.coeffRef(it.row(), it.col());
                };
                cc_.coeffRef(row, row) = aC;
            }
            else if(entry.first.compare("solid") == 0)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
                {
                    fc_.coeffRef(it.row(), it.value()) = fcg_diff_aC_.coeffRef(it.row(), it.value()) *
                                                        gamma_f_[it.value()];
                    cc_.coeffRef(it.row(), it.col()) = fcg_diff_aF_.coeffRef(it.row(), it.value()) *
                                                    gamma_f_[it.value()];
                    row = it.row();
                };
                // boundary
                make<make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>::map_str::iterator it_domain = scheme_ref.mesh.bid.find("solid");
                if(it_domain != scheme_ref.mesh.bid.end())
                {
                    make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int::iterator it_bound = scheme_ref.mesh.bid["solid"].find(row);
                    make<std::pair<std::string, int>>::vec fc_bound;
                    if(it_bound != scheme_ref.mesh.bid[entry.first].end())
                    {
                        for(std::pair<int, make<std::pair<std::string, int>>::vec> entry_bound : scheme_ref.mesh.bid["solid"][row])
                        {
                            fc_bound = scheme_ref.mesh.bid["solid"][row][entry_bound.first];
                            for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                            {
                                this->calc_bound_lhs(row, entry_bound.first, j->first, j->second, scheme_ref, "solid", u_ref, v_ref, w_ref);
                            };
                        };
                    };
                };
                for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; it++)
                {
                    row = it.row();
                    aC += fc_.coeffRef(it.row(), it.col());
                };
                cc_.coeffRef(row, row) = aC;
            }
            else if(entry.first.compare("conj") == 0)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
                {
                    fc_.coeffRef(it.row(), it.value()) = 0;
                    cc_.coeffRef(it.row(), it.col()) = 0;
                    row = it.row();
                };
                for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; it++)
                {
                    row = it.row();
                    aC += fc_.coeffRef(it.row(), it.col());
                };
                cc_.coeffRef(row, row) = aC;
            };
        };
    };
};
void energy::calc_rhs(cfdscheme::scheme& const scheme_ref, momentum& const u_ref, momentum& const v_ref, momentum& const w_ref)
{
    for(std::pair<std::string, make<double>::sp_mat> entry : this->rhs_cc)
    {
        make<double>::sp_mat& fc_ = this->rhs_fc[entry.first];
        make<double>::sp_mat& cc_ = this->rhs_cc[entry.first];
        make<double>::sp_mat& gc_ = scheme_ref.mesh.constants["gc"][entry.first];
        make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc[entry.first];
        make<double>::sp_mat& fcg_conv_aC_ = scheme_ref.mesh.constants["g_conv_aC"][entry.first];
        make<double>::sp_mat& fcg_conv_aF_ = scheme_ref.mesh.constants["g_conv_aF"][entry.first];
        make<double>::sp_mat& fcg_diff_aC_ = scheme_ref.mesh.constants["g_diff_aC"][entry.first];
        make<double>::sp_mat& fcg_diff_aF_ = scheme_ref.mesh.constants["g_diff_aF"][entry.first];
        make<double>::map_int& rho_f_ = scheme_ref.prop.rho["face"];
        make<double>::sp_mat& rho_v_sf_ = scheme_ref.rho_v_sf;
        make<double>::map_int& gamma_f_ = this->gamma["face"][entry.first];
        make<double>::map_int& vol_ = scheme_ref.mesh.size["volume"];
        double gc__;
        axes& Sf_ = scheme_ref.mesh.geom["Sf"][entry.first];
        axes& Ef_ = scheme_ref.mesh.geom["Ef"][entry.first];
        axes& Tf_ = scheme_ref.mesh.geom["Tf"][entry.first];
        axes& eCF_ = scheme_ref.mesh.geom["eCF"][entry.first];
        axes& eCf_ = scheme_ref.mesh.geom["eCf"][entry.first];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"][entry.first];
        coor Sf__; coor Ef__; coor Tf__; coor eCF__; coor eCf__; coor dCf__;
        double cp__;
        coor gradf_itr__; coor gradc__;
        int row; double bC; coor bf_conv; coor bf_diff; double bf_body;
        for(int i = 0; i < cc_.outerSize(); i++)
        {
            bC = 0.0;
            if(entry.first.compare("fluid") == 0)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
                {
                    gc__ = gc_.coeffRef(it.row(), it.value());
                    gradc__ = this->value.cgrad[entry.first][it.row()];
                    gradf_itr__ = gc__ * gradc__ + (1 - gc__) * this->value.cgrad[entry.first][it.col()];
                    Ef__ = Ef_.axes_to_coor(it.row(), it.value());
                    Tf__ = Tf_.axes_to_coor(it.row(), it.value());
                    eCF__ = eCF_.axes_to_coor(it.row(), it.col());
                    dCf__ = dCf_.axes_to_coor(it.row(), it.value());
                    bf_conv = (gradf_itr__.dot(eCF__) * eCF__ - (gradc__ + gradf_itr__)) / 2;
                    bf_diff = gradf_itr__ - (gradf_itr__.dot(eCF__) * eCF__);
                    fc_.coeffRef(it.row(), it.value()) = (bf_conv.dot(dCf__) * rho_v_sf_.coeffRef(it.row(), it.value())) +
                                                        (bf_diff.dot(Ef__) * gamma_f_[it.value()]) + 
                                                        (bf_diff.dot(Tf__) * gamma_f_[it.value()]);
                };
                // boundary
                make<make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>::map_str::iterator it_domain = scheme_ref.mesh.bid.find("fluid");
                if(it_domain != scheme_ref.mesh.bid.end())
                {
                    make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int::iterator it_bound = scheme_ref.mesh.bid["fluid"].find(row);
                    make<std::pair<std::string, int>>::vec fc_bound;
                    if(it_bound != scheme_ref.mesh.bid[entry.first].end())
                    {
                        for(std::pair<int, make<std::pair<std::string, int>>::vec> entry_bound : scheme_ref.mesh.bid["fluid"][row])
                        {
                            fc_bound = scheme_ref.mesh.bid["fluid"][row][entry_bound.first];
                            for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                            {
                                this->calc_bound_lhs(row, entry_bound.first, j->first, j->second, scheme_ref, "fluid", u_ref, v_ref, w_ref);
                            };
                        };
                    };
                };
                for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; it++)
                {
                    row = it.row();
                    bC += fc_.coeffRef(it.row(), it.col());
                };
                double miu_c__ = scheme_ref.prop.miu["cell"][row];
                double miut__ = scheme_ref.wall.miut[row];
                double phi_v__ = scheme_ref.phi_v[row];
                double vol__ = vol_[row];
                cc_.coeffRef(row, row) = bC + ((miu_c__ + miut__) * phi_v__ * vol__)/sizeof(this->rhs_cc);
            }
            else if(entry.first.compare("solid") == 0)
            {
                for(Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(cc_fc_, i); it; it++)
                {
                    gc__ = gc_.coeffRef(it.row(), it.value());
                    gradc__ = this->value.cgrad[entry.first][it.row()];
                    gradf_itr__ = gc__ * gradc__ + (1 - gc__) * this->value.cgrad[entry.first][it.col()];
                    Ef__ = Ef_.axes_to_coor(it.row(), it.value());
                    Tf__ = Tf_.axes_to_coor(it.row(), it.value());
                    eCF__ = eCF_.axes_to_coor(it.row(), it.col());
                    bf_diff = gradf_itr__ - (gradf_itr__.dot(eCF__) * eCF__);
                    fc_.coeffRef(it.row(), it.value()) = (bf_diff.dot(Ef__) * gamma_f_[it.value()]) + 
                                                        (bf_diff.dot(Tf__) * gamma_f_[it.value()]);
                };
                // boundary
                make<make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>::map_str::iterator it_domain = scheme_ref.mesh.bid.find("solid");
                if(it_domain != scheme_ref.mesh.bid.end())
                {
                    make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int::iterator it_bound = scheme_ref.mesh.bid["solid"].find(row);
                    make<std::pair<std::string, int>>::vec fc_bound;
                    if(it_bound != scheme_ref.mesh.bid[entry.first].end())
                    {
                        for(std::pair<int, make<std::pair<std::string, int>>::vec> entry_bound : scheme_ref.mesh.bid["solid"][row])
                        {
                            fc_bound = scheme_ref.mesh.bid["solid"][row][entry_bound.first];
                            for(auto j = fc_bound.begin(); j != fc_bound.end(); j++)
                            {
                                this->calc_bound_lhs(row, entry_bound.first, j->first, j->second, scheme_ref, "solid", u_ref, v_ref, w_ref);
                            };
                        };
                    };
                };
                for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(fc_, i); it; it++)
                {
                    row = it.row();
                    bC += fc_.coeffRef(it.row(), it.col());
                };
                double miu_c__ = scheme_ref.prop.miu["cell"][row];
                double miut__ = scheme_ref.wall.miut[row];
                double phi_v__ = scheme_ref.phi_v[row];
                double vol__ = vol_[row];
                cc_.coeffRef(row, row) = bC + ((miu_c__ + miut__) * phi_v__ * vol__)/sizeof(this->rhs_cc); 
            };
        };
    };
};
void energy::calc_bound_lhs(int row, int col, std::string check, int unique, cfdscheme::scheme& const scheme_ref, std::string domain__, momentum& const u_ref, momentum& const v_ref, momentum& const w_ref)
{
    if(check.compare("temp"))
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"][domain__];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"][domain__];
        make<double>::map_int& gamma_f_ = this->gamma["face"][domain__];
        this->lhs_fc[domain__].coeffRef(row, col) = gamma_f_[col] * Sf_.axes_to_val(row, col) / dCf_.axes_to_val(row, col);
    }
    else if(check.compare("hamb"))
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"][domain__];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"][domain__];
        make<double>::map_int& gamma_f_ = this->gamma["face"][domain__];
        // h_sky
        double g_hamb = gamma_f_[col] * Sf_.axes_to_val(row, col) / dCf_.axes_to_val(row, col);
        double T_f__ = this->value.fvalue[domain__][col];
        double T_amb__ = scheme_ref.source.value["hamb"][0]; double T_sky__ = 0.0552 * pow(T_amb__, 1.5);
        double eps_C__ = scheme_ref.prop.eps["cell"][row];
        double h_sky__ = 5.67 * pow(10, -8) * eps_C__ * (T_f__ + T_sky__) * (pow(T_f__, 2) + pow(T_sky__, 2)) * (T_f__ - T_sky__) / (T_f__ - T_amb__);
        // h_free_conv
        double T_film__ = (T_f__ + T_amb__) / 2;
        double miu_film__ = cfdsolver::calc_fluid_prop("miu", T_film__, 1.01325);
        double rho_film__ = cfdsolver::calc_fluid_prop("rho", T_film__, 1.01325);
        double viu_film__ = miu_film__ / rho_film__;
        double k_film__ = cfdsolver::calc_fluid_prop("k", T_film__, 1.01325);
        double alpha_film__ = cfdsolver::calc_fluid_prop("alpha", T_film__, 1.01325);
        double charl__ = pow(scheme_ref.mesh.size["area"][col], 0.5);
        double Ral__ = 9.81 * (T_f__ - T_amb__) * pow(charl__, 3) / (T_film__ * viu_film__ * alpha_film__);
        double NuN__ = 0.0;
        if(Ral__ < pow(10, 7))
        {
            NuN__ = 0.54 * pow(Ral__, 0.25);
        }
        else
        {
            NuN__ = 0.15 * pow(Ral__, 1/3);
        };
        double h_conv__ = NuN__ * k_film__ / charl__;
        double g_hamb__ = gamma_f_[col] * Sf_.axes_to_val(row, col) / dCf_.axes_to_val(row, col);
        this->lhs_fc[domain__].coeffRef(row, col) = g_hamb__ * (1 - (g_hamb__ / (g_hamb__ + (h_sky__ + h_conv__) * Sf_.axes_to_val(row, col))));
    }
    else if(check.compare("out") == 0)
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        // momentum
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = u_ref.value.cgrad["fluid"][row];
        v1_gradc__ = v_ref.value.cgrad["fluid"][row];
        v2_gradc__ = w_ref.value.cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = u_ref.value.cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v_ref.value.cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = w_ref.value.cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(u_ref.axis) = v0_f__; vf__(v_ref.axis) = v1_f__; vf__(w_ref.axis) = v2_f__; 
        this->lhs_fc["fluid"].coeffRef(row, col) = scheme_ref.prop.rho["face"][col] * (vf__.dot(Sf__)); 
    }
    else
    {
        return;
    };
};
void energy::calc_bound_rhs(int row, int col, std::string check, int unique, cfdscheme::scheme& const scheme_ref, std::string domain__, momentum& const u_ref, momentum& const v_ref, momentum& const w_ref)
{
    if(check.compare("temp"))
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"][domain__];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"][domain__];
        make<double>::map_int& gamma_f_ = this->gamma["face"][domain__];
        this->rhs_fc[domain__].coeffRef(row, col) = gamma_f_[col] * Sf_.axes_to_val(row, col) * scheme_ref.source.value[check][unique]/ dCf_.axes_to_val(row, col);
    }
    else if(check.compare("flux"))
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"][domain__];
        this->rhs_fc[domain__].coeffRef(row, col) = scheme_ref.source.value[check][unique] * Sf_.axes_to_val(row, col);
    }
    else if(check.compare("s2s"))
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"][domain__];
        this->rhs_fc[domain__].coeffRef(row, col) = scheme_ref.source.value[check][unique];
    }
    else if(check.compare("hamb"))
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"][domain__];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"][domain__];
        make<double>::map_int& gamma_f_ = this->gamma["face"][domain__];
        // h_sky
        double g_hamb = gamma_f_[col] * Sf_.axes_to_val(row, col) / dCf_.axes_to_val(row, col);
        double T_f__ = this->value.fvalue[domain__][col];
        double T_amb__ = scheme_ref.source.value["hamb"][0]; double T_sky__ = 0.0552 * pow(T_amb__, 1.5);
        double eps_C__ = scheme_ref.prop.eps["cell"][row];
        double h_sky__ = 5.67 * pow(10, -8) * eps_C__ * (T_f__ + T_sky__) * (pow(T_f__, 2) + pow(T_sky__, 2)) * (T_f__ - T_sky__) / (T_f__ - T_amb__);
        // h_free_conv
        double T_film__ = (T_f__ + T_amb__) / 2;
        double miu_film__ = cfdsolver::calc_fluid_prop("miu", T_film__, 1.01325);
        double rho_film__ = cfdsolver::calc_fluid_prop("rho", T_film__, 1.01325);
        double viu_film__ = miu_film__ / rho_film__;
        double k_film__ = cfdsolver::calc_fluid_prop("k", T_film__, 1.01325);
        double alpha_film__ = cfdsolver::calc_fluid_prop("alpha", T_film__, 1.01325);
        double charl__ = pow(scheme_ref.mesh.size["area"][col], 0.5);
        double Ral__ = 9.81 * (T_f__ - T_amb__) * pow(charl__, 3) / (T_film__ * viu_film__ * alpha_film__);
        double NuN__ = 0.0;
        if(Ral__ < pow(10, 7))
        {
            NuN__ = 0.54 * pow(Ral__, 0.25);
        }
        else
        {
            NuN__ = 0.15 * pow(Ral__, 1/3);
        };
        double h_conv__ = NuN__ * k_film__ / charl__;
        double g_hamb__ = gamma_f_[col] * Sf_.axes_to_val(row, col) / dCf_.axes_to_val(row, col);
        this->rhs_fc[domain__].coeffRef(row, col) = (-1) * g_hamb__ * (h_sky__ + h_conv__) * Sf_.axes_to_val(row, col) * T_amb__ / (g_hamb__ + (h_sky__ + h_conv__) * Sf_.axes_to_val(row, col));
    }
    else if(check.compare("out") == 0)
    {
        axes& Sf_ = scheme_ref.mesh.geom["Sf"]["fluid"];
        axes& eCf_ = scheme_ref.mesh.geom["eCf"]["fluid"];
        axes& dCf_ = scheme_ref.mesh.geom["dCf"]["fluid"];
        coor Sf__ = Sf_.axes_to_coor(row, col);
        coor eCf__ = eCf_.axes_to_coor(row, col);
        coor dCf__ = dCf_.axes_to_coor(row, col);
        // momentum
        coor vf__;
        coor v0_gradc__; coor v1_gradc__; coor v2_gradc__;
        coor v0_gradf__; coor v1_gradf__; coor v2_gradf__;
        double v0_f__; double v1_f__; double v2_f__;
        v0_gradc__ = u_ref.value.cgrad["fluid"][row];
        v1_gradc__ = v_ref.value.cgrad["fluid"][row];
        v2_gradc__ = w_ref.value.cgrad["fluid"][row];
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v1_gradf__ = v1_gradc__ - (v1_gradc__.dot(eCf__) * eCf__);
        v0_gradf__ = v0_gradc__ - (v0_gradc__.dot(eCf__) * eCf__);
        v0_f__ = u_ref.value.cvalue["fluid"][row] + (v0_gradf__.dot(dCf__));
        v1_f__ = v_ref.value.cvalue["fluid"][row] + (v1_gradf__.dot(dCf__));
        v2_f__ = w_ref.value.cvalue["fluid"][row] + (v2_gradf__.dot(dCf__));
        vf__(u_ref.axis) = v0_f__; vf__(v_ref.axis) = v1_f__; vf__(w_ref.axis) = v2_f__;
        // turb_e
        coor T_gradc__; coor T_gradf__;
        T_gradc__ = this->value.cgrad["fluid"][row];
        T_gradf__ = T_gradc__ - (T_gradc__.dot(eCf__) * eCf__);
        this->rhs_fc["fluid"].coeffRef(row, col) = (-1) * scheme_ref.prop.rho["face"][col] * (vf__.dot(Sf__)) * (T_gradf__.dot(dCf__)); 
    }
    else
    {
        return;
    };
};
// s2s
void s2s::calc_lhs(cfdscheme::scheme& scheme_ref)
{
    make<double>::sp_mat& cc_ = this->lhs_cc["s2s"];
    make<make<int>::vec>::map_int& s2s_f_l_ = scheme_ref.mesh.fid["s2s"];
    make<double>::sp_mat& view_ = scheme_ref.mesh.constants["view"]["s2s"];
    make<double>::map_int& rho_f_ = scheme_ref.prop.rho["face"];
    make<double>::map_int& area_ = scheme_ref.mesh.size["area"];
    double rho_f__; double area_clust__; make<int>::vec s2s_f_vec__;
    int row;
    for(int i = 0; i < cc_.outerSize(); i++)
    {
        rho_f__ = 0.0; area_clust__ = 0.0;
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(cc_, i); it; it++)
        {
            cc_.coeffRef(it.row(), it.row()) = 1;
            s2s_f_vec__ = s2s_f_l_[it.row()];
            for(auto j = s2s_f_vec__.begin(); j != s2s_f_vec__.end(); j++)
            {
                rho_f__ += rho_f_[*j] * area_[*j];
                area_clust__ += area_[*j];
            };
            rho_f__ = rho_f__ / area_clust__;
            cc_.coeffRef(it.row(), it.col()) = rho_f__ * view_.coeffRef(it.row(), it.col());
        };
    };
};
void s2s::calc_rhs(cfdscheme::scheme& scheme_ref, energy& const energy_ref)
{
    make<double>::sp_mat& cc_ = this->rhs_cc["fluid"];
    make<make<int>::vec>::map_int& s2s_f_l_ = scheme_ref.mesh.fid["s2s"];
    make<double>::sp_mat& view_ = scheme_ref.mesh.constants["view"]["s2s"];
    make<double>::map_int& eps_f_ = scheme_ref.prop.eps["face"];
    make<double>::map_int& area_ = scheme_ref.mesh.size["area"];
    double eps_f__; double T_f__; double area_clust__; make<int>::vec s2s_f_vec__;
    int row;
    for(int i = 0; i < cc_.outerSize(); i++)
    {
        eps_f__ = 0.0; T_f__ = 0.0; area_clust__ = 0.0;
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(cc_, i); it; it++)
        {
            cc_.coeffRef(it.row(), it.row()) = 1;
            s2s_f_vec__ = s2s_f_l_[it.row()];
            for(auto j = s2s_f_vec__.begin(); j != s2s_f_vec__.end(); j++)
            {
                eps_f__ += eps_f_[*j] * area_[*j];
                T_f__ += energy_ref.value.fvalue["fluid"][*j] * area_[*j];
                area_clust__ += area_[*j];
            };
            eps_f__ = eps_f__ / area_clust__; T_f__ = T_f__ / area_clust__;
            cc_.coeffRef(it.row(), it.col()) = eps_f__ * 5.67 * pow(10, -8) * pow(T_f__, 4);
        };
    };
};
};
#endif