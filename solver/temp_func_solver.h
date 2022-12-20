#include"temp_struct_solver.h"
#include"temp_struct_user.h"

#ifdef SOLVERCFD_H
namespace cfdsolver
{
// solver func
template <class V>
void make_prev_lhs(solver<V> solv__, cfdscheme::scheme& scheme_ref, int step_length__, bool is_init)
{
    make<double>::map_int& rho_c_ = scheme_ref.prop.rho["cell"];
    make<double>::map_int& vol_ = scheme_ref.mesh.size["volume"];
    if(is_init)
    {
        // rho * vol / step length
        for(int i = 0; i < solv__.lhs.outerSize(); i++)
        {
            for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(solv__.lhs.outerSize(), i); it; it++)
            {
                if(it.row() == it.col())
                {
                    solv__.lhs.coeffRef(it.row(), it.row()) = rho_c_[it.row()] * vol_[it.row()] / step_length__;
                };
            };
        };
    }
    else
    {
        // rho * vol / (2*step length)
        for(int i = 0; i < solv__.lhs.outerSize(); i++)
        {
            for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(solv__.lhs.outerSize(), i); it; it++)
            {
                if(it.row() == it.col())
                {
                    solv__.lhs.coeffRef(it.row(), it.row()) = rho_c_[it.row()] * vol_[it.row()] / (2 * step_length__);
                };
            };
        };
    };
};
template <class V>
void make_prev_rhs(solver<V> solv__, cfdscheme::scheme& scheme_ref, int step_length__, bool is_init)
{
    make<double>::map_int& rho_c_ = scheme_ref.prop.rho["cell"];
    make<double>::map_int& vol_ = scheme_ref.mesh.size["volume"];
    make<double>::map_int& cvalue_ = solv__.eq.value.cvalue["fluid"];
    if(is_init)
    {
        // rho * vol * converged value / step length
        for(int i = 0; i < solv__.rhs.outerSize(); i++)
        {
            for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(solv__.rhs.outerSize(), i); it; it++)
            {
                make<double>::map_int::iterator it_cell = cvalue_.find(it.row());
                if(it_cell != cvalue_.end())
                {
                    solv__.rhs.coeffRef(it.row(), 0) = rho_c_[it.row()] * vol_[it.row()] * cvalue_[it.row()] / step_length__;
                };
                break;
            };
        };
    }
    else
    {
        // rho * vol * converged value / (2*step length)
        for(int i = 0; i < solv__.rhs.outerSize(); i++)
        {
            for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(solv__.rhs.outerSize(), i); it; it++)
            {
                make<double>::map_int::iterator it_cell = cvalue_.find(it.row());
                if(it_cell != cvalue_.end())
                {
                    solv__.rhs.coeffRef(it.row(), 0) = rho_c_[it.row()] * vol_[it.row()] * cvalue_[it.row()] / (2 * step_length__);
                };
                break;
            };
        };
    };
};
template <class V>
void make_transient_lhs(solver<V> solv__, double under_relax__)
{
    int ctd = 0;
    make<make<double>::sp_mat>::std::string& lhs_cc_ = solv__.eq.lhs_cc;
    make<double>::sp_mat& lhs_ = solv__.lhs;
    make<double>::sp_mat& prev_lhs_ = solv__.prev_lhs;
    for(std::pair<std::string, make<double>::sp_mat> entry_cc : lhs_cc_)
    {
        if(ctd == 0)
        {
            lhs_ = entry_cc.second + prev_lhs_;
        }
        else
        {
            lhs_ += entry_cc.second;
        };
        ctd += 1;
    };
    // under-relaxation
    for(int i = 0; i < lhs_.outerSize(); i++)
    {
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(lhs_, i); it; it++)
        {
            if(it.row() == it.col())
            {
                lhs_.coeffRef(it.row(), it.row()) = it.value() / under_relax__;
            };
        };
    };
};
template <class V>
void make_transient_rhs(solver<V> solv__, double under_relax__)
{
    int ctd = 0;
    make<make<double>::sp_mat>::std::string& rhs_cc_ = solv__.eq.rhs_cc;
    make<double>::sp_mat& rhs_ = solv__.rhs;
    make<double>::sp_mat& prev_rhs_ = solv__.prev_rhs;
    for(std::pair<std::string, make<double>::sp_mat> entry_cc : eq__.rhs_cc)
    {
        if(ctd == 0)
        {
            rhs_ = entry_cc.second + prev_rhs_;
        }
        else
        {
            rhs_ += entry_cc.second;
        };
        ctd += 1;
    };
    // under-relaxation
    for(std::pair<std::string, make<double>::sp_mat> entry_rhs : rhs_cc_)
    {
        if(entry_rhs.first.compare("conj") != 0)
        {
            for(int i = 0; i < entry_rhs.second.outerSize(); i++)
            {
                for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(entry_rhs.second, i); it; it++)
                {
                    if(it.row() == it.col())
                    {
                        rhs_.coeffRef(it.row(), 0) += (1 - under_relax__) * rhs_.coeffRef(it.row(), it.row()) *
                                                        solv__.eq.value.cvalue[entry_rhs.first][it.row()] / under_relax__;

                    };
                };
            };
        };
    };
}
template <class V>
Eigen::VectorXf get_new_values(solver<V> solv__)
{
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, RowMajor>> solver;
    make<double>::sp_mat lhs__ = solv__.lhs;
    Eigen::VectorXf rhs__ = solv__.rhs.column(0);
    solver.compute(lhs__);
    Eigen::VectorXf x = solver.solve(rhs__);
    return x;
};
template <class V>
void update_values(solver<V> solv__, cfdscheme::scheme& scheme_ref)
{
    interpolate<double> inter;
    Eigen::VectorXf new_cvalues__ = get_new_values(solv__);
    // i.e. gradient computation
    // update cvalue -> cgrad (by iteration) -> fgrad -> fvalue
    for(std::pair<std::string, make<double>::sp_mat> entry : solv__.eq.lhs_cc)
    {
        make<double>::map_int& cvalue_ = eq__.value.cvalue[entry.first];
        make<int>::sp_mat& cc_fc_ = scheme_ref.mesh.cc_fc[entry.first];
        for(std::pair<int, double> entry : eq__.value.cvalue[entry.first])
        {
            entry.second = new_cvalues__(entry.first);
        };
        inter.least_square_itr(eq__, scheme_ref, entry.first);
        for(int i = 0; i < cc_fc_.outerSize(); i++)
        {
            for(make<int>::sp_mat::InnerIterator it(cc_fc_, i); it; it++)
            {
                eq__.value.fgrad[entry.first][it.value()] = inter.quick_face_itr(it.row(), it.col(), eq__, scheme_, "grad");
                eq__.value.fvalue[entry.first][it.value()] = inter.quick_face_itr(it.row(), it.col(), eq__, scheme_, "face");
            };
        };
    };
};
// scphghe
void scphghe::make_scphghe(cfdscheme::scheme& const scheme_ref, double step_length_in, double under_relax_in, int max_iter_in,
                           double min_residual_in)
{
    this->step_length = step_length_in;
    this->under_relax = under_relax_in;
    this->max_iter = max_iter_in;
    this->min_residual = min_residual_in;
    this->current_time = 0;
    cfdlinear::momentum u_ = cfdlinear::momentum(scheme_ref, 0);
    this->solv_u = solver(u_, scheme_ref, step_length_in);
    cfdlinear::momentum v_ = cfdlinear::momentum(scheme_ref, 1);
    this->solv_v = solver(v_, scheme_ref, step_length_in);
    cfdlinear::momentum w_ = cfdlinear::momentum(scheme_ref, 2);
    this->solv_w = solver(w_, scheme_ref, step_length_in);
    cfdlinear::pcorrect pcor_ = cfdlinear::pcorrect(scheme_ref);
    this->solv_pcor = solver(pcor_, scheme_ref, step_length_in);
    cfdlinear::turb_k k_turb_ = cfdlinear::turb_k(scheme_ref);
    this->solv_k = solver(k_turb_, scheme_ref, step_length_in);
    cfdlinear::turb_e e_turb_ = cfdlinear::turb_e(scheme_ref);
    this->solv_e = solver(e_turb_, scheme_ref, step_length_in);
    cfdlinear::energy energy_ = cfdlinear::energy(scheme_ref);
    this->solv_energy = solver(energy_, scheme_ref, step_length_in);
    cfdlinear::s2s s2s_ = cfdlinear::s2s(scheme_ref);
    this->solv_s2s = solver(s2s_, scheme_ref, step_length_in);
};
int scphghe::SIMPLE_loop(cfdscheme::scheme& scheme_ref, bool is_init)
{
    // momentum - pcorrect loop
    cfdlinear::momentum& u_ = this->solv_u.eq;
    cfdlinear::momentum& v_ = this->solv_v.eq;
    cfdlinear::momentum& w_ = this->solv_w.eq;
    cfdlinear::turb_k& k_ = this->solv_k.eq;
    cfdlinear::pcorrect& pcor_ = this->solv_pcor.eq;
    int ctrl = 0;
    int passes = 0;
    while(ctrl < this->max_iter)
    {
        passes += 1;
        // momentum
        u_.update_linear(scheme_ref, k_, v_, w_);
        make_transient_lhs(this->solv_u, this->under_relax);
        make_transient_rhs(this->solv_u, this->under_relax);
        v_.update_linear(scheme_ref, k_, u_, w_);
        make_transient_lhs(this->solv_v, this->under_relax);
        make_transient_rhs(this->solv_v, this->under_relax);
        w_.update_linear(scheme_ref, k_, u_, v_);
        make_transient_lhs(this->solv_w, this->under_relax);
        make_transient_rhs(this->solv_w, this->under_relax);
        update_values(this->solv_u, scheme_ref);
        update_values(this->solv_v, scheme_ref);
        update_values(this->solv_w, scheme_ref);
        u_.calc_wall(scheme_ref, v_, w_);
        v_.calc_wall(scheme_ref, u_, w_);
        w_.calc_wall(scheme_ref, u_, v_);
        // pcor
        pcor_.update_linear(scheme_ref, u_, v_, w_);
        make_transient_lhs(this->solv_pcor, this->under_relax);
        make_transient_rhs(this->solv_pcor, this->under_relax);
        update_values(this->solv_pcor, scheme_ref);
        this->solv_pcor.eq.update_correction(scheme_ref, u_, v_, w_);
        // new cvalue, old lhs - rhs
        if(check_convergence(this->solv_u, this->min_residual) &&
           check_convergence(this->solv_v, this->min_residual) &&
           check_convergence(this->solv_w, this->min_residual))
        {
            if(passes < 2)
            {
                make_prev_lhs<cfdlinear::momentum>(this->solv_u, scheme_ref, this->step_length, is_init);
                make_prev_lhs<cfdlinear::momentum>(this->solv_v, scheme_ref, this->step_length, is_init);
                make_prev_lhs<cfdlinear::momentum>(this->solv_w, scheme_ref, this->step_length, is_init);
                make_prev_lhs<cfdlinear::pcorrect>(this->solv_pcor, scheme_ref, this->step_length, is_init);
                make_prev_rhs<cfdlinear::momentum>(this->solv_u, scheme_ref, this->step_length, is_init);
                make_prev_rhs<cfdlinear::momentum>(this->solv_v, scheme_ref, this->step_length, is_init);
                make_prev_rhs<cfdlinear::momentum>(this->solv_w, scheme_ref, this->step_length, is_init);
                make_prev_rhs<cfdlinear::pcorrect>(this->solv_pcor, scheme_ref, this->step_length, is_init);
                return ctrl;
            }
            else
            {
                passes = 0;
                make_prev_lhs<cfdlinear::momentum>(this->solv_u, scheme_ref, this->step_length, is_init);
                make_prev_lhs<cfdlinear::momentum>(this->solv_v, scheme_ref, this->step_length, is_init);
                make_prev_lhs<cfdlinear::momentum>(this->solv_w, scheme_ref, this->step_length, is_init);
                make_prev_lhs<cfdlinear::pcorrect>(this->solv_pcor, scheme_ref, this->step_length, is_init);
                ctrl += 1;
            };
        };
    };
};
int scphghe::turb_loop(cfdscheme::scheme& scheme_ref, bool is_init)
{
    cfdlinear::momentum& u_ = this->solv_u.eq;
    cfdlinear::momentum& v_ = this->solv_v.eq;
    cfdlinear::momentum& w_ = this->solv_w.eq;
    cfdlinear::turb_k& k_ = this->solv_k.eq;
    cfdlinear::turb_e& e_ = this->solv_e.eq;
    int ctrl = 0;
    int passes = 0;
    while (ctrl < this->max_iter)
    {
        passes += 1;
        // turb_k
        k_.update_linear(scheme_ref, e_, u_, v_, w_);
        make_transient_lhs(this->solv_k, this->under_relax);
        make_transient_rhs(this->solv_k, this->under_relax);
        update_values(this->solv_k, scheme_ref);
        k_.calc_wall(scheme_ref, e_);
        // turb_e
        e_.update_linear(scheme_ref, k_, u_, v_, w_);
        make_transient_lhs(this->solv_e, this->under_relax);
        make_transient_rhs(this->solv_e, this->under_relax);
        update_values(this->solv_e, scheme_ref);
        e_.calc_wall(scheme_ref);
        if(check_convergence(this->solv_k, this->min_residual) &&
           check_convergence(this->solv_e, this->min_residual))
        {
            if(passes < 2)
            {
                make_prev_lhs<cfdlinear::turb_k>(this->solv_k, scheme_ref, this->step_length, is_init);
                make_prev_lhs<cfdlinear::turb_e>(this->solv_e, scheme_ref, this->step_length, is_init);
                make_prev_rhs<cfdlinear::turb_k>(this->solv_k, scheme_ref, this->step_length, is_init);
                make_prev_rhs<cfdlinear::turb_e>(this->solv_e, scheme_ref, this->step_length, is_init);
                return ctrl;
            }
            else
            {
                passes = 0;
                make_prev_lhs<cfdlinear::turb_k>(this->solv_k, scheme_ref, this->step_length, is_init);
                make_prev_lhs<cfdlinear::turb_e>(this->solv_e, scheme_ref, this->step_length, is_init);
                ctrl += 1;
            };
        };
    };
};
int scphghe::energy_loop(cfdscheme::scheme& scheme_ref, double AH, bool is_init)
{
    cfdlinear::momentum& u_ = this->solv_u.eq;
    cfdlinear::momentum& v_ = this->solv_v.eq;
    cfdlinear::momentum& w_ = this->solv_w.eq;
    cfdlinear::energy& energy_ = this->solv_energy.eq;
    cfdlinear::s2s& s2s_ = this->solv_s2s.eq;
    int ctrl = 0;
    int passes = 0;
    while(ctrl < this->max_iter)
    {
        passes += 1;
        // s2s
        s2s_.update_linear(scheme_ref, energy_);
        make_transient_lhs(this->solv_s2s, this->under_relax);
        make_transient_rhs(this->solv_s2s, this->under_relax);
        update_values(this->solv_s2s, scheme_ref);
        s2s_.update_source_s2s(scheme_ref);
        // energy
        energy_.update_linear(scheme_ref, u_, v_, w_, AH);
        make_transient_lhs(this->solv_energy, this->under_relax);
        make_transient_rhs(this->solv_energy, this->under_relax);
        update_values(this->solv_energy, scheme_ref);
        energy_.calc_wall(scheme_ref, AH);
        if(check_convergence(this->solv_energy, this->min_residual) &&
           check_convergence(this->solv_s2s, this->min_residual))
        {
            if(passes < 2)
            {
                update_fluid_prop(scheme_ref, energy_, AH);
                make_prev_lhs<cfdlinear::energy>(this->solv_energy, scheme_ref, this->step_length, is_init);
                make_prev_rhs<cfdlinear::energy>(this->solv_energy, scheme_ref, this->step_length, is_init);
                return ctrl;
            }
            else
            {
                passes = 0;
                make_prev_lhs<cfdlinear::energy>(this->solv_energy, scheme_ref, this->step_length, is_init);
                ctrl += 1;
            };
        };
    };
};
void scphghe::iterate(cfdscheme::scheme& scheme_ref, user& user_ref, bool is_init)
{
    int ctrl_iter = 0;
    int ctrl_inner = 0;
    user_ref.update_source(this->current_time, scheme_ref);
    double AH = this->absolute_humidity;
    while(ctrl_iter < this->max_iter)
    {
        ctrl_inner += SIMPLE_loop(scheme_ref, is_init);
        ctrl_inner += turb_loop(scheme_ref, is_init);
        ctrl_inner += energy_loop(scheme_ref, AH, is_init);
        if(ctrl_inner < 4)
        {
            // export solution at time t
            this->current_time += 1;
        }
        else
        {
            ctrl_iter += 1;
        };
    };
};
};
#endif