#include"temp_struct_solver.h"

#ifdef SOLVERCFD_H
namespace cfdsolver
{
// user func
void user::read_init_csv()
{

};
void user::read_source_csv()
{

};
void user::read_solid_prop_csv()
{

};
void user::update_source(int current_time, cfdscheme::scheme& scheme_ref)
{
    for(std::pair<std::string, make<double>::map_int> entry_face : scheme_ref.source.face_value)
    {
        for(std::pair<int, double> entry_face_int : entry_face.second)
        {
            entry_face_int.second = this->face_source[entry_face.first][entry_face_int.first][current_time];
        };
    };
    for(std::pair<std::string, make<double>::map_str> entry_cell : scheme_ref.source.cell_value)
    {
        for(std::pair<std::string, double> entry_cell_str : entry_cell.second)
        {
            entry_cell_str.second = this->cell_source[entry_cell.first][entry_cell_str.first][current_time];
        };
    };
};
// solver func
template <class V>
void make_prev_lhs(solver<V> solv__, cfdscheme::scheme& scheme_ref, bool is_init, bool is_s2s)
{
    // include if s2s karena beda sendiri (kali area, not volume)
    if(is_init)
    {
        if(is_s2s)
        {
            
        }
        else
        {

        };
    }
    else
    {
        if(is_s2s)
        {

        };
        else
        {

        };
    };
};
template <class V>
void make_prev_rhs(solver<V> solv__, cfdscheme::scheme& scheme_ref, bool is_init, bool is_s2s)
{
    // include if s2s karena beda sendiri (kali area, not volume)
    if(is_init)
    {
        if(is_s2s)
        {

        }
        else
        {

        };
    }
    else
    {
        if(is_s2s)
        {
            
        };
        else
        {

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
    Eigen::VectorXf rhs__ = solv__.rhs.row(0);
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
void scphghe::iterate(cfdscheme::scheme& scheme_ref, user& user_ref, bool is_init)
{
    // momentum - pcorrect loop, loop 1
    // turb_k - turb_e loop, loop 2
    // loop 1 - loop 2 loop, loop 12
    // s2s - energy loop, loop 3
    // loop 12 - loop 3, loop 123
    cfdlinear::momentum& u_ = this->solv_u.eq;
    cfdlinear::momentum& v_ = this->solv_v.eq;
    cfdlinear::momentum& w_ = this->solv_w.eq;
    cfdlinear::pcorrect& pcor_ = this->solv_pcor.eq;
    cfdlinear::turb_k& k_ = this->solv_k.eq;
    cfdlinear::turb_e& e_ = this->solv_e.eq;
    cfdlinear::energy& energy_ = this->solv_energy.eq;
    cfdlinear::s2s& s2s_ = this->solv_s2s.eq;
    // loop 1, anchor check convergence -> pcor
    int ctrl = 0;
    int passes = 0;
    while(ctrl < this->max_iter)
    {
        passes += 1;
        // pcor
        make<make<double>::map_int>::map_str prev_iter_cvalue = pcor_.value.cvalue;
        pcor_.update_linear(scheme_ref, u_, v_, w_);
        make_transient_lhs(this->solv_pcor, this->under_relax);
        make_transient_rhs(this->solv_pcor, this->under_relax);
        update_values(this->solv_pcor, scheme_ref);
        this->solv_pcor.eq.update_correction(scheme_ref, u_, v_, w_);
        make<make<double>::map_int>::map_str current_iter_cvalue = pcor_.value.cvalue;
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
        if(check_convergence(current_iter_cvalue, prev_iter_cvalue, this->min_residual))
        {
            passes += 1;
            // loop 2, anchor check convergence -> turb_k
            make<make<double>::map_int>::map_str prev_iter_cvalue = k_.value.cvalue;
            k_.update_linear(scheme_ref, e_, u_, v_, w_);
            make_transient_lhs(this->solv_k, this->under_relax);
            make_transient_rhs(this->solv_k, this->under_relax);
            update_values(this->solv_k, scheme_ref);
            make<make<double>::map_int>::map_str current_iter_cvalue = k_.value.cvalue;
            e_.update_linear(scheme_ref, k_, u_, v_, w_);
            make_transient_lhs(this->solv_e, this->under_relax);
            make_transient_rhs(this->solv_e, this->under_relax);
            update_values(this->solv_e, scheme_ref);
            if(check_convergence(current_iter_cvalue, prev_iter_cvalue, this->min_residual))
            {
                passes += 1;
                // loop 3, anchor check convergence -> energy
                s2s_.update_linear(scheme_ref, energy_);
                make_transient_lhs(this->solv_s2s, this->under_relax);
                make_transient_rhs(this->solv_s2s, this->under_relax);
                update_values(this->solv_s2s, scheme_ref);
                s2s_.update_source_s2s(scheme_ref);
                make<make<double>::map_int>::map_str prev_iter_cvalue = energy_.value.cvalue;
                energy_.update_linear(scheme_ref, u_, v_, w_);
                make_transient_lhs(this->solv_energy, this->under_relax);
                make_transient_rhs(this->solv_energy, this->under_relax);
                update_values(this->solv_energy, scheme_ref);
                make<make<double>::map_int>::map_str current_iter_cvalue = energy_.value.cvalue;
                if(check_convergence(current_iter_cvalue, prev_iter_cvalue, this->min_residual))
                {
                    update_fluid_prop(scheme_ref, energy_);
                    update_wall(scheme_ref, u_, v_, w_, k_, e_, energy_);
                    ctrl += 1;
                    if(passes < 5) // max. 4 passes
                    {
                        make_prev_rhs<cfdlinear::momentum>(this->solv_u, scheme_ref, is_init, false);
                        make_prev_rhs<cfdlinear::momentum>(this->solv_v, scheme_ref, is_init, false);
                        make_prev_rhs<cfdlinear::momentum>(this->solv_w, scheme_ref, is_init, false);
                        make_prev_rhs<cfdlinear::pcorrect>(this->solv_pcor, scheme_ref, is_init, false);
                        make_prev_rhs<cfdlinear::turb_k>(this->solv_k, scheme_ref, is_init, false);
                        make_prev_rhs<cfdlinear::turb_e>(this->solv_e, scheme_ref, is_init, false);
                        make_prev_rhs<cfdlinear::energy>(this->solv_energy, scheme_ref, is_init, false);
                        make_prev_rhs<cfdlinear::s2s>(this->solv_s2s, scheme_ref, is_init, true);
                        this->current_time += 1;
                        user_ref.update_source(this->current_time, scheme_ref);
                        export_converge();
                        return;
                    }
                    else
                    {
                        // update lhs prev transient (aC)
                        make_prev_lhs<cfdlinear::momentum>(this->solv_u, scheme_ref, is_init, false);
                        make_prev_lhs<cfdlinear::momentum>(this->solv_v, scheme_ref, is_init, false);
                        make_prev_lhs<cfdlinear::momentum>(this->solv_w, scheme_ref, is_init, false);
                        make_prev_lhs<cfdlinear::pcorrect>(this->solv_pcor, scheme_ref, is_init, false);
                        make_prev_lhs<cfdlinear::turb_k>(this->solv_k, scheme_ref, is_init, false);
                        make_prev_lhs<cfdlinear::turb_e>(this->solv_e, scheme_ref, is_init, false);
                        make_prev_lhs<cfdlinear::energy>(this->solv_energy, scheme_ref, is_init, false);
                        make_prev_lhs<cfdlinear::s2s>(this->solv_s2s, scheme_ref, is_init, true);
                    };
                };
            };
        };
    };
};
};
#endif