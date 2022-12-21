#include"temp_util_solver.h"
#include"temp_struct_user.h"

#ifdef SOLVERCFD_H
namespace cfdsolver
{
template <class V>
struct solver
{
    V eq;
    make<double>::sp_mat lhs;
    make<double>::sp_mat rhs;
    make<double>::sp_mat prev_lhs;
    make<double>::sp_mat prev_rhs;
    solver() {};
    solver(V eq, cfdscheme::scheme& scheme_ref, int step_length)
    {
        this->eq = eq;
        make<double>::sp_mat lhs__; make<double>::sp_mat rhs__;
        make<double>::sp_mat prev_lhs__; make<double>::sp_mat prev_rhs__;
        make<sparse_input>::vec prev_lhs_input__; make<sparse_input>::vec prev_rhs_input__;
        make<double>::map_int rho_c_ = scheme_ref.prop.rho["cell"];
        make<double>::map_int vol_ = scheme_ref.mesh.size["volume"];
        int ctd = 0;
        for(std::pair<std::string, make<double>::sp_mat> entry : this->eq.lhs_cc)
        {
            if(ctd == 0)
            {
                lhs__ = entry.second;
                rhs__ = this->eq.rhs_cc[entry.first];
                if(entry.first.compare("fluid") == 0)
                {
                    for(int i = 0; i < lhs__.outerSize(); i++)
                    {
                        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(lhs__, i); it; it++)
                        {
                            if(it.row() == it.col())
                            {
                                double rho_c__ = rho_c_[it.row()];
                                double vol__ = vol_[it.row()];
                                sparse_input temp_lhs__(it.row(), it.row(), rho_c__ * vol__ / step_length);
                                sparse_input temp_rhs__(it.row(), 0, rho_c__ * vol__ * eq.value.cvalue["fluid"][it.row()]/ step_length);
                                prev_lhs_input__.push_back(temp_lhs__); prev_rhs_input__.push_back(temp_rhs__);
                            };
                        };
                    };
                };
            }
            else
            {
                lhs__ += entry.second;
                rhs__ += this->eq.rhs_cc[entry.first];
                if(entry.first.compare("fluid") == 0)
                {
                    for(int i = 0; i < lhs__.outerSize(); i++)
                    {
                        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(lhs__, i); it; it++)
                        {
                            if(it.row() == it.col())
                            {
                                double rho_c__ = rho_c_[it.row()];
                                double vol__ = vol_[it.row()];
                                sparse_input temp_lhs__(it.row(), it.row(), rho_c__ * vol__ / step_length);
                                sparse_input temp_rhs__(it.row(), 0, rho_c__ * vol__ * eq.value.cvalue["fluid"][it.row()]/ step_length);
                                prev_lhs_input__.push_back(temp_lhs__); prev_rhs_input__.push_back(temp_rhs__);
                            };
                        };
                    };
                };
            };
            ctd += 1;
        };
        prev_lhs__.setFromTriplets(prev_lhs_input__.begin(), prev_lhs__input__.end());
        prev_rhs__.setFromTriplets(prev_rhs_input__.begin(), prev_rhs_input__.end());
        this->lhs = lhs__; this->rhs = rhs__;
        this->prev_lhs = prev_lhs__; this->prev_rhs = prev_rhs__;
    };
};
class scphghe
{
    public:
    double step_length;
    double under_relax;
    double min_residual;
    int max_iter;
    int current_time;
    solver<cfdlinear::momentum> solv_u;
    solver<cfdlinear::momentum> solv_v;
    solver<cfdlinear::momentum> solv_w;
    solver<cfdlinear::pcorrect> solv_pcor;
    solver<cfdlinear::turb_k> solv_k;
    solver<cfdlinear::turb_e> solv_e;
    solver<cfdlinear::energy> solv_energy;
    solver<cfdlinear::s2s> solv_s2s;
    scphghe() {};
    private:
    void make_scphghe(cfdscheme::scheme&, user&, double, double, int, double);
    int SIMPLE_loop(cfdscheme::scheme&, bool);
    int turb_loop(cfdscheme::scheme&, bool);
    int energy_loop(cfdscheme::scheme&, double, bool);
    void iterate(cfdscheme::scheme&, user&, bool);
};
};
#endif