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
void user::update_source()
{
    // exclude s2s since already hard coded into its update_linear
    // marching time updates

};
// solver func
template <class V>
void solver<V>::make_solver()
{

};
template <class V>
void solver<V>::make_init_lhs()
{

};
template <class V>
void solver<V>::make_init_rhs()
{

};
template <class V>
void solver<V>::make_transient_lhs()
{
    solver& self = *this;
    make<double>::sp_mat lhs__;
    int ctd = 0;
    for(std::pair<std::string, make<double>::sp_mat> entry_cc : this->eq.lhs_cc)
    {
        if(ctd == 0)
        {
            this->lhs = entry_cc.second + this->prev_lhs;
        }
        else
        {
            this->lhs += entry_cc.second;
        };
        ctd += 1;
    };
    // under-relaxation
    for(int i = 0; i < this->lhs.outerSize(); i++)
    {
        for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(this->lhs, i); it; it++)
        {
            if(it.row() == it.col())
            {
                this->lhs.coeffRef(it.row(), it.row()) = it.value() / this->under_relax;
            };
        };
    };
};
template <class V>
void solver<V>::make_transient_rhs()
{
    make<double>::sp_mat rhs__;
    int ctd = 0;
    for(std::pair<std::string, make<double>::sp_mat> entry_cc : eq.rhs_cc)
    {
        if(ctd == 0)
        {
            this->rhs = entry_cc.second + this->prev_rhs;
        }
        else
        {
            this->rhs += entry_cc.second;
        };
        ctd += 1;
    };
    // under-relaxation
    for(std::pair<std::string, make<double>::sp_mat> entry_lhs : this->eq.lhs_cc)
    {
        if(entry_lhs.first.compare("conj") != 0)
        {
            for(int i = 0; i < entry_lhs.second.outerSize(); i++)
            {
                for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(entry_lhs.second, i); it; it++)
                {
                    if(it.row() == it.col())
                    {
                        this->rhs.coeffRef(it.row(), 0) += (1 - this->under_relax) * this->lhs.coeffRef(it.row(), it.row()) *
                                                           this->eq.value.cvalue[entry_lhs.first][it.row()] / this->under_relax;

                    };
                };
            };
        };
    };
}
template <class V>
void solver<V>::update_values()
{
    
};
template <class V>
void solver<V>::initialize()
{

};
template <class V>
void solver<V>::iterate()
{

};


};
#endif