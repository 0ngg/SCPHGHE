#include"temp_util_solver.h"

#ifdef SOLVERCFD_H
namespace cfdsolver
{
class user
{
    public:
    make<make<make<double>::map_str>::map_int>::map_str face_init; // face name - unique - u, v, w, T, k, e
    make<make<make<double>::map_str>::map_str>::map_str cell_init; // cell domain - domain name - u, v, w, T, k, e
    make<make<make<double>::map_int>::map_int>::map_str face_source; // .... - time step - value
    make<make<make<double>::map_int>::map_str>::map_str cell_source; // .... - time step - value
    user() {};
    private:
    void read_init_csv();
    void read_source_csv();
    void read_solid_prop_csv();
    void update_source();
};
template <class V>
class solver
{
    public:
    double step_length;
    double under_relax;
    int max_iter;
    double min_residual;
    int current_time;
    make<double>::sp_mat& lhs;
    make<double>::sp_mat& rhs;
    make<double>::sp_mat& prev_lhs;
    make<double>::sp_mat& prev_rhs;
    V eq;
    solver() {};
    private:
    void make_solver();
    void make_init_lhs();
    void make_init_rhs();
    void make_transient_lhs();
    void make_transient_rhs();
    void update_values();
    void initialize();
    void iterate();
};

};
#endif