#include"temp_util_export.h"

#ifdef CFDEXPORT_H

template <class V>
struct compexport
{
    typedef make<make<V>::vec>::map_str comp_str;
};
class export
{
    public:
    // common information
    std::string mesh_name;
    int number_of_cells;
    /* map_str names: P, u, v, w, k, e, T, s2s */
    // converged values
    make<compexport<double>::comp_str>::map_int cell_export;
    make<compexport<double>::comp_str>::map_int face_export;
    /* map_str names: u, v, w, k, e, T, s2s */
    // final numerical errors (at converged iteration)
    // BiCGSTAB estimated numerical errors
    compexport<double>::comp_str err_export;
    // converged residuals for each time step of each variable of interest
    compexport<double>::comp_str res_export;
    // time till convergence for each time step
    make<long long>::vec time_export;
    // number of overarching loops (all, SIMPLE-turbulence-energy) needed till all equation converged
    make<int>::vec iter_export;
    export() {};
    export(std::string filename, cfdscheme::scheme& scheme_ref)
    {
        this->mesh_name = filename;
        int n_cells = 0;
        compexport<double>::comp_str empty_err_res;
        empty_err_res.insert({"u", make<double>::vec()});
        empty_err_res.insert({"v", make<double>::vec()});
        empty_err_res.insert({"w", make<double>::vec()});
        empty_err_res.insert({"k", make<double>::vec()});
        empty_err_res.insert({"e", make<double>::vec()});
        empty_err_res.insert({"T", make<double>::vec()});
        empty_err_res.insert({"s2s", make<double>::vec()});
        empty_err_res.insert({"pcor", make<double>::vec()});
        compexport<double>::comp_str empty_converged(empty_err_res);
        empty_converged.insert({"P", make<double>::vec()});       
        for(std::pair<std::string, make<make<int>::vec>::map_str> entry1 : scheme_ref.mesh.cid)
        {
            for(std::pair<std::string, make<int>::vec> entry2 : entry1.second)
            {
                for(auto i = entry2.second.begin(); i != entry2.second.end(); i++)
                {
                    n_cells += 1;
                    this->cell_export.insert({*i, empty_converged});
                    this->face_export.insert({*i, empty_converged});
                };
            };
        };
        this->number_of_cells = n_cells;
        this->err_export = compexport<double>::comp_str(empty_err_res);
        this->res_export = compexport<double>::comp_str(empty_err_res);
        this->time_export = make<long long>::vec();
        this->iter_export = make<int>::vec();
    };
    void update_export(cfdsolver::scphghe, cfdscheme::scheme&, make<double>::map_str,
                       make<double>::map_str, long long, int);
    private:
    void export_to_sql();
    // void export_to_vtk();
};

#endif