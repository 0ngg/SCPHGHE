#include"temp_util_export.h"

#ifdef CFDEXPORT_H

struct cfexport
{
    // cell data export
    int id;
    std::string name; // domain or boundary
    // converged values
    make<double>::vec P;
    make<double>::vec u;
    make<double>::vec v;
    make<double>::vec w;
    make<double>::vec T;
    make<double>::vec k_turb;
    make<double>::vec e_turb;
};
template <typename V>
struct compexport
{
    // computation report (estimated error, residuals,
    // time at (till) converged iteration, overarching iteration
    // until all equation converged)
    make<V>::vec u;
    make<V>::vec v;
    make<V>::vec w;
    make<V>::vec T;
    make<V>::vec k_turb;
    make<V>::vec e_turb;
};
class export
{
    make<cfexport>::map_int cell_export;
    make<cfexport>::map_int face_export;
    // final numerical errors (at converged iteration)
    // BiCGSTAB estimated numerical errors
    make<compexport<double>>::map_str err_export;
    // converged residuals for each time step of each variable of interest
    make<compexport<double>>::map_str res_export;
    // time till convergence for each time step
    make<compexport<double>>::map_int time_export;
    // number of overarching loops (SIMPLE-turbulence-energy) needed till all equation converged
    make<compexport<int>>::map_int iter_export;
};
// TO DO
// due dilligence
// grid convergence study https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html
// analytical errors (python) VARIABELNYA APA AJA YG DICARI?!
// temporal convergence gausah wkwk

#endif