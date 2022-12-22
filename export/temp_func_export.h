#include"temp_struct_export.h"

#ifdef CFDEXPORT_H

void export::update_export(cfdsolver::scphghe scphghe_ref, cfdscheme::scheme& scheme_ref, make<double>::map_str err,
                           make<double>::map_str res, long long time_iter, int outer_iter)
{
    cfdlinear::momentum& u_ = scphghe_ref.solv_u.eq;
    cfdlinear::momentum& v_ = scphghe_ref.solv_v.eq;
    cfdlinear::momentum& w_ = scphghe_ref.solv_w.eq;
    cfdlinear::turb_k& k_ = scphghe_ref.solv_k.eq;
    cfdlinear::turb_e& e_ = scphghe_ref.solv_e.eq;
    cfdlinear::energy& energy_ = scphghe_ref.solv_energy.eq;
    make<compexport<double>::comp_str>::map_int& cell_export_ = this->cell_export;
    make<compexport<double>::comp_str>::map_int& face_export_ = this->face_export;
    // converged cell values
    for(std::pair<int, double> entry : scheme_ref.pressure.cvalue["fluid"])
    {
        cell_export_[entry.first]["P"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : u_.value.cvalue["fluid"])
    {
        cell_export_[entry.first]["u"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : v_.value.cvalue["fluid"])
    {
        cell_export_[entry.first]["v"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : w_.value.cvalue["fluid"])
    {
        cell_export_[entry.first]["w"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : k_.value.cvalue["fluid"])
    {
        cell_export_[entry.first]["k"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : e_.value.cvalue["fluid"])
    {
        cell_export_[entry.first]["e"].push_back(entry.second);
    };
    for(std::pair<std::string, make<double>::map_int> entry_str : energy_.value.cvalue)
    {
        for(std::pair<int, double> entry_int : entry_str.second)
        {
            cell_export_[entry_int.first]["T"].push_back(entry_int.second);
        };
    };
    // converged face values
    for(std::pair<int, double> entry : scheme_ref.pressure.fvalue["fluid"])
    {
        face_export_[entry.first]["P"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : u_.value.fvalue["fluid"])
    {
        face_export_[entry.first]["u"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : v_.value.fvalue["fluid"])
    {
        face_export_[entry.first]["v"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : w_.value.fvalue["fluid"])
    {
        face_export_[entry.first]["w"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : k_.value.fvalue["fluid"])
    {
        face_export_[entry.first]["k"].push_back(entry.second);
    };
    for(std::pair<int, double> entry : e_.value.fvalue["fluid"])
    {
        face_export_[entry.first]["e"].push_back(entry.second);
    };
    for(std::pair<std::string, make<double>::map_int> entry_str : energy_.value.fvalue)
    {
        for(std::pair<int, double> entry_int : entry_str.second)
        {
            face_export_[entry_int.first]["T"].push_back(entry_int.second);
        };
    };
    // BiCGSTAB estimated numerical errors and residuals
    for(std::pair<std::string, double> entry : err)
    {
        this->err_export[entry.first].push_back(entry.second);
    };
    for(std::pair<std::string, double> entry : res)
    {
        this->res_export[entry.first].push_back(entry.second);
    };
    // time and total iter
    this->time_export.push_back(time_iter);
    this->iter_export.push_back(outer_iter);
};
void export::export_to_sql()
{

};
/*
void export::export_to_vtk()
{

};
*/

#endif