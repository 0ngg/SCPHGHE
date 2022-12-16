#ifndef SOLVERCFD_H
#define SOLVERCFD_H

#include"temp_util.h"
#include"temp_struct_linear.h"

namespace cfdsolver
{
template <class V>
struct interpolate
{
    V linear_face_itr(V valC__, V valF__, double gc__, bool is_gamma)
    {
        if(is_gamma)
        {
            return (gc__ * valF__ + (1 - gc__) * valC__) / (valC__ * valF__);
        }
        else
        {
            return gc__ * valC__ + (1 - gc__) * valF__;
        };
    };
    double quick_face_itr(double valC__, coor gradC__, coor gradf__, coor dCf__)
    {
        coor grad__ = gradC__ + gradf__;
        return valC__ + gradf__.dot(dCf__) / 2;
    };
    coor quick_grad_itr(double valC__, double valF__, coor gradC__, coor gradF__, double gc__, double coor dCF__, coor eCF__)
    {
        double valdCF__ = sqrt_sum(dCF__);
        coor gradf__ = interpolate::linear_face_itr(gradC__, gradF__, gc__, false);
        return gradf__ + ((valF__ - valC__) / valdCF__) * eCF__ - (gradf__.dot(eCF__) * eCF__);
    };
};
double calc_fluid_prop(std::string what, double P, double T)
{
    // rho, miu, cp
};
// scheme object update functions outside cfdscheme namespace
void update_wall(cfdscheme::scheme& const scheme_ref, cfdlinear::momentum& const u_ref, cfdlinear::momentum& const v_ref, cfdlinear::momentum& const w_ref,
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
void update_fluid_prop(cfdscheme::scheme& const scheme_ref, cfdlinear::energy& const energy_ref)
{
    cfdscheme::vinfo& pressure_ = scheme_ref.pressure;
    double pval__; double Tval__;
    for(std::pair<std::string, make<double>::map_int> entry_str : scheme_ref.prop.rho)
    {
        for(std::pair<int, double> entry_int : entry_str.second)
        {
            if(entry_str.first.compare("cell") == 0)
            {
                pval__ = pressure_.cvalue["fluid"][entry_int.first];
                Tval__ = energy_ref.value.cvalue["fluid"][entry_int.first];
            }
            else
            {
                pval__ = pressure_.fvalue["fluid"][entry_int.first];
                Tval__ = energy_ref.value.fvalue["fluid"][entry_int.first];
            };
            scheme_ref.prop.rho[entry_str.first][entry_int.first] = calc_fluid_prop("rho", pval__, Tval__);
            scheme_ref.prop.miu[entry_str.first][entry_int.first] = calc_fluid_prop("miu", pval__, Tval__);
            scheme_ref.prop.cp[entry_str.first][entry_int.first] = calc_fluid_prop("cp", pval__, Tval__);
        };
    };
};
};
#endif