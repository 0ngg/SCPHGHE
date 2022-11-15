/*
cfdlinear
    Just build for SCPHGHE use case dumbass, dont complicate things just yet.
*/

#ifndef CFDLINEAR_H
#define CFDLINEAR_H

#include<iostream>
#include<vector>
#include<map>
#include<string>
#include<utility> // pair
#include<algorithm>
#include<memory>
#include<numeric>
#include<functional>
#include<Eigen/Dense>
#include"cfdgeom.h"

namespace cfdlinear
{

typedef unsigned int uint;
typedef Eigen::Vector3d coor;

struct clin // cell info for linearizing
{
    // neigh cell id, coor value
    std::map<uint, coor> csf;
    std::map<uint, coor> cef;
    std::map<uint, double> cdCF;

    // neigh face id, coor value
    std::map<uint, coor> fsf;
    std::map<uint, double> fdCf;
};

struct sclin // surface cluster for s2s
{
    std::vector<uint> cmember;
    double scF; // view factor
    coor scsf;
};

struct prop
{
    double k;
    double cp;
    double rho;
    double miu;
    double alpha;
};

struct coef
{
    Eigen::MatrixXd* a_mat;
    Eigen::MatrixXd* b_mat;
};

// coef value build with class scheme

void gen_init()
{



};

class scheme
{
    public:
    std::unique_ptr<std::map<std::string, coef>> lin_coef; // pointer matrix 0 -> a, 1 -> b
    std::unique_ptr<std::vector<std::vector<prop>>> lin_prop; // always current iter values
    
    scheme (cfdgeom::mesh minfo, Eigen::MatrixXd* Pinits, Eigen::MatrixXd* Tinits,
            std::vector<std::string> lphys, bool if_turb, bool if_s2s)
    {
        // lphys specify if momentum and/or energy; u, v, z, T, P
        // if_turb specify if turbulence is on (k-e and RANS); turb_k, turb_e
        // if_s2s specify if s2s is on; s2s
        // pressure correction is always on
        // constructor with c_mat equal sets of initial values

        // initialize vars
        // lin_cell
        gen_lincell(lin_cell, minfo);

        // lin_sc
        gen_linsc(lin_sc, minfo);

        // lin_prop
        update_prop(lin_prop, Pinits, Tinits);

        // new a_mat, b_mat
        for(auto i = lphys.begin(); i != lphys.end(); i++)
        {
            switch(*i)
            {
                case "momentum":
                {   
                    Eigen::MatrixXd* x_mom_a = new Eigen::MatrixXd;
                    Eigen::MatrixXd* y_mom_a = new Eigen::MatrixXd;
                    Eigen::MatrixXd* z_mom_a = new Eigen::MatrixXd;
                    Eigen::MatrixXd* p_a = new Eigen::MatrixXd;

                    Eigen::MatrixXd* x_mom_b = new Eigen::MatrixXd;
                    Eigen::MatrixXd* y_mom_b = new Eigen::MatrixXd;
                    Eigen::MatrixXd* z_mom_b = new Eigen::MatrixXd;
                    Eigen::MatrixXd* p_b = new Eigen::MatrixXd;

                    coef xtemp;
                    xtemp.a_mat = x_mom_a;
                    xtemp.b_mat = x_mom_b;
                    *mat.insert({"u", xtemp});

                    coef ytemp;
                    ytemp.a_mat = y_mom_a;
                    ytemp.b_mat = y_mom_b;
                    *mat.insert({"v", ytemp});

                    coef ztemp;
                    ztemp.a_mat = z_mom_a;
                    ztemp.b_mat = z_mom_b;
                    *mat.insert({"2", ztemp});

                    coef Ptemp;
                    Ptemp.a_mat = p_a;
                    Ptemp.b_mat = p_b;
                    *mat.insert({"P", Ptemp});
                }; // momentum

                case "energy":
                {
                    Eigen::MatrixXd* en_a = new Eigen::MatrixXd;
                    Eigen::MatrixXd* en_b = new Eigen::MatrixXd;

                    coef entemp;
                    entemp.a_mat = en_a;
                    entemp.b_mat = en_b;
                    *mat.insert({"T", entemp});
                };

                default:
                cout << "please specify variables to solve. Aborting..." << endl;
                return;
            }; // switch
        }; // for

        if(if_turb == true)
        {
            Eigen::MatrixXd* k_a = new Eigen::MatrixXd*;
            Eigen::MatrixXd* e_a = new Eigen::MatrixXd*;

            Eigen::MatrixXd* k_b = new Eigen::MatrixXd*;
            Eigen::MatrixXd* e_b = new Eigen::MatrixXd*;

            coef ktemp;
            ktemp.a_mat = k_a;
            ktemp.b_mat = k_b;
            *mat.insert({"k", ktemp});

            coef etemp;
            etemp.a_mat = e_a;
            etemp.b_mat = e_b;
            *mat.insert({"e", etemp});
        };

        if(if_s2s == true)
        {
            Eigen::MatrixXd* s2s_a = new Eigen::MatrixXd*;
            Eigen::MatrixXd* s2s_b = new Eigen::MatrixXd*;

            coef rtemp;
            rtemp.a_mat = s2s_a;
            rtemp.b_mat = s2s_b;
            *mat.insert({"s2s", rtemp});
        };

        // initial coefs
        update_coef(a_mat, b_mat, lin_cell, lin_prop);

    }; // constructor

    // P and T values
    void update_prop(std::vector<std::vector<prop>>*, Eigen::MatrixXd*, Eigen::MatrixXd*);
    
    void update_coef(std::map<std::string, Eigen::MatrixXd*>*, std::map<std::string,
                     Eigen::MatrixXd*>*, std::vector<std::vector<clin>>*,
                     std::vector<std::vector<prop>>*); 
    
    protected:
    std::unique_ptr<std::vector<std::vector<clin>>> lin_cell;
    std::unique_ptr<std::vector<std::vector<sclin>>> lin_sc; // for s2s

    private:
    void gen_lincell(std::vector<std::vector<clin>>*, cfdgeom::mesh);
    void gen_linsc(std::vector<std::vector<sclin>>*, cfdgeom::mesh);

}; // scheme

void scheme::gen_lincell()
{


};

void scheme::gen_linsc()
{


};

void scheme::update_prop(lprop, Pinit, Tinit)
{



}; // update_prop

void scheme::update_coef(lcoef, lclin, lprop)
{
    for(std::pair<K,V> entry: *lcoef)
    {
        switch(entry->first)
        {
            case "v": // self matrix, turb_k matrix, P matrix
            momentum::calc_y_mom(entry->second, lclin, lprop);
            
            case "en": // 
            energy::calc_en(entry->second, lclin, lprop);

            case "s2s":
            energy::calc_s2s(entry->second, lclin, lprop);

            case "k":
            turbulence::calc_k(entry->second, lclin, lprop);

            default:
            momentum::calc_gen_mom(entry->first, lclin, lprop)


        }; // switch
    }; // for


}; // update_coef

class scheme : private momentum
{
    // x y z mom
    // diff, conv, src is always on

    // axis, self coef matrix, P value matrix, turb_k value matrix
    void calc_gen_mom   (std::string, coef, Eigen::MatrixXd*, Eigen::MatrixXd*,
                         std::vector<std::vector<clin>>*, std::vector<std::vector<prop>>*);

    // self matrix, 
    void calc_y_mom(coef, std::vector<std::vector<clin>>*, std::vector<std::vector<prop>>*);

}; // momentum

class scheme : protected p_correct
{
    // diff, conv, src is always on
    void calc_p_correct(coef, std::vector<std::vector<clin>>*, std::vector<std::vector<prop>>*);


};

class scheme : private energy
{
    // diff, src is always on
    // conv is specified
    void calc_en(coef, std::vector<std::vector<clin>>*, std::vector<std::vector<prop>>*);
    void calc_s2s(coef, std::vector<std::vector<clin>>*, std::vector<std::vector<prop>>*);
    
    private:
    void calc_soil();

}; // energy

class scheme : private turbulence
{
    // diff, conv, src is always on
    void calc_turb_k(coef, std::vector<std::vector<clin>>*, std::vector<std::vector<prop>>*);
    void calc_turb_e(coef, std::vector<std::vector<clin>>*, std::vector<std::vector<prop>>*);

}; // turbulence


}; // namespace cfdlinear

#endif