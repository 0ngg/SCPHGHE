#ifndef CFDLINEAR_H
#define CFDLINEAR_H

#include"cfdutil.h"

namespace cfd
{
// mesh
struct geom
{
    make<std::string>::vec names;
    make<uint>::vec nodes;
    coor centroid;
    double size;
    geom() {};
    geom(make<uint>::vec& nodes_in, coor& centroid_in, double& size_in,
         make<std::string>::vec& names_in): nodes(nodes_in), centroid(centroid_in),
         size(size_in), names(names_in) {};
};
struct finfo : geom
{
    coor normal;
    finfo(make<uint>::vec& nodes_in, coor& centroid_in, double& size_in,
          make<std::string>::vec& names_in, coor& normal_in): normal(normal_in)
    {
        geom(nodes_in, centroid_in, size_in, names_in);
    };
    finfo& operator()(const finfo& other)
    {
        this->nodes = other.nodes;
        this->centroid = other.centroid;
        this->size = other.size;
        this->names = other.names;
        this->normal = other.normal;
    };
};
struct cinfo : geom
{
    make<uint>::vec faces;
    make<std::pair<uint, uint>>::vec neigh_c;
    make<uint>::vec neigh_b;
    bool is_wall;
    coor parallel;
    cinfo(make<uint>::vec& nodes_in, coor& centroid_in, double& size_in,
          make<std::string>::vec& names_in, make<uint>::vec& faces_in,
          make<std::pair<uint, uint>>::vec& neigh_c_in,
          make<uint>::vec& neigh_b_in, bool& is_wall_in, coor& parallel_in):
          faces(faces_in), neigh_c(neigh_c_in), neigh_b(neigh_b_in), is_wall(is_wall_in),
          parallel(parallel_in)
    {
        geom(nodes_in, centroid_in, size_in, names_in);
    };
    cinfo& operator()(const cinfo& other)
    {
        this->nodes = other.nodes;
        this->faces = other.faces;
        this->neigh_c = other.neigh_c;
        this->neigh_b = other.neigh_b;
        this->is_wall = other.is_wall;
        this->centroid = other.centroid;
        this->size = other.size;
        this->names = other.names;
    };
};
struct scinfo
{
    make<uint>::vec faces;
    make<uint>::vec size;
    scinfo() {};
    scinfo(make<uint>::vec& faces_in, double& size_in): faces(faces_in),
           size(size_in) {};
};
struct minfo
{
    axes Sf;
    axes Ef;
    axes Tf;
    axes eCf;
    axes eCF;
    axes dCf;
    axes dCF;
    minfo() {};
    minfo(axes Sf_in, axes Ef_in, axes Tf_in, axes eCf_in, axes eCF_in, axes dCf_in,
          axes dCF_in): Sf(Sf_in), Ef(Ef_in), Tf(Tf_in), eCf(eCf_in), eCF(eCF_in),
          dCf(dCf_in), dCF(dCF_in) {};
};
class mesh
{
    public:
    make<coor>::vec nodes;
    make<finfo>::vec faces;
    make<cinfo>::vec cells;
    make<scinfo>::vec clusts;
    minfo info;
    make<double>::sp_mat view;
    mesh() {};
    mesh(make<coor>::vec& nodes_in, make<finfo>::vec& faces_in,
         make<cinfo>::vec& cells_in, make<scinfo>::vec& clusts_in,
         minfo& info_in, make<double>::sp_mat& view_in): nodes(nodes_in), faces(faces_in),
         cells(cells_in), clusts(clusts_in), info(info_in), view(view_in) {};
};
// scheme
struct vinfo
{
    make<double>::sp_mat cvalue;
    make<double>::sp_mat fvalue;
    make<double>::sp_mat prev_cvalue;
    make<double>::sp_mat prev_fvalue;
    axes cgrad;
    axes fgrad;
    axes prev_cgrad;
    axes prev_fgrad;
    vinfo() {};
    vinfo(make<double>::sp_mat& cvalue_in, make<double>::sp_mat& fvalue_in,
          make<double>::sp_mat& prev_cvalue_in, make<double>::sp_mat& prev_fvalue_in,
          axes& cgrad_in, axes& fgrad_in, axes& prev_cgrad_in, axes& prev_fgrad_in):
          cvalue(cvalue_in), fvalue(fvalue_in), prev_cvalue(prev_cvalue_in),
          prev_fvalue(prev_fvalue_in), cgrad(cgrad_in), fgrad(fgrad_in),
          prev_cgrad(prev_cgrad_in), prev_fgrad(prev_fgrad_in) {};
};
struct ginfo
{
    make<double>::sp_mat g_conv;
    make<double>::sp_mat g_diff;
    make<double>::sp_mat g;
    make<double>::sp_mat g_src_conv;
    make<double>::sp_mat g_src_diff;
    make<double>::sp_mat g_src;
    make<uint>::sp_mat cc_fc;
    make<double>::sp_mat rho_v_sf;
    make<double>::sp_mat prop_f;
    ginfo() {};
    ginfo(make<double>::sp_mat& g_diff_in, make<double>::sp_mat& g_conv_in,
          make<double>::sp_mat& g_src_diff_in, make<double>::sp_mat& g_src_conv_in,
          make<uint>::sp_mat& cc_fc_in, make<double>::sp_mat& rho_v_sf_in,
          make<double>::sp_mat& prop_f_in): g_diff(g_diff_in), g_conv(g_conv_in),
          g_src_diff(g_src_diff_in), g_src_conv(g_src_conv_in), cc_fc(cc_fc_in),
          rho_v_sf(rho_v_sf_in), prop_f(prop_f_in) {};
};
template <class U>
class scheme
{
    public:
    U eq;
    vinfo value;
    ginfo g;
    scheme() {};
    scheme(U& eq_in, vinfo& val_in, ginfo& g_in): eq(eq_in), value(val_in), g(g_in) {};
};
// solver
class solver
{
    public:
    make<double>::sp_mat lhs;
    make<double>::sp_mat rhs;
    double res = 0.0;
    solver() {};
    solver(make<double>::sp_mat& lhs_in, make<double>::sp_mat& rhs_in) : lhs(lhs_in),
           rhs(rhs_in)
    {
        this->res = 0.0;
    };
};
}
#endif