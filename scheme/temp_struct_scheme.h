#include"temp_util.h"
#include"temp_util_scheme.h"

#ifdef SCHEMECFD_H
namespace cfdscheme
{
struct finfo
{
    make<int>::vec fnode;
    coor fcentroid;
    coor fnormal;
    coor fparallel;
    double farea;
};
struct cinfo
{
    make<int>::vec cnode;
    make<int>::vec cface;
    coor ccentroid;
    double cvolume;
};
class minfo
{  
    public:
    make<make<axes>::map_str>::map_str geom; // Sf, Ef, Tf, eCf, eCF, dCf, dCF, parallel = -Sf.norm() from fluid perspective
    make<make<double>::vec>::map_str size; // area, volume
    make<make<make<int>::vec>::map_int>::map_str fid; // store common (ns, inlet, outlet, conj) index 0, temp, flux, and s2s index unique id
    make<make<make<int>::vec>::map_str>::map_str cid; // store domain name unique id (derived from fluid/solid)
    make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int bid; // store cell-face-<boundary name, unique-constant/iter value>
    make<make<int>::sp_mat>::map_str cc_fc; // store neighbor cell to face alias id;
    make<make<double>::sp_mat>::map_str fc; // store face-cell sp_mat templates for fluid, solid, and conj (energy use)
    make<make<double>::sp_mat>::map_str cc; // store cell-cell sp_mat templates for fluid, solid, conj, and s2s (energy use)
    make<make<make<double>::sp_mat>::map_str>::map_str constants; // g_conv_aC, g_conv_aF, g_diff_aC, g_diff_aF, gc, view (cc, s2s)
    minfo() {};
    minfo(mshio::MshSpec spec)
    {
        this->make_minfo(spec);
    };
    private:
    void make_minfo(mshio::MshSpec);
    void make_id(mshio::MshSpec, make<cinfo>::map_int&);
    void make_template(make<coor>::map_int&, make<finfo>::map_int&, make<cinfo>::map_int&);
    void make_size(make<coor>::map_int&, make<finfo>::map_int&, make<cinfo>::map_int&);
    void make_geom(make<finfo>::map_int&, make<cinfo>::map_int&);
    void make_constants();
    void make_clust(make<finfo>::map_int&);
};
class pinfo
{
    public:
    make<make<double>::map_int>::map_str rho; // cell, fluid
    make<make<double>::map_int>::map_str miu; // cell, fluid 
    make<make<double>::map_int>::map_str cp; // cell, fluid
    make<make<double>::map_int>::map_str k; // cell, solid
    make<double>::map_int eps; // cell, s2s
    make<double>::map_int alpha; // cell, s2s
    pinfo() {};
    pinfo(minfo& mesh_)
    {
        this->make_pinfo(mesh_);
    };
    private:
    void make_pinfo(minfo&);
};
class winfo
{
    public:
    make<double>::map_int ts; // cell, fluid
    make<double>::map_int utau; // cell, fluid
    make<double>::map_int miut; // cell, fluid
    winfo() {};
    winfo(minfo& mesh_)
    {
        this->make_winfo(mesh_);
    };
    private:
    void make_winfo(minfo&);
};
struct vinfo
{
    make<make<double>::map_int>::map_str value; // cc
    make<make<double>::map_int>::map_str prev_value; // cc
    make<make<coor>::map_int>::map_str grad; // cc
    make<make<coor>::map_int>::map_str prev_grad; // cc
    vinfo() {};
    vinfo(make<std::string>::vec which, minfo& mesh_)
    {
        this->make_vinfo(which, mesh_);
    };
    private:
    void make_vinfo(make<std::string>::vec, minfo&);
};
struct binfo
{
    /*momentum, pcorrect, turb_k, and turb_e*/ // -> fluid only
    // inlet value, environment velocity magnitude
    /*energy*/
    // isotemp(temp) value (unique), environment/ambient temperature
    // isoflux(flux) value (unique), environment irr q
    // s2s(s2s) value (unique), s2s irr q from solver
    // hamb(hamb) value, ambient temperature
    make<make<double>::map_int>::map_str value; // boundary name-unique-constant/iter value
    binfo() {};
    binfo(make<make<double>::map_int>::map_str& value_in): value(value_in) {};
};
class scheme
{
    public:
    minfo mesh;
    pinfo prop;
    winfo wall;
    vinfo pressure;
    binfo source;
    scheme() {};
    scheme(std::string msh_file)
    {
        mshio::MshSpec spec = mshio::load_msh(msh_file);
        minfo mesh_(spec); pinfo prop_(mesh_); winfo wall_(mesh_);
        vinfo pressure_(make<std::string>::vec{"fluid"}, mesh_);
        this->mesh = mesh_;
        this->prop = prop_;
        this->wall = wall_;
        this->pressure = pressure_;
    };
};
};
#endif