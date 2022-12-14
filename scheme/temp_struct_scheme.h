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
    make<coor>::map_int nodes;
    make<finfo>::map_int faces;
    make<cinfo>::map_int cells;
    make<make<axes>::map_str>::map_str geom; // Sf, Ef, Tf, eCf, eCF, dCf, dCF, parallel = -Sf.norm() from fluid perspective
    make<make<double>::map_int>::map_str size; // area, volume
    make<make<make<int>::vec>::map_int>::map_str fid; // store common (ns, inlet, outlet, conj) index 0, temp, flux, and s2s index unique id
    make<make<make<int>::vec>::map_str>::map_str cid; // store domain name unique id (derived from fluid/solid)
    make<make<make<make<std::pair<std::string, int>>::vec>::map_int>::map_int>::map_str bid; // store cell-face-<boundary name, unique-constant/iter value>
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
    make<make<double>::map_int>::map_str rho; // cell-face, fluid
    make<make<double>::map_int>::map_str miu; // cell-face, fluid 
    make<make<double>::map_int>::map_str cp; // cell-face, fluid
    make<make<double>::map_int>::map_str k; // cell-face, solid -> constant
    make<make<double>::map_int>::map_str eps; // face, solid(s2s -> glass/absorber, hamb -> glass/soil ids -> inline! (0.8)) -> constant
    pinfo() {};
    pinfo(minfo& mesh_, user& user_ref)
    {
        this->make_pinfo(mesh_, user_ref);
    };
    private:
    void make_pinfo(minfo&, user&);
};
class winfo
{
    public:
    make<double>::map_int wall_dist; // dCf wall
    make<coor>::map_int wall_parallel; // cell parallel transformator
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
    make<make<double>::map_int>::map_str cvalue;
    make<make<double>::map_int>::map_str fvalue;
    make<make<double>::map_int>::map_str prev_cvalue;
    make<make<double>::map_int>::map_str prev_fvalue;
    make<make<coor>::map_int>::map_str cgrad;
    make<make<coor>::map_int>::map_str fgrad;
    make<make<coor>::map_int>::map_str prev_cgrad;
    make<make<coor>::map_int>::map_str prev_fgrad;
    vinfo() {};
    vinfo(make<std::string>::vec which, minfo& mesh_, double init_value)
    {
        this->make_vinfo(which, mesh_, init_value);
    };
    private:
    void make_vinfo(make<std::string>::vec, minfo&, double);
};
struct binfo
{
    /*momentum, pcorrect, turb_k, and turb_e*/ // -> fluid only
    // inlet value, environment velocity magnitude
    /*energy*/
    // isotemp(temp) value (unique), environment/ambient temperature
    // isoflux(flux) value (unique), environment irr q
    // s2s(flux) value (unique), s2s solver q // create inline program
    // hamb(hamb) value, ambient temperature
    make<make<double>::map_int>::map_str face_value; // boundary name-unique-constant/iter value
    make<make<double>::map_str>::map_str cell_value;
    binfo() {};
    binfo(make<make<double>::map_int>::map_str& face_value_in, make<make<double>::map_str>::map_str& cell_value_in):
          face_value(face_value_in), cell_value(cell_value_in) {};
};
class scheme
{
    public:
    minfo mesh;
    pinfo prop;
    winfo wall;
    vinfo pressure;
    binfo source;
    make<double>::sp_mat rho_v_sf; // fluid, fc
    make<double>::map_int phi_v; // fluid, f
    scheme() {};
    scheme(std::string msh_file, user& user_ref)
    {
        mshio::MshSpec spec = mshio::load_msh(msh_file);
        minfo mesh_(spec); pinfo prop_(mesh_, user_ref); winfo wall_(mesh_);
        vinfo pressure_(make<std::string>::vec{"fluid"}, mesh_, user_ref.P_init);
        make<make<double>::map_int>::map_str face_value__;
        for(std::pair<std::string, make<make<double>::map_int>::map_int> entry_face1 : user_ref.face_source)
        {
            make<make<double>::map_int>::map_str::iterator it = face_value__.find(entry_face1.first);
            if(it == face_value__.end())
            {
                face_value__.insert({entry_face1.first, make<double>::map_int()});
            };
            for(std::pair<int, make<double>::map_int> entry_face2 : entry_face1.second)
            {
                face_value__[entry_face1.first].insert({entry_face2.first, entry_face2.second[0]});                  
            };
        };
        make<make<double>::map_str>::map_str cell_value__;
        for(std::pair<std::string, make<make<double>::map_int>::map_str> entry_cell1 : user_ref.cell_source)
        {
            make<make<double>::map_str>::map_str::iterator it = cell_value__.find(entry_cell1.first);
            if(it == cell_value__.end())
            {
                cell_value__.insert({entry_cell1.first, make<double>::map_str()});
            };
            for(std::pair<std::string, make<double>::map_int> entry_cell2 : entry_cell1.second)
            {
                cell_value__[entry_cell1.first].insert({entry_cell2.first, entry_cell2.second[0]});                  
            };
        };
        this->mesh = mesh_;
        this->prop = prop_;
        this->wall = wall_;
        this->pressure = pressure_;
        this->source = binfo(face_value__, cell_value__);
    };
};
};
#endif