#include"temp_util.h"
#include"temp_util_scheme.h"
#include"temp_struct_scheme.h"

#ifdef SCHEMECFD_H
namespace cfdscheme
{
// minfo
struct quad
{
    make<coor>::vec centroid;
    make<double>::vec area;
};
void make_tri(finfo& face_, make<coor>::map_int& const nodes_)
{
    const make<int>::vec& fnode_ = face_.fnode;
    coor centroid(0.0, 0.0, 0.0);
    coor v1 = nodes_[fnode_[1]] - nodes_[fnode_[0]];
    coor v2 = nodes_[fnode_[2]] - nodes_[fnode_[0]];
    coor vcross = v1.cross(v2);
    for(auto i = fnode_.begin(); i != fnode_.end(); i++)
    {
        centroid += nodes_[*i];
    };
    centroid = centroid/3;
    vcross(0) = pow(vcross(0), 2);
    vcross(1) = pow(vcross(1), 2);
    vcross(2) = pow(vcross(2), 2);
    double area = 0.5 * pow(vcross.block(0, 0, 3, 1).sum(), 0.5);
    coor vnorm = v1.cross(v2);
    vnorm.norm();
    face_.fcentroid(centroid);
    face_.fnormal(vnorm);
    face_.farea = area;
};
void make_quad(finfo& face_, make<coor>::map_int& nodes_)
{
    quad fquad;
    make<int>::vec& fnode_ = face_.fnode;
    std::pair<int, int> diag1;
    coor acc; double len; double rhs;
    for(int i = 0; i < 3; i++)
    {
        for(int j = i; j < 4; j++)
        {
            acc = nodes_[fnode_[i]] + nodes_[fnode_[j]];
            rhs = sqrt_sum(acc);
            if(rhs > len)
            {
                len = rhs;
                diag1.first = fnode_[i];
                diag1.second = fnode_[j];
            };
        };
    };
    double area; make<int>::vec diag2;
    for(auto i = fnode_.begin(); i != fnode_.end(); i++)
    {
        if(*i != diag1.first && *i != diag1.second)
        {
            coor centroid(0.0, 0.0, 0.0);
            centroid += nodes_[diag1.first] + nodes_[diag1.second] + nodes_[*i];
            centroid = centroid/3;
            coor v1 = nodes_[diag1.first] - nodes_[*i];
            coor v2 = nodes_[diag1.second] - nodes_[*i];
            coor vcross = v1.cross(v2);
            area = 0.5 * sqrt_sum(vcross);
            fquad.centroid.push_back(centroid);
            fquad.area.push_back(area);
        }
        else
        {
            diag2.push_back(*i);
        };
    };
    area = std::accumulate(fquad.area.begin(), fquad.area.end(), 0.0);
    coor centroid(0.0, 0.0, 0.0);
    for(size_t i = 0; i < sizeof(fquad.centroid); i++)
    {
        centroid += fquad.centroid[i] * fquad.area[i]; 
    };
    centroid = centroid / area;
    coor v1 =  nodes_[0] - nodes_[diag1.first];
    coor v2 =  nodes_[0] - nodes_[diag1.second];
    coor vnorm = v1.cross(v2);
    vnorm.norm();
    face_.fcentroid(centroid);
    face_.fnormal(vnorm);
    face_.farea = area;
};
void make_cell(minfo& mesh_, cinfo& cell_, make<finfo>::map_int& faces_)
{
    coor centroid(0.0, 0.0, 0.0);
    double area = 0.0;
    for(std::pair<int, finfo> entry : faces_)
    {
        centroid += entry.second.fcentroid * entry.second.farea;
        area += entry.second.farea;
    };
    centroid = centroid/area;
    double vol = 0.0; coor parallel;
    for(std::pair<int, finfo> entry : faces_)
    {
        coor norm = entry.second.fnormal;
        coor vinner = entry.second.fcentroid - centroid;
        if(norm.dot(vinner) < 0)
        {
            norm = norm * -1;
        };
        vol += centroid.dot(norm) * entry.second.farea / 3;
    };
    cell_.ccentroid(centroid);
    cell_.cvolume = vol;
};
void minfo::make_id(mshio::MshSpec spec)
{
    // physical names to boundary and domain id
    make<make<make<int>::vec>::map_int>::map_str fid__;
    make<make<make<int>::vec>::map_str>::map_str cid__;
    make<make<make<int>::vec&>::vec>::map_int it_id_tag;
    for(const auto& group : spec.physical_groups)
    {
        make<std::string>::vec parsed__;
        std::stringstream ssin(group.name);
        std::string token;
        make<make<int>::vec&>::vec it_id_l;
        switch(group.dim)
        {
            case 2:
            while (ssin >> token)
            {
                if(is_unique(token))
                {
                    std::string first = token.substr(0, sizeof(token) - 2);
                    char num(token[sizeof(token) - 1]);
                    int second(num);
                    make<make<make<int>::vec>::map_int>::map_str::iterator it_str = fid__.find(first);
                    if(it_str == fid__.end())
                    {
                        fid__.insert({first, make<make<int>::vec>::map_int()});
                        fid__[first].insert({second, make<int>::vec()});
                    }
                    else
                    {
                        make<make<int>::vec>::map_int::iterator it_int = fid__[first].find(second);
                        if(it_int == fid__[first].end())
                        {
                            fid__[first].insert({second, make<int>::vec()});
                        };
                    };
                    it_id_l.push_back(fid__[first][second]);
                }
                else
                {  
                    make<make<make<int>::vec>::map_int>::map_str::iterator it_str = fid__.find(token);
                    if(it_str == fid__.end())
                    {
                        fid__.insert({token, make<make<int>::vec>::map_int()});
                        fid__[token].insert({0, make<int>::vec()});
                    };
                    it_id_l.push_back(fid__[token][0]);
                };
            };
            case 3:
            while (ssin >> token)
            {
                std::string first = token.substr(0, 5);
                std::string second = token.substr(6, sizeof(token) - 6);
                make<make<make<int>::vec>::map_str>::map_str::iterator it_str = cid__.find(token);
                if(it_str == cid__.end())
                {
                    cid__.insert({first, make<make<int>::vec>::map_str()});
                    cid__[first].insert({second, make<int>::vec()});
                }
                else
                {
                    make<make<int>::vec>::map_str::iterator it_str = cid__[first].find(token);
                    if(it_str == cid[first].end())
                    {
                        cid__[first].insert({second, make<int>::vec()});
                    };
                };
                it_id_l.push_back(cid__[first][second]);
            };
            default:
            return;
        };
    };
    fid__.insert({"none", make<make<int>::vec>::map_int()});
    fid__["none"].insert({0, make<int>::vec()});
    int tag; int type; int n;
    for(int i = 0; i < spec.elements.num_entity_blocks; i++)
    {
        mshio::ElementBlock& block = spec.elements.entity_blocks[i];
        tag = block.entity_tag;
        type = block.element_type;
        n = mshio::nodes_per_element(type);
        make<make<make<int>::vec&>::vec>::map_int::iterator it = it_id_tag.find(tag);
        if(it != it_id_tag.end())
        {
            for(auto j = it_id_tag[tag].begin(); j != it_id_tag[tag].end(); j++)
            {
                make<int>::vec& temp_ref = *j;
                for(int k = 0; k < block.num_elements_in_block; k++)
                {
                    temp_ref.push_back(block.data[k*(n+1)]);
                };
            };
        }
        else
        {
            for(int k = 0; k < block.num_elements_in_block; k++)
            {
                fid__["none"][0].push_back(block.data[k*(n+1)]);
            };
        };
    };
    this->fid = fid__; this->cid = cid__;
};
void minfo::make_template(make<coor>::map_int& nodes, make<finfo>::map_int& faces, make<cinfo>::map_int& cells)
{
    make<make<sparse_input>::vec>::map_str fc_input;
    make<make<sparse_input>::vec>::map_str cc_input;
    make<make<sparse_input>::vec>::map_str cc_fc_input;
    // ids
    make<make<make<int>::vec>::map_int>::map_str& fid_ = this->fid;
    make<make<make<int>::vec>::map_str>::map_str& cid_ = this->cid;
    // fluid
    cc_fc_input.insert({"fluid", make<sparse_input>::vec()});
    fc_input.insert({"fluid", make<sparse_input>::vec()});
    cc_input.insert({"fluid", make<sparse_input>::vec()});
    int fneigh;
    for(std::pair<std::string, make<int>::vec> entry_fluid : cid_["fluid"])
    {
        for(int i = 0; i < sizeof(entry_fluid.second) - 1; i++)
        {
            for(auto j = cells[entry_fluid.second[i]].cface.begin(); j != cells[entry_fluid.second[i]].cface.end(); j++)
            {
                fc_input["fluid"].push_back(sparse_input(entry_fluid.second[i], *j, 0.0));
            };
            cc_input["fluid"].push_back(sparse_input(entry_fluid.second[i], entry_fluid.second[i], 0.0));
            for(int j = i; j < sizeof(entry_fluid.second); j++)
            {
                fneigh = is_neighbor(cells[entry_fluid.second[i]].cface, cells[entry_fluid.second[j]].cface);
                if(fneigh >= 0)
                {
                    cc_fc_input["fluid"].push_back(sparse_input(entry_fluid.second[i], entry_fluid.second[j], fneigh));
                    cc_fc_input["fluid"].push_back(sparse_input(entry_fluid.second[j], entry_fluid.second[i], fneigh));
                    cc_input["fluid"].push_back(sparse_input(entry_fluid.second[i], entry_fluid.second[j], 0.0));
                    cc_input["fluid"].push_back(sparse_input(entry_fluid.second[j], entry_fluid.second[i], 0.0));
                };
            };
        };
        int tail = *entry_fluid.second.end();
        for(auto i = cells[tail].cface.begin(); i != cells[tail].cface.end(); i++)
        {
            fc_input["fluid"].push_back(sparse_input(tail, *i, 0.0));
        };
        cc_input["fluid"].push_back(sparse_input(tail, tail, 0.0));
    };
    // solid + conj
    make<make<make<int>::vec>::map_str>::map_str::iterator it_solid = cid_.find("solid");
    if(it_solid != cid_.end())
    {
        // solid
        cc_fc_input.insert({"solid", make<sparse_input>::vec()});
        fc_input.insert({"solid", make<sparse_input>::vec()});
        cc_input.insert({"solid", make<sparse_input>::vec()});
        // conj
        cc_fc_input.insert({"conj", make<sparse_input>::vec()});
        fc_input.insert({"conj", make<sparse_input>::vec()});
        cc_input.insert({"conj", make<sparse_input>::vec()});
        // wall
        cid_.insert({"misc", make<make<int>::vec>::map_str()});
        cid_["misc"].insert({"wall", make<int>::vec()});
        // iter
        make<std::string>::vec done;
        for(std::pair<std::string, make<int>::vec> entry_solid1 : cid_["solid"])
        {
            done.push_back(entry_solid1.first);
            // self connection
            for(std::size_t i = 0; i < sizeof(entry_solid1.second) - 1; i++)
            {
                for(auto k = cells[entry_solid1.second[i]].cface.begin(); k != cells[entry_solid1.second[i]].cface.end(); k++)
                {
                    fc_input["solid"].push_back(sparse_input(entry_solid1.second[i], *k, 0.0));
                };
                for(std::size_t j = i; j < sizeof(entry_solid1.second); j++)
                {
                    fneigh = is_neighbor(cells[entry_solid1.second[i]].cface, cells[entry_solid1.second[j]].cface);
                    if(fneigh >= 0)
                    {
                        cc_fc_input["solid"].push_back(sparse_input(entry_solid1.second[i], entry_solid1.second[j], fneigh));
                        cc_fc_input["solid"].push_back(sparse_input(entry_solid1.second[j], entry_solid1.second[i], fneigh));
                        cc_input["solid"].push_back(sparse_input(entry_solid1.second[i], entry_solid1.second[j], 0.0));
                        cc_input["solid"].push_back(sparse_input(entry_solid1.second[j], entry_solid1.second[i], 0.0));
                    };
                };
                int tail = *entry_solid1.second.end();
                for(auto k = cells[tail].cface.begin(); k != cells[tail].cface.end(); k++)
                {
                    fc_input["solid"].push_back(sparse_input(tail, *k, 0.0));
                };
                cc_input["solid"].push_back(sparse_input(0.0, tail, tail));
            };
            // solid-solid connection
            for(std::pair<std::string, make<int>::vec> entry_solid2 : cid_["solid"])
            {
                make<std::string>::vec::iterator it_done = std::find(done.begin(), done.end(), entry_solid2);
                if(it_done == done.end())
                {
                    for(std::size_t i = 0; i < sizeof(entry_solid1.second); i++)
                    {
                        for(std::size_t j = 0; j < sizeof(entry_solid2.second); j++)
                        {
                            fneigh = is_neighbor(cells[entry_solid1.second[i]].cface, cells[entry_solid2.second[j]].cface);
                            if(fneigh >= 0)
                            {
                                cc_fc_input["solid"].push_back(sparse_input(entry_solid1.second[i], entry_solid2.second[j], fneigh));
                                cc_fc_input["solid"].push_back(sparse_input(entry_solid2.second[j], entry_solid1.second[i], fneigh));
                                fc_input["solid"].push_back(sparse_input(entry_solid1.second[i], fneigh, 0.0));
                                fc_input["solid"].push_back(sparse_input(entry_solid2.second[j], fneigh, 0.0));
                                cc_input["solid"].push_back(sparse_input(entry_solid1.second[i], entry_solid2.second[j], 0.0));
                                cc_input["solid"].push_back(sparse_input(entry_solid2.second[j], entry_solid1.second[i], 0.0));
                            };
                        };
                    };
                };
            };
            // conj connection
            for(std::pair<std::string, make<int>::vec> entry_fluid : cid_["fluid"])
            {
                for(std::size_t i = 0; i < sizeof(entry_solid1.second); i++)
                {
                    for(std::size_t j = 0; j < sizeof(entry_fluid.second); j++)
                    {
                        fneigh = is_neighbor(cells[entry_solid1.second[i]].cface, cells[entry_fluid.second[j]].cface);
                        if(fneigh >= 0)
                        {
                            cc_fc_input["conj"].push_back(sparse_input(entry_solid1.second[i], entry_fluid.second[j], fneigh));
                            cc_fc_input["conj"].push_back(sparse_input(entry_fluid.second[j], entry_solid1.second[i], fneigh));
                            fc_input["conj"].push_back(sparse_input(entry_solid1.second[i], fneigh, 0.0));
                            fc_input["conj"].push_back(sparse_input(entry_fluid.second[j], fneigh, 0.0));
                            cc_input["conj"].push_back(sparse_input(entry_solid1.second[i], entry_fluid.second[j], 0.0));
                            cc_input["conj"].push_back(sparse_input(entry_fluid.second[j], entry_solid1.second[i], 0.0));
                            cid_["misc"]["wall"].push_back(entry_fluid.second[j]);
                        };
                    };
                };
            };
        };
    };
    // template
    make<make<double>::sp_mat>::map_str fc__;
    for(std::pair<std::string, make<sparse_input>::vec> entry : fc_input)
    {
        fc__.insert({entry.first, make<double>::sp_mat()});
        fc__[entry.first].setFromTriplets(entry.second.begin(), entry.second.end());
    };
    this->fc = fc__;
    make<make<double>::sp_mat>::map_str cc__;
    for(std::pair<std::string, make<sparse_input>::vec> entry : cc_input)
    {
        cc__.insert({entry.first, make<double>::sp_mat()});
        cc__[entry.first].setFromTriplets(entry.second.begin(), entry.second.end());
    };
    this->cc = cc__;
    make<make<int>::sp_mat>::map_str cc_fc__;
    for(std::pair<std::string, make<sparse_input>::vec> entry : cc_fc_input)
    {
        cc_fc__.insert({entry.first, make<int>::sp_mat()});
        cc_fc__[entry.first].setFromTriplets(entry.second.begin(), entry.second.end());
    };
    this->cc_fc = cc_fc__;
};
void minfo::make_size(make<coor>::map_int& nodes, make<finfo>::map_int& faces, make<cinfo>::map_int& cells)
{
    for(std::pair<int, finfo> entry : faces)
    {
        finfo& face_ = entry.second;
        switch(sizeof(face_.fnode))
        {
            case 3:
            make_tri(face_, nodes);
            case 4:
            make_quad(face_, nodes);
            default:
            return;
        };
    };
    for(std::pair<int, cinfo> entry : cells)
    {
        minfo& mesh_ = *this;
        cinfo& cell_ = entry.second;
        make_cell(mesh_, cell_, faces);
    };
    make<make<double>::vec>::map_str size__;
    size__.insert({"area", make<double>::vec()}); size__.insert({"volume", make<double>::vec()});
    for(std::pair<int, finfo> entry : faces)
    {
        size__["area"].push_back(entry.second.farea);
    };
    for(std::pair<int, cinfo> entry : cells)
    {
        size__["volume"].push_back(entry.second.cvolume);
    };
    this->size = size__;
};
void minfo::make_geom(make<finfo>::map_int& faces, make<cinfo>::map_int& cells)
{
    // refs
    make<make<double>::sp_mat>::map_str& cc_ = this->cc;
    make<make<double>::sp_mat>::map_str& fc_ = this->fc;
    make<make<int>::sp_mat>::map_str& cc_fc_ = this->cc_fc;
    // cc-centered
    make<axes>::map_str eCF_; coor eCF__;
    make<axes>::map_str dCF_; coor dCF__;
    for(std::pair<std::string, make<double>::sp_mat> entry : cc_)
    {
        eCF_.insert({entry.first, axes()});
        dCF_.insert({entry.first, axes()});
        for(int i = 0; i < cc_[entry.first].outerSize(); i++)
        {
            for(make<double>::sp_mat::InnerIterator it(cc_[entry.first], i); it; it++)
            {
                cinfo& row_current = cells[it.row()];
                cinfo& col_current = cells[it.col()];
                coor dCF__ = row_current.ccentroid - col_current.ccentroid;
                dCF_[entry.first](it.row(), it.col(), dCF__);
                coor eCF__(dCF__.norm());
                eCF_[entry.first](it.row(), it.col(), eCF__);
            };
        };
    };
    // fc-centered
    make<axes>::map_str Sf_; coor Sf__;
    make<axes>::map_str eCf_; coor eCf__;
    make<axes>::map_str dCf_; coor dCf__;
    for(std::pair<std::string, make<double>::sp_mat> entry : fc_)
    {
        Sf_.insert({entry.first, axes()});
        eCf_.insert({entry.first, axes()});
        dCf_.insert({entry.first, axes()});
        for(int i = 0; i < fc_[entry.first].outerSize(); i++)
        {
            for(make<double>::sp_mat::InnerIterator it(fc_[entry.first], i); it; it++)
            {
                cinfo& row_current = cells[it.row()];
                finfo& col_current = faces[it.col()];
                dCf__ = col_current.fcentroid - row_current.ccentroid;
                eCf__(dCf__.norm());
                Sf__ = col_current.fnormal * col_current.farea;
                if(Sf__.dot(eCf__) < 0)
                {
                    Sf__ = Sf__ * (-1);
                };
                Sf_[entry.first](it.row(), it.col(), Sf__);
                eCf_[entry.first](it.row(), it.col(), eCf__);
                dCf_[entry.first](it.row(), it.col(), dCf__);
            };
        };
    };
    // using cc_fc
    make<axes>::map_str Ef_; coor Ef__;
    make<axes>::map_str Tf_; coor Tf__;
    for(std::pair<std::string, make<int>::sp_mat> entry : cc_fc_)
    {
        Ef_.insert({entry.first, axes()});
        Tf_.insert({entry.first, axes()});
        for(int i = 0; i < cc_fc_[entry.first].outerSize(); i++)
        {
            for(make<int>::sp_mat::InnerIterator it(cc_fc_[entry.first], i); it; it++)
            {
                Sf__(Sf_[entry.first].axes_to_coor(it.row(), it.value()));
                eCF__(eCF_[entry.first].axes_to_coor(it.row(), it.col()));
                Ef__(eCF__ * (Sf__.dot(Sf__) / eCF__.dot(Sf__)));
                Tf__(Sf__ - Ef__);
                Ef_[entry.first](it.row(), it.value(), Ef__);
                Tf_[entry.first](it.row(), it.value(), Tf__);
            };
        };
    };
    make<make<axes>::map_str>::map_str geom__;
    geom__.insert({"Sf", make<axes>::map_str(Sf_)}); geom__.insert({"Ef", make<axes>::map_str(Ef_)});
    geom__.insert({"Tf", make<axes>::map_str(Tf_)}); geom__.insert({"eCf", make<axes>::map_str(eCf_)});
    geom__.insert({"eCF", make<axes>::map_str(eCF_)}); geom__.insert({"dCf", make<axes>::map_str(dCf_)});
    geom__.insert({"dCF", make<axes>::map_str(dCF_)});
    this->geom = geom__;
};
void minfo::make_constants()
{
    // refs
    make<make<int>::sp_mat>::map_str& cc_fc_ = this->cc_fc;
    make<make<double>::sp_mat>::map_str& fc_ = this->fc;
    make<make<double>::sp_mat>::map_str& cc_ = this->cc;
    make<axes>::map_str& Sf_ = this->geom["Sf"]; coor Sf__;
    make<axes>::map_str& eCF_ = this->geom["eCF"]; coor eCF__;
    make<axes>::map_str& dCf_ = this->geom["dCf"]; coor dCf__;
    make<axes>::map_str& dCF_ = this->geom["dCF"]; double dCF__;
    // fcg_conv_aC = ( 1 - ( eCF.dot(dCf) / (2*sqrt_sum(dCF)) ) ), fc
    // fcg_conv_aF = (eCF.dot(dCf) / (2*sqrt_sum(dCF))), cc
    // fcg_diff_aC = (eCF.dot(Sf) / sqrt_sum(dCF))
    // fcg_diff_aF = - (eCF.dot(Sf) / sqrt_sum(dCF))
    // gc = sqrt_sum(dFf) / (sqrt_sum(dFf) + sqrt_sum(dCf)) , cc
    make<make<double>::sp_mat>::map_str g_conv_aC_(fc_);
    make<make<double>::sp_mat>::map_str g_conv_aF_(cc_);
    make<make<double>::sp_mat>::map_str g_diff_aC_(fc_);
    make<make<double>::sp_mat>::map_str g_diff_aF_(cc_);
    make<make<double>::sp_mat>::map_str gc_(cc_);
    double g_conv__; double g_diff__; double gc__;
    for(std::pair<std::string, make<int>::sp_mat> entry : cc_fc_)
    {
        for(int i = 0; i < cc_fc_[entry.first].outerSize(); i++)
        {
            for(make<int>::sp_mat::InnerIterator it(cc_fc_[entry.first], i); it; it++)
            {
                eCF__ = eCF_[entry.first].axes_to_coor(it.row(), it.col());
                dCf__ = dCf_[entry.first].axes_to_coor(it.row(), it.value());
                dCF__ = dCF_[entry.first].axes_to_val(it.row(), it.col());
                g_conv__ = eCF__.dot(dCf__) / (2*dCF__);
                g_diff__ = eCF__.dot(Sf__) / dCF__;
                gc__ = dCf_[entry.first].axes_to_val(it.col(), it.value()) /
                       (dCf_[entry.first].axes_to_val(it.col(), it.value()) +
                       dCf_[entry.first].axes_to_val(it.row(), it.value()));
                g_conv_aC_[entry.first].coeffRef(it.row(), it.value()) = 1 - g_conv__;
                g_conv_aF_[entry.first].coeffRef(it.row(), it.col()) = g_conv__;
                g_diff_aC_[entry.first].coeffRef(it.row(), it.value()) = g_diff__;
                g_diff_aF_[entry.first].coeffRef(it.row(), it.col()) = (-1) * g_diff__;
                gc_[entry.first].coeffRef(it.row(), it.col()) = gc__;
            };
        };
    };
    make<make<make<double>::sp_mat>::map_str>::map_str constants__;
    constants__.insert({"g_conv_aC", make<make<double>::sp_mat>::map_str(g_conv_aC_)});
    constants__.insert({"g_conv_aF", make<make<double>::sp_mat>::map_str(g_conv_aF_)});
    constants__.insert({"g_diff_aC", make<make<double>::sp_mat>::map_str(g_diff_aC_)});
    constants__.insert({"g_diff_aC", make<make<double>::sp_mat>::map_str(g_diff_aC_)});
    constants__.insert({"gc", make<make<double>::sp_mat>::map_str(gc_)});
    this->constants = constants__;
};
void minfo::make_clust(make<finfo>::map_int& faces)
{
    make<make<int>::vec>::map_int& fs2s_ = this->fid["s2s"];
    make<sparse_input>::vec fcs2s_input;
    make<sparse_input>::vec ccs2s_input;
    for(int i = 0; i < sizeof(fs2s_) - 1; i++)
    {
        for(int j = i; j < sizeof(fs2s_); j++)
        {
            ccs2s_input.push_back(sparse_input(i, j, 0.0));
        };
    };
    int tail = sizeof(fs2s_) - 1;
    ccs2s_input.push_back(sparse_input(0, 0, 0.0));
    ccs2s_input.push_back(sparse_input(tail, tail, 0.0));
    for(std::pair<int, make<int>::vec> entry : fs2s_)
    {
        for(auto i = entry.second.begin(); i != entry.second.end(); i++)
        {
            fcs2s_input.push_back(sparse_input(entry.first, *i, 0.0));
        };
    };
    make<double>::sp_mat fc__; fc__.setFromTriplets(fcs2s_input.begin(), fcs2s_input.end());
    make<double>::sp_mat cc__; cc__.setFromTriplets(ccs2s_input.begin(), ccs2s_input.end());
    // constants g_s2s (view)
    make<double>::sp_mat view_(cc__);
    coor norm1__; coor norm2__; coor v1__; coor v2__; double cos1__; double cos2__; double sq_sum__;
    for(int i = 0; i < view_.outerSize(); i++)
    {
        for(make<double>::sp_mat::InnerIterator it(view_, i); it; it++)
        {
            double view__ = 0.0;
            if(it.row() != it.col())
            {
                double area_tot = 0.0; 
                for(auto k = fs2s_[it.row()].begin(); k != fs2s_[it.row()].end(); k++)
                {
                    area_tot += faces[*k].farea;
                    for(auto l = fs2s_[it.col()].begin(); l != fs2s_[it.col()].end(); l++)
                    {
                        norm1__ = faces[*k].fnormal; norm2__ = faces[*l].fnormal;
                        v1__ = faces[*l].fcentroid - faces[*k].fcentroid;
                        v2__ = faces[*k].fcentroid - faces[*l].fcentroid;
                        if(norm1__.dot(v1__) < 0)
                        {
                            norm1__ = norm1__ * (-1);
                        };
                        if(norm2__.dot(v2__) < 0)
                        {
                            norm2__ = norm2__ * (-1);  
                        };
                        cos1__ = norm1__.dot(v1__) / sqrt_sum(v1__);
                        cos2__ = norm2__.dot(v2__) / sqrt_sum(v2__);
                        sq_sum__ = std::accumulate(v1__.begin(), v1__.end(), 0.0, square_sum);
                        view__ += (cos1__ * cos2__) / (*pi * sq_sum__);
                    };
                };
                view_.coeffRef(it.row(), it.col()) = view__ / area_tot;
            };
        };
    };
    for(int i = 0; i < sizeof(fs2s_); i++)
    {
        cc__.coeffRef(i, i) = 1.0;
    };
    this->fc.insert({"s2s", make<double>::sp_mat(fc__)});
    this->cc.insert({"s2s", make<double>::sp_mat(cc__)});
    this->constants.insert({"s2s", make<make<double>::sp_mat>::map_str()});
    this->constants["s2s"].insert({"view", make<double>::sp_mat(view_)});
};
void minfo::make_minfo(mshio::MshSpec spec)
{
    make<coor>::map_int nodes;
    make<finfo>::map_int faces;
    make<cinfo>::map_int cells;
    for(int i = 0; i < spec.nodes.num_entity_blocks; i++)
    {
        mshio::NodeBlock& block = spec.nodes.entity_blocks[i];
        for(int j = 0; j < block.num_nodes_in_block; j++)
        {
            coor node__(block.data[j*3+0], block.data[j*3+1], block.data[j*3+2]);
            nodes.insert({block.tags[j] - 1, node__});
        };
    };
    make<int>::map_int type_to_nface{{5,6}, {6,5}};
    int tag; int type; int n;
    for(int i = 0; i < spec.elements.num_entity_blocks; i++)
    {
        mshio::ElementBlock& block = spec.elements.entity_blocks[i];
        tag = block.entity_tag;
        type = block.element_type;
        n = mshio::nodes_per_element(type);
        switch(block.entity_dim)
        {
            case 2:
            for(int j = 0; j < block.num_elements_in_block; j++)
            {
                finfo face__;
                make<int>::vec fnode__{int(block.data[j * (n+1) - 1]), int(block.data[j * (n+1)]), int(block.data[j * (n+1) + 1])};
            };
            case 3:
            for(int j = 0; j < block.num_elements_in_block; j++)
            {
                cinfo cell__;
                make<int>::vec cnode__; for(int k = 0; k < n; k++) {cnode__.push_back(block.data[j * (n+1) + k - 1]);};
                make<int>::vec cface__; 
                for(std::pair<int, finfo> entry : faces)
                {
                    int ctd = 0;
                    for(auto k = entry.second.fnode.begin(); k != entry.second.fnode.end(); k++)
                    {
                        if(std::find(cnode__.begin(), cnode__.end(), *k) != cnode__.end())
                        {
                            ctd += 1;
                        };
                    };
                    if(ctd == sizeof(entry.second.fnode))
                    {
                        cface__.push_back(entry.first);
                    }
                    else if(sizeof(cface__) == type_to_nface[type])
                    {
                        break;
                    };
                };
                cell__.cnode = cnode__; cell__.cface = cface__; 
                cells.insert({block.data[j * (n+1)] - 1, cell__});
            };
            default:
            return;
        };
    };
    this->make_id(spec);
    this->make_template(nodes, faces, cells);
    this->make_size(nodes, faces, cells);
    this->make_geom(faces, cells);
    this->make_constants();
    make<make<make<int>::vec>::map_int>::map_str::iterator it_s2s = this->fid.find("s2s");
    if(it_s2s != this->fid.end())
    {
       this->make_clust(faces);
    };
};
// pinfo
void pinfo::make_pinfo(minfo& mesh_)
{
    make<make<make<int>::vec>::map_str>::map_str& cid_ = mesh_.cid;
    make<make<double>::sp_mat>::map_str& cc_ = mesh_.cc;
    for(std::pair<std::string, make<make<int>::vec>::map_str> entry_cid : cid_)
    {
        if(entry_cid.first.compare("fluid") == 0)
        {
            make<make<double>::map_int>::map_str rho_;
            make<make<double>::map_int>::map_str miu_;
            make<make<double>::map_int>::map_str cp_;
            for(std::pair<std::string, make<int>::vec> entry_fluid : entry_cid.second)
            {
                rho_.insert({entry_fluid.first, make<double>::map_int()});
                miu_.insert({entry_fluid.first, make<double>::map_int()});
                cp_.insert({entry_fluid.first, make<double>::map_int()});
                for(auto i = entry_fluid.second.begin(); i != entry_fluid.second.end(); i++)
                {
                    rho_[entry_fluid.first].insert({*i, 0.0});
                    miu_[entry_fluid.first].insert({*i, 0.0});
                    cp_[entry_fluid.first].insert({*i, 0.0});
                };
            };
            this->rho = rho_;
            this->miu = miu_;
            this->cp = cp_;
        }
        else if(entry_cid.first.compare("solid") == 0)
        {
            make<make<double>::map_int>::map_str k_;
            for(std::pair<std::string, make<int>::vec> entry_solid : entry_cid.second)
            {
                k_.insert({entry_solid.first, make<double>::map_int()});
                for(auto i = entry_solid.second.begin(); i != entry_solid.second.end(); i++)
                {
                    k_[entry_solid.first].insert({*i, 0.0});
                };
            };
            this->k = k_;
        };
    };
    make<make<double>::sp_mat>::map_str::iterator it_s2s = cc_.find("s2s");
    if(it_s2s != cc_.end())
    {
        make<double>::map_int eps_;
        make<double>::map_int alpha_;
        for(int i = 0; i < cc_["s2s"].rows(); i++)
        {
            eps_.insert({i, 0.0});
            alpha_.insert({i, 0.0});
        };
        this->eps = eps_;
        this->alpha = alpha_;
    };
};
// winfo
void winfo::make_winfo(minfo& mesh_)
{
    make<make<make<int>::vec>::map_str>::map_str& cid_ = mesh_.cid;
    make<make<make<int>::vec>::map_str>::map_str::iterator it_misc = cid_.find("misc");
    if(it_misc != cid_.end())
    {
        make<make<int>::vec>::map_str::iterator it_wall = cid_["misc"].find("wall");
        if(it_wall != cid_["misc"].end())
        {
            make<double>::map_int ts_;
            make<double>::map_int utau_;
            make<double>::map_int miut_;
            for(auto i = cid_["misc"]["wall"].begin(); i != cid_["misc"]["wall"].end(); i++)
            {
                ts_.insert({*i, 0.0}); utau_.insert({*i, 0.0}); miut_.insert({*i, 0.0});
            };
            this->ts = ts_;
            this->utau = utau_;
            this->miut = miut_;
        };
    };
};
// vinfo
void vinfo::make_vinfo(make<std::string>::vec which, minfo& mesh_)
{
    make<make<double>::map_int>::map_str value_;
    make<make<double>::map_int>::map_str prev_value_;
    make<make<coor>::map_int>::map_str grad_;
    make<make<coor>::map_int>::map_str prev_grad_;
    for(auto i = which.begin(); i != which.end(); i++)
    {
        make<double>::sp_mat& cc_ = mesh_.cc[*i];
        value_.insert({*i, make<double>::map_int()});
        grad_.insert({*i, make<coor>::map_int()});
        for(int j = 0; j < cc_.rows(); j++)
        {
            value_[*i].insert({j, 0.0});
            grad_[*i].insert({j, coor(0.0, 0.0, 0.0)});
        };
    };
    this->value = value_;
    this->prev_value = value_;
    this->grad = grad_;
    this->prev_grad = grad_;
};
};
#endif