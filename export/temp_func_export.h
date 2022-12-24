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
void export::export_to_sql(cfdscheme::scheme& scheme_ref)
{
    sqlite3* database_p;
    int exit = 0;
    char* message_err;
    // open
    exit = sqlite3_open(this->output_name.c_str(), &database_p);
    // sql table templates
    // time step-wise values stored as string with comma limiter
    std::string commons_template =  "CREATE TABLE COMMONS("
                                    "NAME       TEXT            NOT NULL,"
                                    "N_CELLS    INT            NOT NULL);";
    exit = sqlite3_exec(database_p, commons_template.c_str(), NULL, 0, &message_err);
    std::string node_template = "CREATE TABLE NODE("
                                "ID         INT     PRIMARY KEY NOT NULL,"
                                "X          REAL               NOT NULL,"
                                "Y          REAL               NOT NULL,"
                                "Z          REAL               NOT NULL);";
    exit = sqlite3_exec(database_p, node_template.c_str(), NULL, 0, &message_err);
    std::string face_template = "CREATE TABLE FACE("
                                "ID         INT     PRIMARY KEY NOT NULL,"
                                "BOUNDARY   TEXT                NOT NULL,"
                                "NODES      TEXT                NOT NULL,"
                                "AREA       REAL                NOT NULL,"
                                "P          TEXT                NOT NULL,"
                                "U          TEXT                NOT NULL,"
                                "V          TEXT                NOT NULL,"
                                "W          TEXT                NOT NULL,"
                                "K          TEXT                NOT NULL,"
                                "E          TEXT                NOT NULL,"
                                "T          TEXT                NOT NULL);";
    exit = sqlite3_exec(database_p, face_template.c_str(), NULL, 0, &message_err);
    std::string cell_template = "CREATE TABLE CELL("
                                "ID         INT     PRIMARY KEY NOT NULL,"
                                "DOMAIN     TEXT                NOT NULL,"
                                "NODES      TEXT                NOT NULL,"
                                "VOLUME     REAL                NOT NULL,"
                                "P          TEXT                NOT NULL,"
                                "U          TEXT                NOT NULL,"
                                "V          TEXT                NOT NULL,"
                                "W          TEXT                NOT NULL,"
                                "K          TEXT                NOT NULL,"
                                "E          TEXT                NOT NULL,"
                                "T          TEXT                NOT NULL);";
    exit = sqlite3_exec(database_p, cell_template.c_str(), NULL, 0, &message_err);
    std::string err_template =  "CREATE TABLE ERROR("
                                "U          TEXT                NOT NULL,"
                                "V          TEXT                NOT NULL,"
                                "W          TEXT                NOT NULL,"
                                "K          TEXT                NOT NULL,"
                                "E          TEXT                NOT NULL,"
                                "T          TEXT                NOT NULL);";
    exit = sqlite3_exec(database_p, err_template.c_str(), NULL, 0, &message_err);
    std::string res_template =  "CREATE TABLE RESIDUAL("
                                "U          TEXT                NOT NULL,"
                                "V          TEXT                NOT NULL,"
                                "W          TEXT                NOT NULL,"
                                "K          TEXT                NOT NULL,"
                                "E          TEXT                NOT NULL,"
                                "T          TEXT                NOT NULL);";
    exit = sqlite3_exec(database_p, res_template.c_str(), NULL, 0, &message_err);
    std::string time_iter_template = "CREATE TABLE TIME_ITER("
                                     "TIME       TEXT                NOT NULL,"
                                     "ITER       TEXT                NOT NULL);";
    exit = sqlite3_exec(database_p, time_iter_template.c_str(), NULL, 0, &message_err);
    // commons
    std::string commons_input("INSERT INTO COMMONS VALUES(");
    commons_input += "'" + this->mesh_name + "'" + ", ";
    commons_input += std::to_string(this->number_of_cells) + ");";
    exit = sqlite3_exec(database_p, commons_input.c_str(), NULL, 0, &message_err);
    // nodes
    for(std::pair<int, coor> entry : scheme_ref.mesh.nodes)
    {
        std::string nodes_input("INSERT INTO NODES VALUES(");
        nodes_input += std::to_string(entry.first) + ", ";
        nodes_input += std::to_string(entry.second(0)) + ", ";
        nodes_input += std::to_string(entry.second(1)) + ", ";
        nodes_input += std::to_string(entry.second(2)) + ");";
        exit = sqlite3_exec(database_p, nodes_input.c_str(), NULL, 0, &message_err);
    };
    // face
    for(std::pair<std::string, make<make<int>::vec>::map_int> entry_str : scheme_ref.mesh.fid)
    {
        for(std::pair<int, make<int>::vec> entry_int : entry_str.second)
        {
            for(auto i = entry_int.second.begin(); i != entry_int.second.end(); i++)
            {
                std::string faces_input("INSERT INTO FACES VALUES(");
                faces_input += std::to_string(*i) + ", ";
                faces_input += "'" + entry_str.first + "_" + std::to_string(entry_int.first) + "'" + ", ";
                faces_input += vec_to_str<int>(scheme_ref.mesh.faces[*i].fnode) + ", ";
                faces_input += std::to_string(scheme_ref.mesh.faces[*i].farea) + ", ";
                faces_input += vec_to_str<double>(this->face_export[*i]["P"]) + ", ";
                faces_input += vec_to_str<double>(this->face_export[*i]["u"]) + ", ";
                faces_input += vec_to_str<double>(this->face_export[*i]["v"]) + ", ";
                faces_input += vec_to_str<double>(this->face_export[*i]["w"]) + ", ";
                faces_input += vec_to_str<double>(this->face_export[*i]["k"]) + ", ";
                faces_input += vec_to_str<double>(this->face_export[*i]["e"]) + ", ";
                faces_input += vec_to_str<double>(this->face_export[*i]["T"]) + ");";
                exit = sqlite3_exec(database_p, faces_input.c_str(), NULL, 0, &message_err);
            };
        };
    };
    // cell
    for(std::pair<std::string, make<make<int>::vec>::map_str> entry_str1 : scheme_ref.mesh.cid)
    {
        for(std::pair<std::string, make<int>::vec> entry_str2 : entry_str1.second)
        {
            for(auto i = entry_str2.second.begin(); i != entry_str2.second.end(); i++)
            {
                std::string cells_input("INSERT INTO CELLS VALUES(");
                cells_input += std::to_string(*i) + ", ";
                cells_input += "'" + entry_str2.first + "'" + ", ";
                cells_input += vec_to_str<int>(scheme_ref.mesh.cells[*i].cnode) + ", ";
                cells_input += std::to_string(scheme_ref.mesh.cells[*i].cvolume) + ", ";
                cells_input += vec_to_str<double>(this->cell_export[*i]["P"]) + ", ";
                cells_input += vec_to_str<double>(this->cell_export[*i]["u"]) + ", ";
                cells_input += vec_to_str<double>(this->cell_export[*i]["v"]) + ", ";
                cells_input += vec_to_str<double>(this->cell_export[*i]["w"]) + ", ";
                cells_input += vec_to_str<double>(this->cell_export[*i]["k"]) + ", ";
                cells_input += vec_to_str<double>(this->cell_export[*i]["e"]) + ", ";
                cells_input += vec_to_str<double>(this->cell_export[*i]["T"]) + ");";
                exit = sqlite3_exec(database_p, cells_input.c_str(), NULL, 0, &message_err);
            };
        };
    };
    // err
    std::string err_input("INSERT INTO ERROR VALUES(");
    err_input += vec_to_str<double>(this->err_export["u"]) + ", ";
    err_input += vec_to_str<double>(this->err_export["v"]) + ", ";
    err_input += vec_to_str<double>(this->err_export["w"]) + ", ";
    err_input += vec_to_str<double>(this->err_export["k"]) + ", ";
    err_input += vec_to_str<double>(this->err_export["e"]) + ", ";
    err_input += vec_to_str<double>(this->err_export["T"]) + ");";
    exit = sqlite3_exec(database_p, err_input.c_str(), NULL, 0, &message_err);
    // res
    std::string res_input("INSERT INTO RESIDUAL VALUES(");
    res_input += vec_to_str<double>(this->res_export["u"]) + ", ";
    res_input += vec_to_str<double>(this->res_export["v"]) + ", ";
    res_input += vec_to_str<double>(this->res_export["w"]) + ", ";
    res_input += vec_to_str<double>(this->res_export["k"]) + ", ";
    res_input += vec_to_str<double>(this->res_export["e"]) + ", ";
    res_input += vec_to_str<double>(this->res_export["T"]) + ");";
    exit = sqlite3_exec(database_p, res_input.c_str(), NULL, 0, &message_err);
    // time_iter
    std::string time_iter_input("INSERT INTO TIME_ITER VALUES(");
    time_iter_input += vec_to_str<long long>(this->time_export) + ", ";
    time_iter_input += vec_to_str<int>(this->iter_export) + ");";
    exit = sqlite3_exec(database_p, time_iter_input.c_str(), NULL, 0, &message_err);
    // close
    sqlite3_close(database_p);
};

#endif