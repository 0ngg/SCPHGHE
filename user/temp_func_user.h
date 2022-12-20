#include"temp_struct_user.h"

#ifdef CFDUSER_H

// user func
void user::read_source_csv(std::string source_file)
{
    make<make<std::string>::vec>::vec parsed = parse_csv(source_file);
    // make<make<make<double>::map_int>::map_int>::map_str face_source; // .... - time step - value
    // make<make<make<double>::map_int>::map_str>::map_str cell_source; // .... - time step - value
    for(auto i = parsed.begin(); i != parsed.end(); i++)
    {
        make<std::string>::vec temp = *i;
        char test_ifdigit = temp[1][0];
        if(std::isdigit(test_ifdigit))
        {
            // face_source
            make<make<make<double>::map_int>::map_int>::map_str::iterator ifface = this->face_source.find(temp[0]);
            if(ifface == this->face_source.end())
            {
                this->face_source.insert({temp[0], make<make<double>::map_int>::map_int()});
                make<double>::map_int value_vec;
                int ctd = 0;
                for(int j = 2; j < temp.size(); j++)
                {
                    value_vec.insert({ctd, std::stod(temp[j])});
                    ctd += 1;
                };
                this->face_source[temp[0]].insert({std::stoi(temp[1]), make<double>::map_int(value_vec)});
            }  
            else
            {
                make<double>::map_int value_vec;
                int ctd = 0;
                for(int j = 2; j < temp.size(); j++)
                {
                    value_vec.insert({ctd, std::stod(temp[j])});
                    ctd += 1;
                };
                this->face_source[temp[0]].insert({std::stoi(temp[1]), make<double>::map_int(value_vec)});
            };
        }
        else
        {
            // cell_source
            make<make<make<double>::map_int>::map_str>::map_str::iterator ifcell = this->cell_source.find(temp[0]);
            if(ifcell == this->cell_source.end())
            {
                this->cell_source.insert({temp[0], make<make<double>::map_int>::map_str()});
                make<double>::map_int value_vec;
                int ctd = 0;
                for(int j = 2; j < temp.size(); j++)
                {
                    value_vec.insert({ctd, std::stod(temp[j])});
                    ctd += 1;
                };
                this->cell_source[temp[0]].insert({temp[1], make<double>::map_int(value_vec)});
            }
            else
            {
                make<double>::map_int value_vec;
                int ctd = 0;
                for(int j = 2; j < temp.size(); j++)
                {
                    value_vec.insert({ctd, std::stod(temp[j])});
                    ctd += 1;
                };
                this->cell_source[temp[0]].insert({temp[1], make<double>::map_int(value_vec)});  
            };
        };
    };
};
void user::read_solid_prop_csv(std::string prop_file)
{
    make<make<std::string>::vec>::vec parsed = parse_csv(prop_file);
    for(auto i = parsed.begin(); i != parsed.end(); i++)
    {
        make<std::string>::vec temp = *i;
        char test_ifdigit = temp[1][0];
        if(std::isdigit(test_ifdigit))
        {
            // s2s_eps
            this->s2s_eps.insert({std::stoi(temp[1]), std::stod(temp[2])});
        }
        else
        {
            // solid_k
            this->solid_k.insert({temp[1], std::stod(temp[2])});
        };
    };
};
void user::update_source(int current_time, cfdscheme::scheme& scheme_ref)
{
    for(std::pair<std::string, make<double>::map_int> entry_face : scheme_ref.source.face_value)
    {
        for(std::pair<int, double> entry_face_int : entry_face.second)
        {
            entry_face_int.second = this->face_source[entry_face.first][entry_face_int.first][current_time];
        };
    };
    for(std::pair<std::string, make<double>::map_str> entry_cell : scheme_ref.source.cell_value)
    {
        for(std::pair<std::string, double> entry_cell_str : entry_cell.second)
        {
            entry_cell_str.second = this->cell_source[entry_cell.first][entry_cell_str.first][current_time];
        };
    };
};
#endif
