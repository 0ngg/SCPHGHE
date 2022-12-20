#include"temp_util_user.h"
#include"temp_struct_scheme.h"

#ifdef CFDUSER_H

class user
{
    public:
    double P_init;
    double T_init;
    double W_init;
    make<double>::map_str solid_k;
    make<double>::map_int s2s_eps;
    make<make<make<double>::map_int>::map_int>::map_str face_source; // .... - time step - value
    make<make<make<double>::map_int>::map_str>::map_str cell_source; // .... - time step - value
    user() {};
    user(double P_init_in, double T_init_in, double W_init_in, std::string source_file,
         std::string prop_solid_file): P_init(P_init_in), T_init(T_init_in), W_init(W_init_in)
    {
        this->read_source_csv(source_file);
        this->read_solid_prop_csv(prop_solid_file);
    };
    void update_source(int, cfdscheme::scheme&);
    private:
    void read_source_csv(std::string);
    void read_solid_prop_csv(std::string);
};

#endif