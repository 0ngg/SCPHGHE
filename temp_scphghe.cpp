#include"temp_util.h"
#include"temp_struct_user.h"
#include"temp_struct_scheme.h"
#include"temp_struct_solver.h"
#include"temp_struct_export.h"

int main(int argc, char** argv)
{
    // args
    // mesh name, output name, source file, prop_solid file
    // double P_init, double T_init, double W_init
    // double step_length, double under_relax, double min_residual, int max_iter;
    std::cout << "Initializing..." << std::endl;
    user main_user(std::stod(std::string(argv[4])), std::stod(std::string(argv[5])), std::stod(std::string(argv[6])), std::string(argv[2]), std::string(argv[3]));
    std::cout << "User obj. successfully built." << std::endl;
    cfdscheme::scheme main_scheme(std::string(argv[0]), main_user);
    std::cout << "Scheme obj. successfully built." << std::endl;
    exports main_export(std::string(argv[1]), std::string(argv[0]), main_scheme);
    std::cout << "Export obj. successfully built." << std::endl;
    cfdsolver::scphghe main_scphghe(main_scheme, main_user, std::stod(std::string(argv[7])), std::stod(std::string(argv[8])), std::stod(std::string(argv[9])),
                                    std::stoi(std::string(argv[10])));
    std::cout << "SCPHGHE obj. successfully built." << std::endl;
    main_scphghe.iterate(main_scheme, main_export, main_user);
    std::cout << "Done." << std::endl;
    return 0;
};