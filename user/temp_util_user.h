#ifndef CFDUSER_H
#define CFDUSER_H

#include"temp_util.h"

make<make<std::string>::vec>::vec parse_csv(std::string filename)
{
    // output strings, revert to double later
    std::ifstream data(filename);
    std::string line;
    make<make<std::string>::vec>::vec parsed;
    while(std::getline(data, line))
    {
        std::stringstream lineStream(line);
        std::string element;
        make<std::string>::vec parsed_row;
        while(std::getline(lineStream, element, ','))
        {
            parsed_row.push_back(element);
        };
        parsed.push_back(parsed_row);
    };
    return parsed;
};

#endif