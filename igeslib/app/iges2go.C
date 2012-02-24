#include "GoTools/igeslib/IGESconverter.h"
#include <iostream>
#include <iomanip>

using namespace Go;

int main()
{
    IGESconverter conv;
    conv.readIGES(std::cin);
    std::cout << std::setprecision(15) << std::flush;
    conv.writego(std::cout);
}



