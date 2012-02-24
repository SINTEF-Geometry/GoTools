#include "GoTools/igeslib/IGESconverter.h"
#include <iostream>
#include <iomanip>

using namespace Go;

int main()
{
    IGESconverter conv;
    conv.readdisp(std::cin);
    conv.writego(std::cout);
}



