#include "GoTools/igeslib/IGESconverter.h"
#include <iostream>
#include <iomanip>


int main()
{
    IGESconverter conv;
    conv.readdisp(std::cin);
    conv.writego(std::cout);
}



