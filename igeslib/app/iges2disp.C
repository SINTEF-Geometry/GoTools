#include "GoTools/igeslib/IGESconverter.h"
#include <iostream>

using namespace Go;

int main()
{
    IGESconverter conv;
    conv.readIGES(std::cin);
    conv.writedisp(std::cout);
}



