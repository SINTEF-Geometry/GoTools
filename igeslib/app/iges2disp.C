#include "GoTools/igeslib/IGESconverter.h"
#include <iostream>

int main()
{
    IGESconverter conv;
    conv.readIGES(std::cin);
    conv.writedisp(std::cout);
}



