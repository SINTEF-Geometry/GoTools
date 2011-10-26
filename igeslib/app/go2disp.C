#include "GoTools/igeslib/IGESconverter.h"
#include <iostream>


using std::cin;
using std::cout;


int main()
{
    IGESconverter conv;
    conv.readgo(cin);
    conv.writedisp(cout);
}



