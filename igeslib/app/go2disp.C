#include "GoTools/igeslib/IGESconverter.h"
#include <iostream>

using namespace Go;

using std::cin;
using std::cout;


int main()
{
    IGESconverter conv;
    conv.readgo(cin);
    conv.writedisp(cout);
}



