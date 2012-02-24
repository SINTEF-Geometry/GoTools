#include "GoTools/igeslib/IGESconverter.h"
#include <iostream>

using namespace Go;
using namespace std;

int main()
{
    IGESconverter conv;
    conv.readgo(cin);
    conv.writeIGES(cout);
}



