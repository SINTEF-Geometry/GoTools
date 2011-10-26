//==========================================================================
//                                                                          
// File: test_errormacros.C                                                
//                                                                          
// Created:     
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision:
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/utils/errormacros.h"
#include <iostream>


using namespace std;
// using namespace Go;


int main()
{
    cout << "This is a MESSAGE:" << endl;
    MESSAGE("MESSAGE text");
    cout << "That was a MESSAGE." << endl;


    return 0;
}
