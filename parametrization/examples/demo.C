/*----------------------------------------------*/
/* Test program for creating a PR type triangulation from
   given node and triangle data, parametrizing it

   Michael Floater, September 1997.                          
   Revised, Michael Floater, Nov. 1998 
   Revised to suit the Go library, Michael Floater, April. 2001.
   Added Mean Value Parameterization, April 2002.    */
/*----------------------------------------------*/

#include "GoTools/parametrization/PrFastUnorganized_OP.h" // krull
#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrParametrizeInt.h"
#include "GoTools/parametrization/PrParametrizeBdy.h"
#include "GoTools/parametrization/PrPrmShpPres.h"
#include "GoTools/parametrization/PrPrmUniform.h"
#include "GoTools/parametrization/PrPrmLeastSquare.h"
#include "GoTools/parametrization/PrPrmEDDHLS.h"
#include "GoTools/parametrization/PrPrmMeanValue.h"
#include <fstream>
#include <string>
#include <cstring>
#include <memory>
// #include <boost/tokenizer.hpp>

using std::cerr;
using std::endl;
using std::cout;
using std::string;
// using namespace boost;

string suffix_of(char* filename);

int main(int argc, const char** argv)
{
    if (argc < 2) {
	cerr << "This sample program generates a parametrization for a given triangulation, and " << endl;
	cerr << "writes the result to disk." << endl;
	cerr << "Usage: demo <-option argument>, where options/arguments are: " << endl;
	cerr << " -i <filename>      : specify file containing the input triangulation." << endl;
	cerr << "                      If the filename suffix is 'pcloud', it is assumed to " << endl;
	cerr << "                      be a pointcloud.  Otherwise it is assumed to be a " << endl;
	cerr << "                      triangulation." << endl;
	cerr << " -intpar <integer>  : specify method for parametrizing interior of triangulation" << endl;
	cerr << "                      Valid options are: 1 - shape-preserving parametrization" << endl;
	cerr << "                                         2 - uniform parametrization" << endl;
	cerr << "                                         3 - least-square based parametrization" << endl;
	cerr << "                                         4 - the EDDHLS parametrization" << endl;
	cerr << "                                         5 - mean value parametrization" << endl;
	cerr << " -bdypar <integer>  : specify method for parametrizing boundary of triangulation " << endl;
	cerr << "                      Valid options are: 1 - chord length parametrization " << endl;
	cerr << "                                         2 - centripetal parametrization " << endl;
	cerr << "                                         3 - uniform boundary parametrization " << endl;
	cerr << " -bdymeth <integer> : specify shape of parameter domain " << endl;
	cerr << "                      Valid options are: 1 - circle" << endl;
	cerr << "                                         2 - rectangle, based on shape of 3D " << endl;
	cerr << "                                             geometry" << endl;
	cerr << "                                         3 - rectangle, based on extent of current" << endl;
	cerr << "                                             parametrization.  Don't use this if " << endl;
	cerr << "                                             your data is not already " << endl;
	cerr << "                                             parametrized somehow. " << endl;
	cerr << "The result (the parametrized triangulation) will be written to the file" << endl;
	cerr << "\"triangulation\"." << endl;
	return 0;
    }

    bool inf_specified = false;
    char inf[80];
    int intparam_type = 5; // default is mean value parametrization
    int bdyparam_type = 1;
    int bdy_method = 1;

    int undefCmd = 0;

    for (int i=1; i<argc; i=i+2)
	{
	    if(strcmp(argv[i],"-i") == 0) {
		strncpy(inf,argv[i+1],80);
		inf_specified = true;
	    } else if( strcmp(argv[i],"-intpar") == 0 ) {
		intparam_type = atoi(argv[i+1]);
	    } else if( strcmp(argv[i],"-bdypar") == 0 ) {
		bdyparam_type = atoi(argv[i+1]);
	    } else if ( strcmp(argv[i],"-bdymeth") == 0) {
		bdy_method = atoi(argv[i+1]);
	    } else {
		undefCmd = 1;
	    }
	}

    if(undefCmd)
	{
	    cerr << "Error in command line" << endl;
	    return -1;
	}
    if (!inf_specified) {
	cerr << "No input file specified.  Aborting." << endl;
	return -1;
    }
    std::ifstream infile(inf);
    if (!infile) {
	cerr << "Unable to open file '" << inf <<  "'. Aborting." << endl;
	return -1;
    }

    bool is_pointcloud = (suffix_of(inf) == string("pcloud"));

    // Get the data from the infile
    int no_comps = 1, genus = 1;
    shared_ptr<PrOrganizedPoints> pr_triang;
    if (is_pointcloud) {
	shared_ptr<PrFastUnorganized_OP> tmp(new PrFastUnorganized_OP);
	tmp->scanRawData(infile);
	tmp->initNeighbours();
	no_comps = tmp->findNumComponents();
	pr_triang = tmp;
	
    } else {
	shared_ptr<PrTriangulation_OP> tmp(new PrTriangulation_OP);
	tmp->scanRawData(infile);
	tmp->printInfo(cout);
	no_comps = tmp->findNumComponents();
	genus = tmp->findGenus();
	pr_triang = tmp;

    }

    if(no_comps == 1 && genus == 1)
	{
	    PrParametrizeBdy pb;
  
	    switch(bdyparam_type)
		{
		case 1: pb.setParamKind(PrCHORDLENGTHBDY); break;
		case 2: pb.setParamKind(PrCENTRIPETAL); break;
		case 3: pb.setParamKind(PrUNIFBDY); break;
		}

	    pb.attach(pr_triang);

	    cout << "Parametrizing boundary..." << endl;
	    int cor[4];

	    if(bdy_method == 1)
		{
		    pb.parametrize(); // circle
		}
	    else
		{
		    if(bdy_method == 2) {
			pb.findCornersFromXYZ(cor);
		    } else if(bdy_method == 3) {
			pb.findCornersFromUV(cor);
		    } else {
			cerr << "erroneous boundary shape type.  Aborting. " << endl;
			return -1;
		    }

		    for(int k=0; k<4; k++)
			{
			    cout << "corner(i,u,v,x,y,z) = " << cor[k];
			    cout << " " << pr_triang->getU(cor[k]);
			    cout << " " << pr_triang->getV(cor[k]);
			    cout << " " << pr_triang->get3dNode(cor[k]).x();
			    cout << " " << pr_triang->get3dNode(cor[k]).y();
			    cout << " " << pr_triang->get3dNode(cor[k]).z();
			    cout << endl;
			}
		    pb.parametrize(cor[0],cor[1],cor[2],cor[3]); // square boundary
		}

	    cout << "Parametrizing interior..." << endl;

	    PrParametrizeInt *pi;
	    switch(intparam_type)
		{
		case 1: pi = new PrPrmShpPres; break;
		case 2: pi = new PrPrmUniform; break;
		case 3: pi = new PrPrmLeastSquare; break;
		case 4: pi = new PrPrmEDDHLS; break;
		case 5: pi = new PrPrmMeanValue; break;
		}
  
	    pi->attach(pr_triang);
	    pi->parametrize();
	    delete pi;

	    std::ofstream uv_nodes_file("uv_nodes");
	    pr_triang->printUVNodes(uv_nodes_file);

	    //     std::ofstream uv_edges_file("uv_edges");
	    //     pr_triang->printUVEdges(uv_edges_file);
	    //     std::ofstream uv_triangles_file("uv_triangles");

	    //     ofstream krull("uv_triangles_file");
	    //     pr_triang->printUVTriangles(krull);

	    //     Std::ofstream uvxyz_nodes_file("uvxyz_nodes");
	    //     pr_triang->printUVXYZNodes(uvxyz_nodes_file);

	    if (!is_pointcloud) {
		std::ofstream triangulation_file("triangulation");
		dynamic_pointer_cast<PrTriangulation_OP, PrOrganizedPoints>(pr_triang)->print(triangulation_file);
	    }
	}
    else return 0;
}



string suffix_of(char* filename)
{
    // TESTME! @jbt
    string s(filename);
    size_t pos = s.find_last_of(".");
    if (pos > s.size())
        return "";
    return s.substr(pos, s.size()-pos);

    // typedef boost::tokenizer<char_separator<char> > tokenizer;
    // string s(filename);
    // boost::char_separator<char> sep(".");
    // tokenizer tok(s, sep);
    // tokenizer::iterator last;
    // for (tokenizer::iterator it = tok.begin(); it != tok.end(); last = it++);
    // return *last;
}
