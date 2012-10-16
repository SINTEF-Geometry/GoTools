#define BOOST_TEST_MODULE gotools-core/testBoundedSurface
#include <boost/test/unit_test.hpp>

#include <fstream>
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedUtils.h"

using namespace std;
using namespace Go;


struct Config {
public:
    Config()
    {
        datadir = "data/"; // Relative to build/gotools-core

        infiles.push_back("test_bounded_sf_2.g2");
        numloops.push_back(1);
        num_cvs_in_loop.push_back(4);

        GoTools::init();
    }

public:
    string datadir;
    vector<string> infiles;
    vector<int> numloops;
    vector<int> num_cvs_in_loop;
};


BOOST_FIXTURE_TEST_CASE(testBoundedSurface, Config)
{
    int nfiles = infiles.size();
    for (int i = 0; i < nfiles; ++i) {

        string infile = datadir + infiles[i];

        // Read input arguments
        std::ifstream file1(infile.c_str());
        BOOST_CHECK_MESSAGE(file1.good(), "Input file not found or file corrupt");

        std::ofstream file2("outfile.g2");

        ObjectHeader header;
        header.read(file1);

        shared_ptr<BoundedSurface> bs(new BoundedSurface());
        bs->read(file1);

        vector<CurveLoop> loops = bs->allBoundaryLoops();
        int nloops = loops.size();
        BOOST_CHECK_EQUAL(nloops, numloops[i]);

        int valid_state = 0;
        bool is_valid = bs->isValid(valid_state);
        BOOST_MESSAGE("BoundedSurface valid state: " << valid_state);

        Go::BoundedUtils::fixInvalidBoundedSurface(bs);
        is_valid = bs->isValid(valid_state);
        BOOST_MESSAGE("BoundedSurface valid state after fixing: " 
            << valid_state);

    }
}

