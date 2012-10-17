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

        infiles.push_back("test_bounded_sf_3.g2");

        GoTools::init();
    }

public:
    string datadir;
    vector<string> infiles;
    vector<int> numloops;
    vector<int> num_cvs_in_loop;
};


BOOST_AUTO_TEST_CASE(testBoundedSurface)
{
    GoTools::init();
    ObjectHeader header;
    string datadir = "data/";


    // *** test_bounded_sf_2 ***
    {
        string filename = "test_bounded_sf_2.g2";
        string infile = datadir + filename;
        ifstream in(infile.c_str());
        BOOST_CHECK_MESSAGE(in.good(), "Input file not found or file corrupt");
        header.read(in);
        shared_ptr<BoundedSurface> bs(new BoundedSurface());
        bs->read(in);

        vector<CurveLoop> loops = bs->allBoundaryLoops();
        int nloops = loops.size();
        BOOST_CHECK_EQUAL(nloops, 1);

        int valid_state = 0;
        bool is_valid = bs->isValid(valid_state);
        BOOST_CHECK_MESSAGE(is_valid, "BoundedSurface valid state: " << valid_state);

        Go::BoundedUtils::fixInvalidBoundedSurface(bs);
        is_valid = bs->isValid(valid_state);
        BOOST_CHECK_MESSAGE(is_valid, "BoundedSurface valid state after fixing: " 
                            << valid_state);
    }

    // *** test_bounded_sf_3 ***
    {
        string filename = "test_bounded_sf_3.g2";
        string infile = datadir + filename;
        ifstream in(infile.c_str());
        BOOST_CHECK_MESSAGE(in.good(), "Input file not found or file corrupt");
        header.read(in);
        shared_ptr<BoundedSurface> bs(new BoundedSurface());
        bs->read(in);
        vector<CurveLoop> loops = bs->allBoundaryLoops();
        int nloops = loops.size();
        BOOST_CHECK_EQUAL(nloops, 1);

        int valid_state = 0;
        bool is_valid = bs->isValid(valid_state);
        BOOST_CHECK_MESSAGE(is_valid, "BoundedSurface valid state: " << valid_state);

        Go::BoundedUtils::fixInvalidBoundedSurface(bs);
        is_valid = bs->isValid(valid_state);
        BOOST_CHECK_MESSAGE(is_valid, "BoundedSurface valid state after fixing: " 
                            << valid_state);
    }



}

