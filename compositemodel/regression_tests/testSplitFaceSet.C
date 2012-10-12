#define BOOST_TEST_MODULE compositemodel/splitFaceSet
#include <boost/test/unit_test.hpp>

#include <fstream>
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/RegularizeFaceSet.h"

using namespace std;
using namespace Go;


struct Config {
public:
    Config()
    {
        datadir = "../data/step_reader/data2/"; // Relative to build/compositemodel

        infiles.push_back("b25.g2");
        numfaces.push_back(20);

        infiles.push_back("b26.g2");
        numfaces.push_back(16);

        gap = 0.001; // 0.001;
        neighbour = 0.01; // 0.01;
        kink = 0.01;
        approxtol = 0.01;

    }

public:
    string datadir;
    vector<string> infiles;
    vector<int> numfaces;
    double gap;
    double neighbour;
    double kink;
    double approxtol;
};


BOOST_FIXTURE_TEST_CASE(splitFaceSet, Config)
{
    int nfiles = infiles.size();
    for (int i = 0; i < nfiles; ++i) {

        string infile = datadir + infiles[i];

        // Read input arguments
        std::ifstream file1(infile.c_str());
        BOOST_CHECK_MESSAGE(file1.good(), "Input file not found or file corrupt");

        std::ofstream file2("outfile.g2");

        CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

        CompositeModel *model = factory.createFromG2(file1);

        SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);

        if (sfmodel)
        {
            std::vector<shared_ptr<ftSurface> > faces = sfmodel->allFaces();

            //RegularizeFaceSet reg(faces, gap, kink, true);
            RegularizeFaceSet reg(faces, gap, kink, false);
            std::vector<shared_ptr<ftSurface> > sub_faces = reg.getRegularFaces();

            int nsubfaces = sub_faces.size();
            BOOST_CHECK_EQUAL(nsubfaces, numfaces[i]);

            for (size_t ki=0; ki < nsubfaces; ++ki)
            {
                shared_ptr<ParamSurface> surf = sub_faces[ki]->surface();
                surf->writeStandardHeader(file2);
                surf->write(file2);
            }
        }
    }
}

