#include "GoTools/creators/CurveCreators.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/creators/SmoothCurveSet.h"
#include <fstream>


using namespace Go;
using std::shared_ptr;
using std::vector;


int main(int argc, char* argv[])
{
    if (argc != 6) {
	MESSAGE("Usage: inputfile nmb_crvs datafile approxweight outputfile.");
	return 0;
    }

    // Read input arguments
    std::ifstream infile(argv[1]);
    ALWAYS_ERROR_IF(infile.bad(), "Input file not found or file corrupt");

    int nmb_crv = atoi(argv[2]);

    std::ifstream datafile(argv[3]);

    double wgt = atof(argv[4]);

    std::ofstream outfile(argv[5]);

    vector<shared_ptr<SplineCurve> > crvs;
    int ki, kj;
    for (ki=0; ki<nmb_crv; ki++)
      {
	ObjectHeader header;
	header.read(infile);
	shared_ptr<SplineCurve> curr = 
	  shared_ptr<SplineCurve>(new SplineCurve());
	curr->read(infile);
	crvs.push_back(curr);
      }

    // Keep the endpoints of all curves fixed
    vector<vector<int> > coef_known;
    coef_known.resize(nmb_crv);
    for (ki=0; ki<nmb_crv; ki++)
      {
	int kn = crvs[ki]->numCoefs();
	coef_known[ki].resize(kn);
	std::fill(coef_known[ki].begin(), coef_known[ki].end(), 0);
	coef_known[ki][0] = coef_known[ki][kn-1] = 1;
      }

    // No closed seem
    vector<int> seem(nmb_crv, 0);

    // Read data file
    int nmb_approx;
    vector<vector<double> > approx_pos;
    vector<vector<double> > approx_par;
    vector<vector<double> > approx_wgt;
    approx_pos.resize(nmb_crv);
    approx_par.resize(nmb_crv);
    approx_wgt.resize(nmb_crv);
    datafile >> nmb_approx;
    for (ki=0; ki<nmb_approx; ++ki)
      {
	int cv_idx;
	datafile >> cv_idx;
	double tmp;
	for (kj=0; kj<3; ++kj)
	  {
	    datafile >> tmp;
	    approx_pos[cv_idx].push_back(tmp);
	  }
	datafile >> tmp;
	approx_par[cv_idx].push_back(tmp);
	datafile >> tmp;
	approx_wgt[cv_idx].push_back(tmp);
      }
	
    int nmb_inter;
    vector<vector<double> > inter_pos;
    vector<vector<double> > inter_par;
    vector<vector<int> > inter_der;
    inter_pos.resize(nmb_crv);
    inter_par.resize(nmb_crv);
    inter_der.resize(nmb_crv);
    datafile >> nmb_inter;
    for (ki=0; ki<nmb_inter; ++ki)
      {
	int cv_idx;
	datafile >> cv_idx;
	double tmp;
	for (kj=0; kj<3; ++kj)
	  {
	    datafile >> tmp;
	    inter_pos[cv_idx].push_back(tmp);
	  }
	datafile >> tmp;
	inter_par[cv_idx].push_back(tmp);
	int der;
	datafile >> der;
	inter_der[cv_idx].push_back(der);
      }
	
    int nmb_side;
    vector<shared_ptr<cvSetConstraint> > constraints;
    datafile >> nmb_side;
    for (ki=0; ki<nmb_side; ++ki)
      {
	int cv_idx1, cv_idx2;
	double par1, par2;
	datafile >> cv_idx1;
	datafile >> cv_idx2;
	datafile >> par1;
	datafile >> par2;
	shared_ptr<cvSetConstraint> curr = 
	  shared_ptr<cvSetConstraint>(new cvSetConstraint(cv_idx1, par1, 0,
							  cv_idx2, par2, 0,
							  false));
	constraints.push_back(curr);
      }

    SmoothCurveSet smooth;
    smooth.attach(crvs, seem, coef_known, nmb_side);

    double wgt2, wgt3;
    wgt2 = 0.5*(1.0 - wgt);
    wgt3 = 0.5*(1.0 - wgt);
    smooth.setOptimize(0.0, wgt2, wgt3);

    if (nmb_approx > 0)
      smooth.setLeastSquares(approx_pos, approx_par, approx_wgt, wgt);

    int kstat = 0;
    if (nmb_inter > 0)
      smooth.setInterpolationConditions(inter_pos, inter_par, inter_der,
					false, wgt, &kstat);

    if (nmb_side > 0)
	smooth.setCvSetConstraints(constraints, false, wgt);

    vector<shared_ptr<SplineCurve> > out_crvs;
    smooth.equationSolve(out_crvs);

    for (ki=0; ki<(int)out_crvs.size(); ++ki)
      {
	out_crvs[ki]->writeStandardHeader(outfile);
	out_crvs[ki]->write(outfile);
      }
}

    

