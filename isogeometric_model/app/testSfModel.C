#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/isogeometric_model/IsogeometricSfModel.h"
#include "GoTools/isogeometric_model/IsogeometricSfBlock.h"
#include "GoTools/isogeometric_model/BdConditionType.h"
#include "GoTools/isogeometric_model/SfSolution.h"
#include "GoTools/isogeometric_model/SfBoundaryCondition.h"
#include "GoTools/isogeometric_model/BdCondFunctor.h"
#include <fstream>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::shared_ptr;
using std::dynamic_pointer_cast;
using std::min;

using namespace Go;

class DirichletFunctor : public BdCondFunctor
{
public:
  DirichletFunctor(std::shared_ptr<SplineCurve> geom_crv,
		   std::shared_ptr<SplineCurve> cond_crv);
  virtual ~DirichletFunctor();
  virtual Point evaluate(const Point& geom_pos);

private:
  std::shared_ptr<SplineCurve> geom_crv_;
  std::shared_ptr<SplineCurve> cond_crv_;
};

DirichletFunctor::DirichletFunctor(std::shared_ptr<SplineCurve> geom_crv,
				   std::shared_ptr<SplineCurve> cond_crv)
  : BdCondFunctor(), geom_crv_(geom_crv), cond_crv_(cond_crv)
{
}

DirichletFunctor::~DirichletFunctor()
{
}

Point DirichletFunctor::evaluate(const Point& geom_pos)
{
  double par, dist;
  Point close;
  geom_crv_->closestPoint(geom_pos, geom_crv_->startparam(),
			  geom_crv_->endparam(), par, close, dist);
  return cond_crv_->ParamCurve::point(par);
}

int main( int argc, char* argv[] )
{
  if (argc != 2) {
    cout << "Input parameters : Input file on g2 format" << endl;
    exit(-1);
  }

  // Read input arguments
  ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = 
    shared_ptr<CompositeModel>(factory.createFromG2(file1));

  shared_ptr<SurfaceModel> sfmodel = 
    dynamic_pointer_cast<SurfaceModel, CompositeModel>(model);

  vector<int> sol_dim(1);
  sol_dim[0] = 1;

  // Make block structured isogeometric surface model
  shared_ptr<IsogeometricSfModel> isomodel =
    shared_ptr<IsogeometricSfModel>(new IsogeometricSfModel(sfmodel,
							    sol_dim));

  vector<shared_ptr<IsogeometricSfBlock> > sf_blocks;
  isomodel->getIsogeometricBlocks(sf_blocks);

  int nmb_blocks = (int)sf_blocks.size();
  cout << "Number of blocks: " << nmb_blocks << endl;

  int nmb_bd = isomodel->getNmbOfBoundaries();
  CurveLoop bd = isomodel->getOuterBoundary();
  cout << "Number of boundaries: " << nmb_bd << endl;

  ofstream of1("model_bd.g2");
  int nmb_cvs = bd.size();
  int ki, kj;
  for (ki=0; ki<nmb_cvs; ++ki)
    {
      shared_ptr<ParamCurve> cv1 = bd[ki];
      shared_ptr<CurveOnSurface> cv2 = 
	dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv1);
      if (cv2.get())
	{
	  cv2->spaceCurve()->writeStandardHeader(of1);
	  cv2->spaceCurve()->write(of1);
	}
      else
	{
	  cv1->writeStandardHeader(of1);
	  cv1->write(of1);
	}
    }
	  
  vector<Point> joint_pnts = bd.getCorners();
  of1 << "400 1 0 4 255 0 0 255" << endl;
  of1 << joint_pnts.size() << endl;
  for (ki=0; ki<(int)joint_pnts.size(); ++ki)
    of1 << joint_pnts[ki] << endl;

  vector<Point> bd_pnts;
  int bd_stop = min(3, (int)joint_pnts.size());
  for (ki=0; ki<bd_stop; ++ki)
    bd_pnts.push_back(joint_pnts[ki]);
  isomodel->addBoundaryCond(0, bd_pnts, ZERO_NEUMANN, NULL, 0);

  vector<Point> bd_pnts2;
  int bd_stop2 = min(7, (int)joint_pnts.size());
  for (ki=bd_stop; ki<bd_stop2; ++ki)
    bd_pnts2.push_back(joint_pnts[ki]);

  double c_val = 1.0;
  isomodel->addBoundaryCond(0, bd_pnts2, CONSTANT_DIRICHLET, NULL, 0, &c_val);


  vector<Point> bd_pnts3;
  bd_pnts3.push_back(joint_pnts[joint_pnts.size()-1]);
  bd_pnts3.push_back(joint_pnts[0]);
  shared_ptr<SplineCurve> geom_crv = 
    shared_ptr<SplineCurve>(bd[nmb_cvs-1]->geometryCurve());
  Point pos1(1);
  Point pos2(1);
  pos1[0] = 0.5;
  pos2[0] = 1.0;
  shared_ptr<SplineCurve> cond_crv =
    shared_ptr<SplineCurve>(new SplineCurve(pos1, geom_crv->startparam(),
					    pos2, geom_crv->endparam()));
  shared_ptr<DirichletFunctor> dirfunc = 
    shared_ptr<DirichletFunctor>(new DirichletFunctor(geom_crv, cond_crv));
  isomodel->addBoundaryCond(0, bd_pnts3, DIRICHLET, dirfunc.get(), 0);

  int nmb_coef = sf_blocks[0]->nmbCoefs();
  int deg = 5;
  cout << "Block 0, number of coefs before raising to degree " << deg << ": " << nmb_coef << endl;

  isomodel->setMinimumDegree(deg, 0);
  nmb_coef = sf_blocks[0]->nmbCoefs();
  cout << "Block 0, number of coefs after degree raise: " << nmb_coef << endl;

  int nmb_neighbour = sf_blocks[0]->nmbOfNeighbours();
  cout << "Block 0, number of neighbours: " << nmb_neighbour << endl;

  int nmb_bd_cond = sf_blocks[0]->getNmbOfBoundaryConditions();
  cout << "Block 0, number of bd cond: " << nmb_bd_cond << endl;

  shared_ptr<SfSolution> sol = sf_blocks[0]->getSolutionSpace(0);
  cout << endl << "For next lines, Sol = solutions space 0 of block 0" << endl;
  cout << "Sol, solution degree: " << sol->degree(0) << " " << sol->degree(1) << endl;
  cout << "Sol, number of coefs: " << sol->nmbCoefs(0) << " " << sol->nmbCoefs(1) << endl;
  vector<double> knots1 = sol->knots(0);
  cout << "Sol, knots in dir 0: ";
  for (ki=0; ki<(int)knots1.size(); ++ki)
    {
      if (ki > 0)
	cout << ", ";
      cout << knots1[ki];
    }
  cout << endl;

  vector<double> knots2 = sol->distinctKnots(0);
  cout << "Sol, distinct knots in dir 0: ";
  for (ki=0; ki<(int)knots2.size(); ++ki)
    {
      if (ki > 0)
	cout << ", ";
      cout << knots2[ki];
    }
  cout << endl;

  vector<int> knot_intervals;
  knot_intervals.push_back(5);
  sol->insertKnots(knot_intervals, 0);

   knots2 = sol->distinctKnots(0);
  cout << "Sol, distinct knots in dir 0 after inserting knot in interval 5: ";
  for (ki=0; ki<(int)knots2.size(); ++ki)
    {
      if (ki > 0)
	cout << ", ";
      cout << knots2[ki];
    }
  cout << endl;

  vector<double> newknots;
  for (ki=1; ki<(int)knots2.size(); ++ki)
    newknots.push_back(0.5*(knots2[ki-1]+knots2[ki]));
  sol->insertKnots(newknots, 0);
  
   knots2 = sol->distinctKnots(0);
  cout << "Sol, distinct knots in dir 0 after inserting knot in every interval: ";
  for (ki=0; ki<(int)knots2.size(); ++ki)
    {
      if (ki > 0)
	cout << ", ";
      cout << knots2[ki];
    }
  cout << endl;

   isomodel->updateSolutionSplineSpace();

  for (ki=0; ki<nmb_blocks; ++ki)
    {
      cout << endl << "*********** DATA FOR BLOCK " << ki << " ***********" << endl << endl;
      sol = sf_blocks[ki]->getSolutionSpace(0);
      int nmb_bd_cond = sol->getNmbOfBoundaryConditions();
      cout << "Number of boundary conditions: " << nmb_bd_cond << endl;

      for (kj=0; kj<nmb_bd_cond; ++kj)
	{
	  shared_ptr<SfBoundaryCondition> cond = sol->getBoundaryCondition(kj);
	  BdConditionType bd_type = cond->getBdConditionType();
	  cout << "Condition " << kj << ", boundary condition type: " << bd_type << endl;
	  vector<int> enumeration;
	  cond->getCoefficientsEnumeration(enumeration);
	  cout << "Condition " << kj << ", boundary condition enumeration: ";
	  for (size_t kr=0; kr<enumeration.size(); ++kr)
	    {
	      if (kr > 0)
		cout << ", ";
	      cout << enumeration[kr];
	    }
	  cout << endl;

	  vector<pair<int, Point> > cond_coef;
	  cond->getBdCoefficients(cond_coef);
	  cout << "Condition " << kj << ", boundary coefficients:";
	  if (cond_coef.size() == 0)
	    cout << " (none)";
	  cout << endl;
	  for (size_t kr=0; kr<cond_coef.size(); ++kr)
	    cout << "   " << cond_coef[kr].first << " " << cond_coef[kr].second << endl;
	  cout << "Condition " << kj << ", edge number: " << cond->edgeNumber() << endl;

	  shared_ptr<SplineCurve> crv = cond->getSplineApproximation();
	  cout << "Condition " << kj << ", getSplineApproximation() pointer value: " << crv.get() << endl;
	  if (crv.get())
	    {
	      ofstream ofcond("bdcond_crv.g2");
	      crv->writeStandardHeader(ofcond);
	      crv->write(ofcond);
	    }
	}

      vector<vector<double> > knots(2);
      knots[0] = sol->distinctKnots(0);
      knots[1] = sol->distinctKnots(1);
      vector<vector<double> > param(2);
      for (kj=0; kj<2; ++kj)
	{
	  int deg = sol->degree(kj);
	  for (size_t kr=1; kr<knots[kj].size(); ++kr)
	    {
	      double del = (knots[kj][kr] - knots[kj][kr-1])/(double)(deg+1);
	      for (int kh=0; kh<deg; kh++)
		param[kj].push_back(knots[kj][kr-1]+(kh+1)*del);
	    }
	}
      sol->performPreEvaluation(param);
      vector<double> val;
      vector<double> der_u;
      vector<double> der_v;
      sol->getBasisFunctions(1, 0, val, der_u, der_v);
      cout << "Solution basis functions: " << endl;
      for (size_t kr=0; kr<val.size(); ++kr)
	cout << val[kr] << " ";
      cout << endl;
      cout << "Solution basis functions, derivatives in u-direction: " << endl;
      for (size_t kr=0; kr<der_u.size(); ++kr)
	cout << der_u[kr] << " ";
      cout << endl;
      cout << "Solution basis functions, derivatives in v-direction: " << endl;
      for (size_t kr=0; kr<der_v.size(); ++kr)
	cout << der_v[kr] << " ";
      cout << endl;

    for (kj=ki+1; kj<nmb_blocks; ++kj)
      {
	vector<int> edges, edges_other;
	vector<bool> equal_oriented;
	sf_blocks[ki]->getNeighbourInfo(sf_blocks[kj].get(),
					edges,
					edges_other,
					equal_oriented);

	cout << "Number of neighbouring edges of block " << ki << " and block " << kj << ": ";
	cout << edges.size() << endl;
	size_t kr;
	for (kr=0; kr<edges.size(); ++kr)
	  {
	    cout << "Neighbouring edge " << kr << " of block " << ki << " and block " << kj << ":" << endl;
	    cout << "   Block " << ki << ", edge number: " << edges[kr] << endl;
	    cout << "   Block " << kj << ", edge number: " << edges_other[kr] << endl;
	    cout << "   Equally oriented: " << (equal_oriented[kr] ? "Yes" : "No") << endl;
	  }

	if (sf_blocks[ki]->isNeighbour(sf_blocks[kj].get()))
	  {
	    shared_ptr<SfSolution> sol1 = sf_blocks[ki]->getSolutionSpace(0);
	    shared_ptr<SfSolution> sol2 = sf_blocks[kj]->getSolutionSpace(0);
	    bool matching = sol1->matchingSplineSpace(sol2.get());
	    cout << "Matching spline space for block " << ki << " and " << kj << " at solution 0: ";
	    cout << (matching ? "Yes" : "No") << endl;

	    vector<pair<int,int> > enumeration;
	    sol1->getMatchingCoefficients(sol2.get(), enumeration);
	    cout << "Matching coefs block " << ki << " and " << kj << " at solution 0:";
	    if (enumeration.size() == 0)
	      cout << " (none)";
	    cout << endl;
	    for (kr=0; kr<enumeration.size(); ++kr)
	      cout << "   " << enumeration[kr].first << " <-> " << enumeration[kr].second << endl;
	  }
      }
    }

}
