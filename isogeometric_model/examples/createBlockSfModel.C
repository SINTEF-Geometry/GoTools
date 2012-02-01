//===========================================================================
//                                                                           
// File: createBlockSfModel.C
//                                                                           
//===========================================================================

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
using std::min;

using namespace Go;

// This class mimics the behaviour of a class representing
// a dirichlet boundary condition. It is a subclass of
// BdCondFunctor and can be evaluated with respect to
// a point in geometry space.
// This class must be provided by the application using the
// isogeometric surface model.
class DirichletFunctor : public BdCondFunctor
{
public:
  DirichletFunctor(shared_ptr<SplineCurve> geom_crv,
		   shared_ptr<SplineCurve> cond_crv);
  virtual ~DirichletFunctor();
  virtual Point evaluate(const Point& geom_pos);

private:
  // The boundary curve of the surface set corresponding to the given
  // boundary conditions
  shared_ptr<SplineCurve> geom_crv_; 
  // The dirichlet condition itself. The two curves must have corresponding
  // parameterization, but that falls naturally from the fact that the
  // geometry and the solution fields are represented in the same spline
  // space up to possible refinements in the fields.
  shared_ptr<SplineCurve> cond_crv_;
};

DirichletFunctor::DirichletFunctor(shared_ptr<SplineCurve> geom_crv,
				   shared_ptr<SplineCurve> cond_crv)
  : BdCondFunctor(), geom_crv_(geom_crv), cond_crv_(cond_crv)
{
}

DirichletFunctor::~DirichletFunctor()
{
}

// Evaluation is performed by computing the associated parameter value
// to the input point on the geometry curve by a closest point computation.
// The value of the dirichlet condition are found in the same parameter
// value
Point DirichletFunctor::evaluate(const Point& geom_pos)
{
  double par, dist;
  Point close;
  geom_crv_->closestPoint(geom_pos, geom_crv_->startparam(),
			  geom_crv_->endparam(), par, close, dist);
  return cond_crv_->ParamCurve::point(par);
}

//===========================================================================
//                                                                           
// Description:
//  
// This program illustrates parts of the interface of the class 
// IsogeometricSfModel. This class holds a block structured isogeometric
// 2-variate model and provides an interface for use in isogometric
// analysis. 
//
// This class has not yet been used in an isogeometric analysis environment.
// Thus, the usability of the interface has not been proved. Bug reports, 
// feedback on the interface and suggestions for improvements are very
// welcome.
//
// Input/Output
// The surface set is given in the file data/surface_set.g2.
// Results from enquiry functionality are written to file during the 
// execution. 
//                                                                           
//===========================================================================

int main( int argc, char* argv[] )
{
  // Specify input file
  ifstream file1("data/surface_set.g2");

  // Prepare for output files
  std::string outfile1("data/model_bd.g2");
  std::string outfile2("data/bdcond_crv.g2");

  std::cout << "Preparing input file" << std::endl;

  // Topology tolerances for the surface set. See the example file
  // compositemodel/examples/face2splineSet.C for an explanation
  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.001;
  double bend = 0.01;
  double approxtol = 0.01;

  // Create engine for reading the surface set
  CompositeModelFactory factory(approxtol, gap, neighbour, kink, bend);

  // Read model and represent it as a surface model
  shared_ptr<CompositeModel> model = 
    shared_ptr<CompositeModel>(factory.createFromG2(file1));

  shared_ptr<SurfaceModel> sfmodel = 
    dynamic_pointer_cast<SurfaceModel, CompositeModel>(model);
  if (!sfmodel.get())
    {
      std::cout << "The input file did not contain a surface set" << std::endl;
      exit(-1);
    }

  // Check if the input data is appropriate for isogeometric analysis
    bool isOK = sfmodel->allSplines();
    if (!isOK)
      {
	std::cout << "Not all surfaces are splines. Stopping computation" << std::endl;
	exit(-1);
      }

    // Ensure common spline spaces and corresponding coefficients
    sfmodel->makeCommonSplineSpaces();

  std::cout << "Creating the isogeometric surafce model" << std::endl;
    
  // One solution field will be associated this surface model
  vector<int> sol_dim(1);
  // The field is a scalar field, i.e. has dimension 1
  sol_dim[0] = 1;

  // Make block structured isogeometric surface model
  shared_ptr<IsogeometricSfModel> isomodel =
    shared_ptr<IsogeometricSfModel>(new IsogeometricSfModel(sfmodel,
							    sol_dim));

  std::cout << "Enquiry " << std::endl;
  std::ofstream of("data/sfmodel_info.txt");

  // Fetch information about the number of blocks in this model,
  // and fetch the surface blocks
  vector<shared_ptr<IsogeometricSfBlock> > sf_blocks;
  isomodel->getIsogeometricBlocks(sf_blocks);

  int nmb_blocks = (int)sf_blocks.size();
  of << "Number of blocks: " << nmb_blocks << endl;

  // Fetch information about model boundaries
  int nmb_bd = isomodel->getNmbOfBoundaries();
  of << "Number of boundaries: " << nmb_bd << endl;
  vector<CurveLoop> bd(nmb_bd);
  bd[0] = isomodel->getOuterBoundary();

  for (size_t kr=1; kr<nmb_bd; ++kr)
    bd[kr] = isomodel->getBoundary(kr);

  // Fetch the points associated to each model boundary specifying
  // the joints between block boundaries in the model boundary
  vector<vector<Point>> joint_pnts(nmb_bd);
  for (size_t kr=0; kr<nmb_bd; ++kr)
    joint_pnts[kr] = bd[kr].getCorners();

  // Write model boundaries and associated points
  ofstream of1(outfile1.c_str());
  for (size_t kr=0; kr<nmb_bd; ++kr)
    {
      int nmb_cvs = bd[kr].size();
      int ki, kj;
      for (ki=0; ki<nmb_cvs; ++ki)
	{
	  shared_ptr<ParamCurve> cv1 = bd[kr][ki];
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
	  
      of1 << "400 1 0 4 0 255 0 255" << endl;
      of1 << 1 << endl;
      of1 << joint_pnts[kr][0] << endl;
      of1 << "400 1 0 4 255 0 0 255" << endl;
      of1 << joint_pnts[kr].size()-1 << endl;
      for (ki=1; ki<(int)joint_pnts[kr].size(); ++ki)
	of1 << joint_pnts[kr][ki] << endl;
    }

  // Set some boundary conditions. These conditions are intended to illustrate
  // the possibilities, not be realistic in an analysis situation
  // The first boundary condition is along one boundary in the model which
  // is made from two block boundaries and which starts from the second
  // joint point.
  vector<Point> bd_pnts1;
  int bd_stop = min(4, (int)joint_pnts[0].size());
  int ki;
  for (ki=1; ki<bd_stop; ++ki)
    bd_pnts1.push_back(joint_pnts[0][ki]);
  isomodel->addBoundaryCond(0, bd_pnts1, ZERO_NEUMANN, NULL, 0);

  // The next boundary condition is related to the next model boundary,
  // again two blocks. As the joints are only once, we start the sequence 
  // of points with the point with which we ended the last sequence
  vector<Point> bd_pnts2;
  int bd_stop2 = min(6, (int)joint_pnts[0].size());
  for (ki=bd_stop-1; ki<bd_stop2; ++ki)
    bd_pnts2.push_back(joint_pnts[0][ki]);

  // Constant dirichlet with value 1
  double c_val = 1.0;
  isomodel->addBoundaryCond(0, bd_pnts2, CONSTANT_DIRICHLET, NULL, 0, &c_val);

  // The third boundary condition is again along a model boundary. This
  // time a linear Dirichlet condition is applied
  vector<Point> bd_pnts3;
  int bd_stop3 = min(8, (int)joint_pnts[0].size());
  for (ki=bd_stop2-1; ki<bd_stop3; ++ki)
    bd_pnts3.push_back(joint_pnts[0][ki]);
  
  // Feth the boundary curves associated with the current boundary condition
  shared_ptr<SplineCurve> geom_crv1 = 
    shared_ptr<SplineCurve>(bd[0][6]->geometryCurve());
  shared_ptr<SplineCurve> geom_crv2 = 
    shared_ptr<SplineCurve>(bd[0][7]->geometryCurve());

  // Join curves into one, create copy of 1. curve 
  double dist;   // Will be zero as no non-existing continuity is enforced
  shared_ptr<SplineCurve> geom_crv(geom_crv1->clone());
  geom_crv->appendCurve(geom_crv2.get(), 0, dist, false);  // C0 continuity, 
                                                     // no reparametrization

  // Create linear curve representing the Dirichlet condition
  Point pos1(1);
  Point pos2(1);
  pos1[0] = 1.0;
  pos2[0] = 0.0;
  shared_ptr<SplineCurve> cond_crv =
    shared_ptr<SplineCurve>(new SplineCurve(pos1, geom_crv->startparam(),
  					    pos2, geom_crv->endparam()));
  shared_ptr<DirichletFunctor> dirfunc = 
    shared_ptr<DirichletFunctor>(new DirichletFunctor(geom_crv, cond_crv));
  isomodel->addBoundaryCond(0, bd_pnts3, DIRICHLET, dirfunc.get(), 0);

  // The last piece of the outer boundary of the model. Zero Neumann
  // The possible choices can be found in BdConditionType.h. The
  // symmetry condition is not implemented
  vector<Point> bd_pnts4;
  bd_pnts4.push_back(joint_pnts[0][joint_pnts[0].size()-1]);
  bd_pnts4.push_back(joint_pnts[0][0]);
  bd_pnts4.push_back(joint_pnts[0][1]);
  isomodel->addBoundaryCond(0, bd_pnts4, ZERO_NEUMANN, NULL, 0);
  
  // The inner boundary curve is simply set to be zero Dirichlet
  // Note that the corner points fetch from the model boundary is not a
  // closed loop. Thus the first point must be repeated
  vector<Point> bd_pnts5(joint_pnts[1].begin(), joint_pnts[1].end());
  bd_pnts5.push_back(joint_pnts[1][0]);
  isomodel->addBoundaryCond(1, bd_pnts5, ZERO_DIRICHLET, NULL, 0);

  // Raise degree of the solution space correponsing to the surface model
  // Note that then number of coefficient in the geometric block does not
  // change
  int nmb_coef = sf_blocks[0]->nmbCoefs();
  int deg = 5;
  of << "Block 0, number of coefs before raising to degree " << deg << ": " << nmb_coef << endl;

  isomodel->setMinimumDegree(deg, 0);
  nmb_coef = sf_blocks[0]->nmbCoefs();
  of << "Block 0, number of coefs after degree raise: " << nmb_coef << endl;

  // Fetch information regarding the first surface block
  int nmb_neighbour = sf_blocks[0]->nmbOfNeighbours();
  of << "Block 0, number of neighbours: " << nmb_neighbour << endl;

  int nmb_bd_cond = sf_blocks[0]->getNmbOfBoundaryConditions();
  of << "Block 0, number of bd cond: " << nmb_bd_cond << endl;

  shared_ptr<SfSolution> sol = sf_blocks[0]->getSolutionSpace(0);
  of << endl << "For next lines, Sol = solutions space 0 of block 0" << endl;
  of << "Sol, solution degree: " << sol->degree(0) << " " << sol->degree(1) << endl;
  of << "Sol, number of coefs: " << sol->nmbCoefs(0) << " " << sol->nmbCoefs(1) << endl;

  // Fetch information about the spline space of the solution field
  // corresponding to the first surface block
  // All knot values
  vector<double> knots1 = sol->knots(0);
  of << "Sol, knots in dir 0: ";
  for (ki=0; ki<(int)knots1.size(); ++ki)
    {
      if (ki > 0)
  	of << ", ";
      of << knots1[ki];
    }
  of << endl;

  // The distinct knot values
  vector<double> knots2 = sol->distinctKnots(0);
  of << "Sol, distinct knots in dir 0: ";
  for (ki=0; ki<(int)knots2.size(); ++ki)
    {
      if (ki > 0)
  	of << ", ";
      of << knots2[ki];
    }
  of << endl;

  // Insert one new knot in the middle of each knot interval in the first
  // parameter directions. This functionality can also be used to insert
  // knots at arbitrary positions within the knot interval, but this	
  // opportunity should be used with care as it easily leads to large
  // data amounts when the spline spaces are of the solution fields
  // of the various blocks are harmonized
  vector<double> newknots;
  for (ki=1; ki<(int)knots2.size(); ++ki)
    newknots.push_back(0.5*(knots2[ki-1]+knots2[ki]));
  sol->insertKnots(newknots, 0);
  
   knots2 = sol->distinctKnots(0);
  of << "Sol, distinct knots in dir 0 after inserting knot in every interval: ";
  for (ki=0; ki<(int)knots2.size(); ++ki)
    {
      if (ki > 0)
  	of << ", ";
      of << knots2[ki];
    }
  of << endl;

  // Insert one new knot in the second knot interval in the first parameter 
  // direction of the solution space. The knot will automatically be inserted
  // in the middle of the knot interval
  vector<int> knot_intervals;
  knot_intervals.push_back(1);
  sol->insertKnots(knot_intervals, 0);

   knots2 = sol->distinctKnots(0);
  of << "Sol, distinct knots in dir 0 after inserting knot in interval 1 (counting from 0): ";
  for (ki=0; ki<(int)knots2.size(); ++ki)
    {
      if (ki > 0)
  	of << ", ";
      of << knots2[ki];
    }
  of << endl;

  // Ensure that the blocks corresponding to the solution space have consistent
  // spline spaces at common block boundaries. The geometry model is not changed
   isomodel->updateSolutionSplineSpace();

   // Fetch information from each block
   // Curves representing (possibly an approximation) given boundary
   // conditions are written to one file regardless of associated 
   // surface block and bondary condition
   ofstream of2(outfile2.c_str());
  for (ki=0; ki<nmb_blocks; ++ki)
    {
      of << endl << "*********** DATA FOR BLOCK " << ki << " ***********" << endl << endl;
      of << "Boundary condition information" << endl << endl;

      sol = sf_blocks[ki]->getSolutionSpace(0);
      int nmb_bd_cond = sol->getNmbOfBoundaryConditions();
      of << "Number of boundary conditions: " << nmb_bd_cond << endl;

      // For each boundary conditions, fetch the associated information
      for (int kj=0; kj<nmb_bd_cond; ++kj)
  	{
	  // Boundary condition type:
	  // 0 = unknown, 1 = zero Dirichlet, 2 = constant Dirichlet, 
	  // 3 = Dirichlet, 4 = zero Neumann, 5 = Neumann,
	  // 6 = symmetry (not implemented)
  	  shared_ptr<SfBoundaryCondition> cond = sol->getBoundaryCondition(kj);
  	  BdConditionType bd_type = cond->getBdConditionType();
  	  of << "Condition " << kj << ", boundary condition type: " << bd_type << endl;

	  // Fetch the number of the solution coefficients affected by this
	  // boundary condition, i.e. the coefficients at the solution surface
	  // at the part of the boundary where this condition lives.
	  // The coefficient enumeration runs as follows: (x1, y1), (x2, y1), 
	  // ... (xn, y1), (x1, y2), (x2, y2), ..., (x1, ym), ... (xn, ym)
  	  vector<int> enumeration;
  	  cond->getCoefficientsEnumeration(enumeration);
  	  of << "Condition " << kj << ", boundary condition enumeration: ";
  	  for (size_t kr=0; kr<enumeration.size(); ++kr)
  	    {
  	      if (kr > 0)
  		of << ", ";
  	      of << enumeration[kr];
  	    }
  	  of << endl;

	  // Only in the Diriclet case. Fetch the value of the coefficients of
	  // the solution field corresponding to the current boundary 
	  // condition. The enumeration of the coefficient is given first,
	  // then the value of the coefficient (of the solution field)
  	  vector<pair<int, Point> > cond_coef;
  	  cond->getBdCoefficients(cond_coef);
  	  of << "Condition " << kj << ", boundary coefficients:";
  	  if (cond_coef.size() == 0)
  	    of << " (none)";
  	  of << endl;
  	  for (size_t kr=0; kr<cond_coef.size(); ++kr)
  	    of << "   " << cond_coef[kr].first << " " << cond_coef[kr].second << endl;
	  // The edge number of the solution field (and the surface block)
	  // corresponding to this boundary condition. 0=umin, 1=umax,
	  // 2=vmin, 3=vmax
  	  of << "Condition " << kj << ", edge number: " << cond->edgeNumber() << endl;

  	  shared_ptr<SplineCurve> crv = cond->getSplineApproximation();
  	  of << "Condition " << kj << ", getSplineApproximation() pointer value: " << crv.get() << endl;
  	  if (crv.get())
  	    {
  	      crv->writeStandardHeader(of2);
  	      crv->write(of2);
  	    }
  	}

      of << endl << "Basis function evaluation " << endl << endl;

      // Grid evaluation of the basis functions. First create a parameter
      // grid in which to evaluate. In a real application this will be
      // the Gauss points corresponding to the quadrature formula.
      vector<vector<double> > knots(2);
      knots[0] = sol->distinctKnots(0);
      knots[1] = sol->distinctKnots(1);
      vector<vector<double> > param(2);
      for (int kj=0; kj<2; ++kj)
  	{
  	  int deg = sol->degree(kj);
  	  for (size_t kr=1; kr<knots[kj].size(); ++kr)
  	    {
  	      double del = (knots[kj][kr] - knots[kj][kr-1])/(double)(deg+1);
  	      for (int kh=0; kh<deg; kh++)
  		param[kj].push_back(knots[kj][kr-1]+(kh+1)*del);
  	    }
  	}
      // Perform grid avaluation of the solution field corresponding to
      // the current block
      sol->performPreEvaluation(param);

      // Fetch information from the pre evaluated basis functions which are
      // stored together with the solution field
      // The information is fetched from the Gauss point corresponding to 
      // the second given parameter in the first parameter directions and
      // the first in the second parameter direction. 
      // The degree is in this case 5, thus it is 6 non-zero basis functions
      // in each parameter direction. Multiplying, we expect 36 non-zero
      // basis functions
      vector<double> val;
      vector<double> der_u;
      vector<double> der_v;
      sol->getBasisFunctions(1, 0, val, der_u, der_v);
      of << "Solution basis functions in the Gauss point (1,0): " << endl;
      for (size_t kr=0; kr<val.size(); ++kr)
  	of << val[kr] << " ";
      of << endl;
      of << "Solution basis functions, derivatives in u-direction: " << endl;
      for (size_t kr=0; kr<der_u.size(); ++kr)
  	of << der_u[kr] << " ";
      of << endl;
      of << "Solution basis functions, derivatives in v-direction: " << endl;
      for (size_t kr=0; kr<der_v.size(); ++kr)
  	of << der_v[kr] << " ";
      of << endl;

      // Compute the Jacobian in the same Gauss point
      vector<int> index(2);
      index[0] = 1;
      index[1] = 0;
      double jacobian = sol->getJacobian(index);
      of << "Jacobian: " << jacobian << endl;

      of << endl << "Information regarding block adjacency" << endl << endl;
    for (int kj=ki+1; kj<nmb_blocks; ++kj)
      {
	// For each combination of blocks, check if the two blocks are neighbours
	// and in that case along which boundaries. Do also compute whether these
	// boundaries have the same orientations
	// boundaries: the enumeration of the common boundaries in the
	// first block (sf_blocks[ki]). The interpretation is as follows
	// 0 - start parameter in the first parameter direction, 1 - end parameter
	// in the first parameter direction, 2 - start parameter in the
	// second parameter direction, 3 - end parameter in the secocond 
	// parameter direction
	// boundaries_other: the same information with regard to sf_blocks[kj]
	// equal_oriented: whether or not the corresponding boundary curves 
	// have the same orientation in the two surfaces
  	vector<int> boundaries, boundaries_other;
  	vector<bool> equal_oriented;
  	sf_blocks[ki]->getNeighbourInfo(sf_blocks[kj].get(),
  					boundaries,
  					boundaries_other,
  					equal_oriented);

  	of << "Number of adjacent boundaries between block " << ki << " and block " << kj << ": ";
  	of << boundaries.size() << endl;
  	size_t kr;
  	for (kr=0; kr<boundaries.size(); ++kr)
  	  {
  	    of << "Adjacent edge " << kr << " of block " << ki << " and block " << kj << ":" << endl;
  	    of << "   Block " << ki << ", boundary number: " << boundaries[kr] << endl;
  	    of << "   Block " << kj << ", boundary number: " << boundaries_other[kr] << endl;
  	    of << "   Equally oriented: " << (equal_oriented[kr] ? "Yes" : "No") << endl;
  	  }

	// For the current combination of blocks, check if the two 
	// blocks are neighbours. In this case, we could alternatively
	// check that the length of the boundaries vector is larger than zero
  	if (sf_blocks[ki]->isNeighbour(sf_blocks[kj].get()))
  	  {
	    // The two blocks are neighbours, fetch the solution fields 
	    // related to these blocks
  	    shared_ptr<SfSolution> sol1 = sf_blocks[ki]->getSolutionSpace(0);
  	    shared_ptr<SfSolution> sol2 = sf_blocks[kj]->getSolutionSpace(0);

	    // Check if the spline spaces spaces of these solution fields
	    // match along the common boundary
  	    bool matching = sol1->matchingSplineSpace(sol2.get());
  	    of << "Matching spline space for block " << ki << " and " << kj << " at solution 0: ";
  	    of << (matching ? "Yes" : "No") << endl;

	    // Get the coefficient enumeration of the two solution fields at the
	    // common boundary. The enumeration corresponding to each field is 
	    // pairwise represented in the vector enumeration. If the two spline
	    // space do not match, the enumeration vector will have zero length
  	    vector<pair<int,int> > enumeration;
  	    sol1->getMatchingCoefficients(sol2.get(), enumeration);
  	    of << "Matching coefs block " << ki << " and " << kj << " at solution 0:";
  	    if (enumeration.size() == 0)
  	      of << " (none)";
  	    of << endl;
  	    for (kr=0; kr<enumeration.size(); ++kr)
  	      of << "   " << enumeration[kr].first << " <-> " << enumeration[kr].second << endl;
  	  }
      }
    }

}
