//===========================================================================
//                                                                           
// File: testVolModel.C                                                      
//                                                                           
// Created: Wed Jun 15 15:46:41 2011                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include <iostream>
#include <fstream>

#include "GoTools/isogeometric_model/IsogeometricVolModel.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
// #include "GoTools/trivariate/GoToolsTrivariate.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/isogeometric_model/BdCondFunctor.h"


using std::vector;
using std::cout;
using std::endl;

using std::ifstream;
using std::ofstream;
using std::min;
using std::make_pair;
using std::pair;

using namespace Go;

class DirichletFunctor : public BdCondFunctor
{
public:
  DirichletFunctor(shared_ptr<SplineCurve> geom_crv,
		   shared_ptr<SplineCurve> cond_crv);
  virtual ~DirichletFunctor();
  virtual Point evaluate(const Point& geom_pos);

private:
  shared_ptr<SplineCurve> geom_crv_;
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
  ifstream file1(argv[1]); // Volumes.
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  ofstream of1("tmp/model_bd.g2");

  vector<shared_ptr<ftVolume> > volumes;

  while (!file1.eof())
    {
      // Read volume from file
      ObjectHeader head;
      file1 >> head;
      shared_ptr<SplineVolume> vol2;
      vol2 = shared_ptr<SplineVolume>(new SplineVolume());
      vol2->read(file1);

      shared_ptr<ParamVolume> pvol
          = dynamic_pointer_cast<ParamVolume, SplineVolume>(vol2);
      volumes.push_back(shared_ptr<ftVolume>(new ftVolume(pvol)));
      // volumes.push_back(shared_ptr<ftVolume>(new ftVolume(vol2)));

      eatwhite(file1);
    }

  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double bend = 0.05;
  // double approxtol = 0.01;
  shared_ptr<VolumeModel> volmodel =
    shared_ptr<VolumeModel>(new VolumeModel(volumes, gap, neighbour, 
					    kink, bend));

  // GoToolsTrivariate::init();
  // GoToolsTrivariate gotools_tri;
  // gotools_tri.init();

  if (volmodel.get() == NULL)
  {
    MESSAGE("File does not contain volumes, exiting!");
    return -1;
  }

  vector<int> sol_dim(1);
  sol_dim[0] = 1;

  // Make block structured isogeometric volume model
  shared_ptr<IsogeometricVolModel> isomodel =
    shared_ptr<IsogeometricVolModel>(new IsogeometricVolModel(volmodel,
							      sol_dim));

  vector<shared_ptr<IsogeometricVolBlock> > vol_blocks;
  isomodel->getIsogeometricBlocks(vol_blocks);

  size_t nmb_blocks = vol_blocks.size();
  cout << "Number of (volume) blocks: " << nmb_blocks << endl;

  int nmb_bd = isomodel->getNmbOfBoundaries();
  vector<shared_ptr<ParamSurface> > bd_sfs = isomodel->getOuterBoundary();
  cout << "Number of boundaries: " << nmb_bd << endl;

  // If bd was not extracted there is nothing more to do.
  if (bd_sfs.size() == 0)
    {
      MESSAGE("No boundary detected, some code is missing!");
      return 0;
    }

  // We write to file the outer boundary.
  vector<shared_ptr<ParamCurve> > bd_cvs;
  for (int ki = 0; ki < 6; ++ki)
    {
      shared_ptr<ParamSurface> bd_sf = bd_sfs[ki];
      if (bd_sf.get() != NULL)
	{
	  CurveLoop cv_loop = bd_sf->outerBoundaryLoop();
	  bd_cvs.insert(bd_cvs.end(), cv_loop.begin(), cv_loop.end());
	  bd_sf->writeStandardHeader(of1);
	  bd_sf->write(of1);
	}
    }

  int nmb_cvs = (int)bd_cvs.size();

  // We impose some bd constraints, using the corner points of the
  // boundary.
  vector<pair<ParamSurface*, Point> > joint_pts;
  for (size_t ki = 0; ki < bd_sfs.size(); ++ki)
    {
      vector<pair<Point, Point> > corners;
      ParamSurface* sf = bd_sfs[ki].get();
      sf->getCornerPoints(corners);
      for (size_t kj = 0; kj < corners.size(); ++kj)
	  joint_pts.push_back(make_pair(sf, corners[kj].first));
      // ftSurface* face = dynamic_cast<ftSurface*>(corners[ki]->face());
      // assert(face != NULL);
      // shared_ptr<Vertex> vertex = corners[ki]->getVertex(true);
    }
  of1 << "400 1 0 4 255 0 0 255" << endl;
  of1 << joint_pts.size() << endl;
  for (int ki=0; ki<(int)joint_pts.size(); ++ki)
    of1 << joint_pts[ki].second << endl;

  vector<pair<ParamSurface*, Point> > bd_pts;
  int bd_stop = min(3, (int)joint_pts.size());
  for (int ki=0; ki<bd_stop; ++ki)
    bd_pts.push_back(joint_pts[ki]);
  isomodel->addBoundaryCond(bd_pts, ZERO_NEUMANN, NULL, 0);

  vector<pair<ParamSurface*, Point> > bd_pts2;
  int bd_stop2 = min(7, (int)joint_pts.size());
  for (int ki=bd_stop; ki<bd_stop2; ++ki)
    bd_pts2.push_back(joint_pts[ki]);

  double c_val = 1.0;
  isomodel->addBoundaryCond(bd_pts2, CONSTANT_DIRICHLET, NULL, 0, &c_val);

  // Set (non-constant) Dirichlet condition!
  vector<pair<ParamSurface*, Point> > bd_pts3;
  bd_pts3.push_back(joint_pts[joint_pts.size()-1]);
  bd_pts3.push_back(joint_pts[0]);
  shared_ptr<SplineCurve> geom_crv = 
    shared_ptr<SplineCurve>(bd_cvs[nmb_cvs-1]->geometryCurve());
  Point pos1(1);
  Point pos2(1);
  pos1[0] = 0.5;
  pos2[0] = 1.0;
  shared_ptr<SplineCurve> cond_crv =
    shared_ptr<SplineCurve>(new SplineCurve(pos1, geom_crv->startparam(),
					    pos2, geom_crv->endparam()));
  shared_ptr<DirichletFunctor> dirfunc = 
    shared_ptr<DirichletFunctor>(new DirichletFunctor(geom_crv, cond_crv));
  isomodel->addBoundaryCond(bd_pts3, DIRICHLET, dirfunc.get(), 0);

  // We raise the degree (to later verify that constraints still hold).
  int nmb_coef = vol_blocks[0]->nmbCoefs();
  int deg = 5;
  cout << "Block 0, number of coefs before raising to degree " << deg << ": " << nmb_coef << endl;
  isomodel->setMinimumDegree(deg, 0);
  nmb_coef = vol_blocks[0]->nmbCoefs();
  cout << "Block 0, number of coefs after degree raise: " << nmb_coef << endl;

  // We refine spline space (to later verify that constraints still hold).
  shared_ptr<VolSolution> sol = vol_blocks[0]->getSolutionSpace(0);
  if (sol.get() == NULL)
    {
      MESSAGE("No solution exists, code missing I suppose!");
      return -1;
    }

  cout << endl << "For next lines, Sol = solutions space 0 of block 0" << endl;
  cout << "Sol, solution degree: " << sol->degree(0) << " " << sol->degree(1) << endl;
  cout << "Sol, number of coefs: " << sol->nmbCoefs(0) << " " << sol->nmbCoefs(1) << endl;
  vector<double> knots1 = sol->knots(0);
  cout << "Sol, knots in dir 0: ";
  for (int ki=0; ki<(int)knots1.size(); ++ki)
    {
      if (ki > 0)
	cout << ", ";
      cout << knots1[ki];
    }
  cout << endl;

  vector<double> knots2 = sol->distinctKnots(0);
  cout << "Sol, distinct knots in dir 0: ";
  for (int ki=0; ki<(int)knots2.size(); ++ki)
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
  for (int ki=0; ki<(int)knots2.size(); ++ki)
    {
      if (ki > 0)
	cout << ", ";
      cout << knots2[ki];
    }
  cout << endl;

  vector<double> newknots;
  for (int ki=1; ki<(int)knots2.size(); ++ki)
    newknots.push_back(0.5*(knots2[ki-1]+knots2[ki]));
  sol->insertKnots(newknots, 0);
  
   knots2 = sol->distinctKnots(0);
  cout << "Sol, distinct knots in dir 0 after inserting knot in every interval: ";
  for (int ki=0; ki<(int)knots2.size(); ++ki)
    {
      if (ki > 0)
	cout << ", ";
      cout << knots2[ki];
    }
  cout << endl;

  isomodel->updateSolutionSplineSpace();

  return -1;
}
