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
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/isogeometric_model/BdCondFunctor.h"
#include "GoTools/isogeometric_model/VolBoundaryCondition.h"


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

      Utils::eatwhite(file1);
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

  int nmb_blocks = (int)vol_blocks.size();
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
  int ki, kj, kh, kr, kq, kp, km, kk;
  for (ki = 0; ki < 6; ++ki)
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
  for (ki = 0; ki < (int)bd_sfs.size(); ++ki)
    {
      vector<pair<Point, Point> > corners;
      ParamSurface* sf = bd_sfs[ki].get();
      sf->getCornerPoints(corners);
      for (kj = 0; kj < (int)corners.size(); ++kj)
	  joint_pts.push_back(make_pair(sf, corners[kj].first));
      // ftSurface* face = dynamic_cast<ftSurface*>(corners[ki]->face());
      // assert(face != NULL);
      // shared_ptr<Vertex> vertex = corners[ki]->getVertex(true);
    }
  of1 << "400 1 0 4 255 0 0 255" << endl;
  of1 << joint_pts.size() << endl;
  for (ki=0; ki<(int)joint_pts.size(); ++ki)
    of1 << joint_pts[ki].second << endl;

  vector<pair<ParamSurface*, Point> > bd_pts;
  int bd_stop = min(3, (int)joint_pts.size());
  for (ki=0; ki<bd_stop; ++ki)
    bd_pts.push_back(joint_pts[ki]);
  isomodel->addBoundaryCond(bd_pts, ZERO_NEUMANN, NULL, 0);

  vector<pair<ParamSurface*, Point> > bd_pts2;
  int bd_stop2 = min(7, (int)joint_pts.size());
  for (ki=bd_stop; ki<bd_stop2; ++ki)
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

  int nmb_neighbour = vol_blocks[0]->nmbOfNeighbours();
  cout << "Block 0, number of neighbours: " << nmb_neighbour << endl;

  int nmb_bd_cond = vol_blocks[0]->getNmbOfBoundaryConditions();
  cout << "Block 0, number of bd cond: " << nmb_bd_cond << endl;

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
      sol = vol_blocks[ki]->getSolutionSpace(0);
      int nmb_bd_cond = sol->getNmbOfBoundaryConditions();
      cout << "Number of boundary conditions: " << nmb_bd_cond << endl;

      for (kj=0; kj<nmb_bd_cond; ++kj)
      {
	  shared_ptr<VolBoundaryCondition> cond = sol->getBoundaryCondition(kj);
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

	  vector<pair<int, Point> > cond_coef, cond_coef_bd, cond_coef_bd2;
	  cond->getBdCoefficients(cond_coef);
	  cond->getBdCoefficients(cond_coef_bd, cond_coef_bd2);
	  cout << "Condition " << kj << ", boundary coefficients:";
	  if (cond_coef.size() == 0)
	      cout << " (none)";
	  cout << endl;
	  for (size_t kr=0; kr<cond_coef.size(); ++kr)
	      cout << "   " << cond_coef[kr].first << " " << cond_coef[kr].second << endl;
	  cout << "Condition " << kj << ", face number: " << cond->faceNumber() << endl;

	  shared_ptr<SplineSurface> srf = cond->getSplineApproximation();
	  cout << "Condition " << kj << ", getSplineApproximation() pointer value: " << srf.get() << endl;
	  if (srf.get())
	  {
	      ofstream ofcond("bdcond_srf.g2");
	      srf->writeStandardHeader(ofcond);
	      srf->write(ofcond);
	  }
      }

      vector<vector<double> > knots(3);
      knots[0] = sol->distinctKnots(0);
      knots[1] = sol->distinctKnots(1);
      knots[2] = sol->distinctKnots(2);
      vector<vector<double> > param(3);
      for (kj=0; kj<2; ++kj)
      {
	  int deg = sol->degree(kj);
	  for (kr=1; kr<(int)knots[kj].size(); ++kr)
	  {
	      double del = (knots[kj][kr] - knots[kj][kr-1])/(double)(deg+1);
	      for (kh=0; kh<deg; kh++)
		  param[kj].push_back(knots[kj][kr-1]+(kh+1)*del);
	  }
      }
      sol->performPreEvaluation(param);
      // vector<double> val;
      // vector<double> der_u;
      // vector<double> der_v;
      // @@sbr Verify that order_u*order_v*order_w is indeed the
      // correct size of result arrays. Why not the sum since the
      // basis functions are independent?
      vector<double> val;
      vector<double> der_u;
      vector<double> der_v;
      vector<double> der_w;
//      shared_ptr<BasisDerivs> basis_derivs(new BasisDerivs());
      // sol->getBasisFunctions(1, 0, 0, basis_derivs);
      sol->getBasisFunctions(1, 0, 0, val, der_u, der_v, der_w);
      cout << "Solution basis functions: " << endl;
      // vector<double>& val = basis_derivs->basisValues;
      // vector<double>& der_u = basis_derivs->basisDerivs_u;
      // vector<double>& der_v = basis_derivs->basisDerivs_v;
      // vector<double>& der_w = basis_derivs->basisDerivs_w;
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
      cout << "Solution basis functions, derivatives in w-direction: " << endl;
      for (size_t kr=0; kr<der_w.size(); ++kr)
	  cout << der_w[kr] << " ";
      cout << endl;

#ifndef NDEBUG
      // We run through all basis functions, sum the evaluations and
      // derivs, compare against stored values.
      shared_ptr<SplineVolume> sol_vol =  sol->getSolutionVolume();
      const int num_coefs_u = sol_vol->numCoefs(0);
      const int num_coefs_v = sol_vol->numCoefs(1);
      const int num_coefs_w = sol_vol->numCoefs(2);
      const int dim = sol_vol->dimension();
      const int num_gauss_pts_u = (int)param[0].size();
      const int num_gauss_pts_v = (int)param[1].size();
      const int num_gauss_pts_w = (int)param[2].size();
      vector<double> global_val2(num_gauss_pts_u*num_gauss_pts_v*num_gauss_pts_w*dim, 0.0);
      vector<double> global_der2_u(num_gauss_pts_u*num_gauss_pts_v*num_gauss_pts_w*dim, 0.0);
      vector<double> global_der2_v(num_gauss_pts_u*num_gauss_pts_v*num_gauss_pts_w*dim, 0.0);
      for (kr = 0; kr < num_coefs_w; ++kr)
	  for (kq = 0; kq < num_coefs_v; ++kq)
	      for (kp = 0; kp < num_coefs_u; ++kp)
	      {
	      int basis_func_id_u = kp;//5;
	      int basis_func_id_v = kq;//5;
	      int basis_func_id_w = kr;//5;
	      vector<double> val2;
	      vector<double> der2_u;
	      vector<double> der2_v;
	      vector<double> der2_w;
//	      shared_ptr<BasisDerivs> basis_derivs2(new BasisDerivs());
	      vector<int> gauss_pts1, gauss_pts2, gauss_pts3;  // The id of quadrature points.
	      sol->getBasisFunctionValues(basis_func_id_u, basis_func_id_v, basis_func_id_w,
					  gauss_pts1, gauss_pts2, gauss_pts3,
					  val2, der2_u, der2_v, der2_w);
					  // basis_derivs2);
	      // vector<double>& val2 = basis_derivs2->basisValues;
	      // vector<double>& der2_u = basis_derivs2->basisDerivs_u;
	      // vector<double>& der2_v = basis_derivs2->basisDerivs_v;
	      // vector<double>& der2_w = basis_derivs2->basisDerivs_w;
					  // gauss_pts1, gauss_pts2,
					  // val2, der2_u, der2_v);
	      //puts("Done calling sol->getBasisFunctionValues().");

	      // For each gauss point we add the contribution to our
	      // global vector.
	      for (km = 0; km < (int)gauss_pts1.size(); ++km)
	      {
		  int gid_u = gauss_pts1[km];
		  int gid_v = gauss_pts2[km];
		  int gid_w = gauss_pts3[km];
		  for (kk = 0; kk < dim; ++kk)
		  {
		      global_val2[(gid_v*num_gauss_pts_u+gid_u)*dim+kk] += val2[km*dim+kk];
		      global_der2_u[(gid_v*num_gauss_pts_u+gid_u)*dim+kk] += der2_u[km*dim+kk];
		      global_der2_v[(gid_v*num_gauss_pts_u+gid_u)*dim+kk] += der2_v[km*dim+kk];
		  }
	      }
	  }

      // We then extract the values using the gauss point-interface.
      vector<double> global_val3(num_gauss_pts_u*num_gauss_pts_v*num_gauss_pts_w*dim, 0.0);
      vector<double> global_der3_u(num_gauss_pts_u*num_gauss_pts_v*num_gauss_pts_w*dim, 0.0);
      vector<double> global_der3_v(num_gauss_pts_u*num_gauss_pts_v*num_gauss_pts_w*dim, 0.0);
      vector<double> global_der3_w(num_gauss_pts_u*num_gauss_pts_v*num_gauss_pts_w*dim, 0.0);
      for (kr = 0; kr < num_gauss_pts_w; ++kr)
	  for (kq = 0; kq < num_gauss_pts_v; ++kq)
	      for (kp = 0; kp < num_gauss_pts_u; ++kp)
	  {
	      // shared_ptr<BasisDerivs> basis_derivs2(new BasisDerivs());
	      vector<double> val2;
	      vector<double> der2_u;
	      vector<double> der2_v;
	      vector<double> der2_w;
	      sol->getBasisFunctions(kp, kq, kr,
//				     basis_derivs2);
				     val2, der2_u, der2_v, der2_w);
	      // vector<double>& val2 = basis_derivs2->basisValues;
	      // vector<double>& der2_u = basis_derivs2->basisDerivs_u;
	      // vector<double>& der2_v = basis_derivs2->basisDerivs_v;
	      // vector<double>& der2_w = basis_derivs2->basisDerivs_w;
	      int num_vals = (int)val2.size()/dim;
	      for (km = 0; km < num_vals; ++km)
		  for (kk = 0; kk < dim; ++kk)
		  {
		      global_val3[(kq*num_gauss_pts_u+kp)*dim+kk] += val2[km*dim+kk];
		      global_der3_u[(kq*num_gauss_pts_u+kp)*dim+kk] += der2_u[km*dim+kk];
		      global_der3_v[(kq*num_gauss_pts_u+kp)*dim+kk] += der2_v[km*dim+kk];
		      global_der3_w[(kq*num_gauss_pts_u+kp)*dim+kk] += der2_w[km*dim+kk];
		  }
	  }
      double max_dist = 0.0, max_dist_u = 0.0, max_dist_v = 0.0;
      for (kq = 0; kq < (int)global_val2.size(); ++kq) // We do not bother about dim.
      {
	  double dist = fabs(global_val2[kq] - global_val3[kq]);
	  if (dist > max_dist)
	      max_dist = dist;
	  double dist_u = fabs(global_der2_u[kq] - global_der3_u[kq]);
	  if (dist_u > max_dist_u)
	      max_dist_u = dist_u;
	  double dist_v = fabs(global_der2_v[kq] - global_der3_v[kq]);
	  if (dist_v > max_dist_v)
	      max_dist_v = dist_v;
      }
      std::cout << "max_dist: " << max_dist << ", max_dist_u: " << max_dist_u << ", max_dist_v: " << max_dist_v << std::endl;
      puts("Done testing sol->getBasisFunction().");
#endif

      for (kj=ki+1; kj<nmb_blocks; ++kj)
      {
	  vector<int> faces, faces_other;
	  vector<int> orientation;
	  vector<bool> same_dir_order;
	  vol_blocks[ki]->getNeighbourInfo(vol_blocks[kj].get(),
					   faces,
					   faces_other,
					   orientation,
					   same_dir_order);

	  cout << "Number of neighbouring faces of block " << ki << " and block " << kj << ": ";
	  cout << faces.size() << endl;
	  size_t kr;
	  for (kr=0; kr<faces.size(); ++kr)
	  {
	      cout << "Neighbouring edge " << kr << " of block " << ki << " and block " << kj << ":" << endl;
	      cout << "   Block " << ki << ", edge number: " << faces[kr] << endl;
	      cout << "   Block " << kj << ", edge number: " << faces_other[kr] << endl;
	      cout << "   Oriented: " << orientation[kr] << endl;
	      // cout << "   Equally oriented: " << (equal_oriented[kr] ? "Yes" : "No") << endl;
	  }

	  if (vol_blocks[ki]->isNeighbour(vol_blocks[kj].get()))
	  {
	      shared_ptr<VolSolution> sol1 = vol_blocks[ki]->getSolutionSpace(0);
	      shared_ptr<VolSolution> sol2 = vol_blocks[kj]->getSolutionSpace(0);
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


  return -1;
}
