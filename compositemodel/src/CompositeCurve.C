//===========================================================================
//                                                                           
// File: CompositeCurve.C                                                   
//                                                                           
// Created: May 2009
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/compositemodel/CompositeCurve.h"
#include "GoTools/utils/CurvatureUtils.h"
#include "GoTools/tesselator/CurveTesselator.h"
#include "GoTools/tesselator/TesselatorUtils.h"
#include <fstream>

using std::vector;

namespace Go
{


  //===========================================================================
  CompositeCurve::CompositeCurve(double gap,   // Gap between adjacent curves
				 double neighbour,  // Threshold for whether curves are adjacent
				 double kink,  // Kink between adjacent curves
				 double bend, // Intended G1 discontinuity between adjacent curves
				 vector<shared_ptr<ParamCurve> >& curves)
  //===========================================================================
    : CompositeModel(gap, neighbour, kink, bend)
  {
    if (curves.empty())
      return;

    curves_.resize(curves.size());
    perm_.resize(curves.size());
    turned_.resize(curves.size());
    continuity_.resize(curves.size());
    start_parameters_.resize(curves.size());

    for (size_t ki=0; ki<curves.size(); ++ki)
      {
	curves_[ki] = curves[ki];
	perm_[ki] = (int)ki;
	turned_[ki] = false;
	continuity_[ki] = JOINT_NONE;
      }

    orderCurves();
  }

  //===========================================================================
  CompositeCurve::~CompositeCurve()
  //===========================================================================
  {
    // Empty destructor
  }

  //===========================================================================
  CompositeCurve* CompositeCurve::clone() const
  //===========================================================================
  {
    vector<shared_ptr<ParamCurve> > curves(curves_.size());
    for (size_t ki=0; ki<curves_.size(); ++ki)
      curves[ki] = shared_ptr<ParamCurve>(curves_[ki]->clone());

    return new CompositeCurve(toptol_.gap, toptol_.neighbour,
			      toptol_.kink, toptol_.bend, curves);
  }

  //===========================================================================
  int CompositeCurve::nmbEntities() const
  //===========================================================================
  {
      return (int)curves_.size();
  }

  //===========================================================================
  shared_ptr<ParamCurve> CompositeCurve::getCurve(int idx) const
  //===========================================================================
  {
    shared_ptr<ParamCurve> result;
    if (idx >= 0 && idx< (int)curves_.size())
      return curves_[idx];
    else
      return result;
  }

  //===========================================================================
  int CompositeCurve::getIndex(ParamCurve* curve) const
  //===========================================================================
  {
    for (int ki=0; ki<(int)curves_.size(); ++ki)
      {
	if (curves_[ki].get() == curve)
	  return ki;
      }
    return -1;  // Curve not found
  }

  //===========================================================================
  void CompositeCurve::parameterRange(double& start, double& end) const
  //===========================================================================
  {
    start = start_parameters_[perm_[0]];
    size_t ncrv = curves_.size();
    end = start_parameters_[perm_[ncrv-1]] + 
      (curves_[perm_[ncrv-1]]->endparam() - curves_[perm_[ncrv-1]]->startparam());
  }

  //===========================================================================
  double CompositeCurve::getGlobalPar(int idx,  // Index of sub curve
				      double par) const  // Parameter in sub curve
  //===========================================================================
  {
    return start_parameters_[idx] + par - curves_[idx]->startparam();
  }

  //===========================================================================
  void CompositeCurve::getLocalPar(double global_par,  // Parameter in composite curve
				     int& idx,         // Index of corresponding sub curve
				     double& local_par) const  // Parameter in sub curve
  //===========================================================================
  {
    size_t ki;
    for (ki=1; ki<curves_.size(); ++ki)
      if (start_parameters_[perm_[ki]] > global_par)
	break;

    idx = perm_[ki-1];

    local_par = global_par - start_parameters_[idx] + 
      curves_[idx]->startparam();
  }

  //===========================================================================
  void 
  CompositeCurve::evaluate(int idx,      // Index
			   double par[], // Parameter value
			   Point& pnt) const
  //===========================================================================
  {
    curves_[idx]->point(pnt, par[0]);
  }

  //===========================================================================
  void 
  CompositeCurve::evaluate(int idx,      // Index
			   double par[], // Parameter value
			   int nder,     // Number of derivatives to compute, 0=only position
			   std::vector<Point>& der) const
  //===========================================================================
  {
    curves_[idx]->point(der, par[0], nder);
  }

  //===========================================================================
  void 
  CompositeCurve::evaluateCurve(double par, // Parameter value
				Point& pnt) const  // Result
  //===========================================================================
  {
    int idx;
    double local_par;
    getLocalPar(par, idx, local_par);
    curves_[idx]->point(pnt, local_par);
  }

  //===========================================================================
  void 
  CompositeCurve::evaluateCurve(double par, // Parameter value
				int nder,  // Number of derivatives to compute, 0=only position
				std::vector<Point>& der) const  // Result
  //===========================================================================
  {
    int idx;
    double local_par;
    getLocalPar(par, idx, local_par);
    curves_[idx]->point(der, local_par, nder);
  }

//===========================================================================
void 
CompositeCurve::closestPoint(Point& pnt,     // Input point
			     Point& clo_pnt, // Found closest point
			     int& idx,           // Index of curve where the closest point is found
			     double clo_par[],   // Parameter value corresponding to the closest point
			     double& dist)       // Distance between input point and found closest point
//===========================================================================
  {
    PointOnCurve closest = closestPoint(pnt);
    clo_pnt = closest.getPos();
    clo_par[0] = closest.getPar();
    idx = getIndex(closest.getCurve().get());
    dist = clo_pnt.dist(pnt);
    closest_idx_ = idx;
  }

//===========================================================================
PointOnCurve
CompositeCurve::closestPoint(Point& pnt)  
//===========================================================================
  {
    PointOnCurve result;
    Point cp;
    double dist, par;
    double best_dist = 1e100; // A gogool should be enough

    // If a previous closest point is found, compute the closest point corresponding
    // to that curve to get a new best closest point
    if (closest_idx_ >= 0)
      {
	curves_[closest_idx_]->closestPoint(pnt, par, cp, best_dist);
	result = PointOnCurve(curves_[closest_idx_], par);
      }

    // Traverse all curve to find the best closest point
    for (size_t ki=0; ki<curves_.size(); ++ki)
      {
	if ((int)ki == closest_idx_)
	  continue;  // Already computed

	// Check the distance to the current boundingbox
	BoundingBox box = curves_[ki]->boundingBox();
	double box_dist = boxVecDist(box, pnt);
	if (box_dist > best_dist)
	  continue;  // No possibility for a better closest point

	// Compute current closest point
	curves_[ki]->closestPoint(pnt, par, cp, dist);
	if (dist < best_dist)
	  {
	    best_dist = dist;
	    result = PointOnCurve(curves_[ki], par);
	  }
      }
	    
    return result;
  }



//===========================================================================
BoundingBox CompositeCurve::boundingBox() 
//===========================================================================
{
  BoundingBox box;
  if (curves_.size() == 0)
    return box;

  box = curves_[0]->boundingBox();
  for (size_t ki=0; ki<curves_.size(); ++ki)
    {
      BoundingBox bb = curves_[ki]->boundingBox();
      box.addUnionWith(bb);
    }

  return box;
}

//===========================================================================
BoundingBox CompositeCurve::boundingBox(int idx) const
//===========================================================================
{
  return curves_[idx]->boundingBox();
}

//===========================================================================
bool CompositeCurve::isDegenerate(int idx) const
//===========================================================================
{
  if (idx < 0 || idx >= (int)curves_.size())
    return true;  // Illegal index
  return curves_[idx]->isDegenerate(toptol_.gap);
}

//===========================================================================
double CompositeCurve::curvature(int idx, // Index of curve
				 double *par) const
//===========================================================================
{
  vector<Point> der(3), unitder(3);
  evaluate(idx, par, 2, der);
  (void)curvatureRadius(der, unitder);
  double curvature = unitder[2].length();
  return curvature;
}

//===========================================================================
void CompositeCurve::turn(int idx)
//===========================================================================
{
  curves_[idx]->reverseParameterDirection();
  turned_[idx] = (turned_[idx]) ? false : true;
}

//===========================================================================
void CompositeCurve::turn()
//===========================================================================
{
  for (size_t ki=0; ki<curves_.size(); ++ki)
    {
      curves_[ki]->reverseParameterDirection();
      perm_[ki] = (int)ki;
      turned_[ki] = false;
      continuity_[ki] = JOINT_NONE;
    }

    orderCurves();
}

//===========================================================================
void CompositeCurve::append(shared_ptr<ParamCurve> curve)
//===========================================================================
{
  // Add curve to the curve vector
  int nmb_crvs = (int)curves_.size();
  curves_.push_back(curve);
  continuity_.push_back(JOINT_DISC);
  turned_.push_back(false);

  // Find the position in the curve sequence
  int ki, kj;
  vector<Point> pt1(2), pt2(2);
  curve->point(pt1, curve->startparam(), 1);
  curve->point(pt2, curve->endparam(), 1);

  for (ki=0; ki<nmb_crvs; ki++)
    {
      kj = (ki == 0) ? nmb_crvs-1 : ki-1;

      vector<Point> pt3(2), pt4(2);
      double t3 = (turned_[perm_[kj]]) ? curves_[perm_[kj]]->startparam() :
	curves_[perm_[kj]]->endparam();
      double t4 = (turned_[perm_[ki]]) ? curves_[perm_[ki]]->endparam() :
	curves_[perm_[ki]]->startparam();
      curves_[perm_[kj]]->point(pt3, t3, 1);
      curves_[perm_[ki]]->point(pt4, t4, 1);
      
      // If it already is a smooth joint, do not try to place the curve here
      if (continuity_[perm_[kj]] == JOINT_G1)
	continue;

      // Find continuity status at current position
      double d1 = pt3[0].dist(pt1[0]);
      double d2 = pt4[0].dist(pt2[0]);
      double d3 = pt3[0].dist(pt2[0]);
      double d4 = pt4[0].dist(pt1[0]);
      bool turned1 = false, turned2 = false;

      tpJointType joint1=JOINT_DISC, joint2=JOINT_DISC;
      if (d1 <= toptol_.neighbour || d3 <= toptol_.neighbour)
	{
	  joint1 = JOINT_GAP;
	  if (d3 < d1)
	    turned1 = true;
	  if ((d1 <= toptol_.gap && (!turned1)) || 
	      (d3 <= toptol_.gap && turned1))
	    {
	      joint1 = JOINT_G0;
	      double ang = (turned1) ? pt3[1].angle(pt2[1]) :
		pt3[1].angle(pt1[1]);

	      if (ang <= toptol_.kink)
		joint1 = JOINT_G1;
	      else if (ang <= toptol_.bend)
		joint1 = JOINT_KINK;
	    }
	}

      if (d2 <= toptol_.neighbour || d4 <= toptol_.neighbour)
	{
	  joint2 = JOINT_GAP;
	  if (d4 < d2)
	    turned2 = true;
	  if ((d2 <= toptol_.gap && (!turned2)) || 
	      (d4 <= toptol_.gap && turned2))
	    {
	      joint2 = JOINT_G0;
	      double ang = (turned2) ? pt4[1].angle(pt2[1]) :
		pt4[1].angle(pt1[1]);

	      if (ang <= toptol_.kink)
		joint2 = JOINT_G1;
	      else if (ang <= toptol_.bend)
		joint2 = JOINT_KINK;
	    }
	}

      if (!(joint1 <= continuity_[perm_[kj]] && 
	    joint2 <= continuity_[perm_[kj]]))
	continue;  // Not the correct position

      if ((turned1 && (!turned2)) || (turned2 && !(turned1)))
	{
	  // Inconsistent direction. Must update continuity
	  if (joint1 <= joint2)
	    {
	      turned2 = turned1;
	      double dist = (turned2) ? d4 : d2;
	      double ang = (turned2) ? pt4[1].angle(pt2[1]) :
		pt4[1].angle(pt1[1]);
	      if (dist <= toptol_.gap)
		{
		  if (ang <= toptol_.kink)
		    joint2 = JOINT_G1;
		  else if (ang <= toptol_.bend)
		    joint2 = JOINT_KINK;
		  else
		    joint2 = JOINT_G0;
		}
	      else if (dist <= toptol_.neighbour)
		joint2 = JOINT_GAP;
	      else
		joint2 = JOINT_DISC;
		}
	  else
	    {
	      turned1 = turned2;
	      double dist = (turned1) ? d3 : d1;
	      double ang = (turned1) ? pt3[1].angle(pt2[1]) :
		pt3[1].angle(pt1[1]);
	      if (dist <= toptol_.gap)
		{
		  if (ang <= toptol_.kink)
		    joint1 = JOINT_G1;
		  else if (ang <= toptol_.bend)
		    joint1 = JOINT_KINK;
		  else
		    joint1 = JOINT_G0;
		}
	      else if (dist <= toptol_.neighbour)
		joint1 = JOINT_GAP;
	      else
		joint1 = JOINT_DISC;
	    }
	}

      if (!(joint1 <= continuity_[perm_[kj]] && 
	    joint2 <= continuity_[perm_[kj]]))
	continue;  // Not the correct position

      bool insert_here = false;
      if (joint1 < continuity_[perm_[kj]] && 
	  continuity_[perm_[kj]] >= JOINT_DISC)
	insert_here = true;
      else if (joint2 < continuity_[perm_[kj]] && 
	       continuity_[perm_[kj]] >= JOINT_DISC)
	insert_here = true;
      else if (joint1+joint2 < 2*continuity_[perm_[kj]])
	insert_here = true;
	
	if (insert_here)
	  {
	    continuity_[perm_[kj]] = joint1;
	    continuity_[nmb_crvs] = joint2;
	    turned_[nmb_crvs] = turned1;
	    perm_.insert(perm_.begin()+ki,nmb_crvs);
	    break; // No need to search anymore
	  } 
    }

  if (ki == nmb_crvs)
    {
      // Add curve at the end
      turned_[ki] = false;
      continuity_[perm_[ki-1]] = JOINT_DISC;
      continuity_[ki] = JOINT_DISC;
      perm_.push_back(ki);
    }

  // Update start parameters
  start_parameters_.resize(curves_.size());
  double startpar = 0.0;
  for (size_t kh=0; kh<curves_.size(); ++kh)
    {
      double tdel = curves_[perm_[kh]]->endparam() - 
	curves_[perm_[kh]]->startparam();
      start_parameters_[perm_[kh]] = startpar;
      startpar += tdel;
    }
}

//===========================================================================
void CompositeCurve::tesselate(vector<shared_ptr<GeneralMesh> >& meshes) const
//===========================================================================
{
  int res = 100;
  tesselate(&res, meshes);
}

//===========================================================================
  void CompositeCurve::tesselate(int resolution[],
				 vector<shared_ptr<GeneralMesh> >& meshes) const
//===========================================================================
{
  tesselate(curves_, resolution, meshes);
}

//===========================================================================
  void CompositeCurve::tesselate(const vector<shared_ptr<ParamCurve> >& curves,
				 int resolution[],
				 vector<shared_ptr<GeneralMesh> >& meshes) const
//===========================================================================
{
  meshes.clear();
  for (size_t ki=0; ki<curves.size(); ++ki)
    {
      CurveTesselator tesselator(*curves[ki].get());
      tesselator.changeRes(resolution[0]);
      shared_ptr<GeneralMesh> mesh = tesselator.getMesh();
      meshes.push_back(mesh);
    }
}

//===========================================================================
  void CompositeCurve::tesselate(double density,
				 vector<shared_ptr<GeneralMesh> >& meshes) const
//===========================================================================
{
  tesselate(curves_, density, meshes);
}

//===========================================================================
  void CompositeCurve::tesselate(const vector<shared_ptr<ParamCurve> >& curves,
				 double density,
				 vector<shared_ptr<GeneralMesh> >& meshes) const
//===========================================================================
{
    int min_nmb = 5;
    int max_nmb = (int)(1000000.0/(int)curves.size());

    for (size_t ki=0; ki<curves.size(); ++ki)
      {
	double len = curves[ki]->estimatedCurveLength();
	int res = (int)(len/density);
	res = std::max(min_nmb, std::min(res, max_nmb));

	CurveTesselator tesselator(*curves[ki].get());
	tesselator.changeRes(res);
	shared_ptr<GeneralMesh> mesh = tesselator.getMesh();
	meshes.push_back(mesh);
      }

}

  //===========================================================================
  void CompositeCurve::tesselatedCtrPolygon(vector<shared_ptr<LineCloud> >& ctr_pol) const
  //===========================================================================
  {
    tesselatedCtrPolygon(curves_, ctr_pol);
  }

  //===========================================================================
  void CompositeCurve::tesselatedCtrPolygon(const vector<shared_ptr<ParamCurve> >& curves,
					    vector<shared_ptr<LineCloud> >& ctr_pol) const
  //===========================================================================
  {
    for (size_t ki=0; ki<curves.size(); ++ki)
      {
	shared_ptr<LineCloud> curr_pol = TesselatorUtils::getCtrPol(curves[ki].get());
	ctr_pol.push_back(curr_pol);
      }
  }

   //===========================================================================
  void CompositeCurve::orderCurves()
  //===========================================================================
  {
    size_t ki, kj;  // Counters
    size_t last = 0;    // The last ordered curve
    double startpar = 0.0;

    while (last < curves_.size())
      {
	// Look for a curve with a free end
	bool atstart, atend;
	for (ki=last; ki<curves_.size(); ++ki)
	  {
	    atstart = atend = false;
	    Point pt1, pt2;
	    curves_[perm_[ki]]->point(pt1, curves_[perm_[ki]]->startparam());
	    curves_[perm_[ki]]->point(pt2, curves_[perm_[ki]]->endparam());
	    for (kj=last; kj<curves_.size(); ++kj)
	      {
		Point pt3, pt4;
		if (kj == ki)
		  continue;

		curves_[perm_[kj]]->point(pt3, curves_[perm_[kj]]->startparam());
		curves_[perm_[kj]]->point(pt4, curves_[perm_[kj]]->endparam());

		if (pt1.dist(pt3) < toptol_.neighbour || 
		    pt1.dist(pt4) < toptol_.neighbour)
		  atstart = true;

		if (pt2.dist(pt3) < toptol_.neighbour || 
		    pt2.dist(pt4) < toptol_.neighbour)
		  atend = true;

		if (atstart && atend)
		  break;
	      }

	    if (kj < curves_.size())
	      continue;  // This curve is not an end curve

	    // The curve is found
	    if (atstart && !atend)
	      turned_[perm_[ki]] = true;
	    start_parameters_[perm_[ki]] = startpar;
	    startpar += (curves_[perm_[ki]]->endparam() - 
			 curves_[perm_[ki]]->startparam());
	    std::swap(perm_[ki], perm_[last]);

	    break;  // A start curve is found
	  }

	if (!atstart && !atend)
	  {
	    continuity_[perm_[last]] = JOINT_DISC;
	    last++;
	    continue;  // An isolated curve, look for a new start curve
	  }
	    
	for (ki=last; ki<curves_.size(); ++ki)
	  {
	    vector<Point> pt1(2), pt2(2);
	    curves_[perm_[ki]]->point(pt1, curves_[perm_[ki]]->startparam(), 1);
	    curves_[perm_[ki]]->point(pt2, curves_[perm_[ki]]->endparam(), 1);
	    if (turned_[perm_[ki]])
	      std::swap(pt1, pt2);

	    for (kj=ki+1; kj<curves_.size(); ++kj)
	      {
		vector<Point> pt3(2), pt4(2);
		curves_[perm_[kj]]->point(pt3, curves_[perm_[kj]]->startparam(), 1);
		curves_[perm_[kj]]->point(pt4, curves_[perm_[kj]]->endparam(), 1);
		if (turned_[perm_[kj]])
		  std::swap(pt3, pt4);

		if (pt1[0].dist(pt3[0]) <= toptol_.neighbour)
		  {
		    std::ofstream out_file("comp_crv.g2");
		    for (size_t k4=0; k4<curves_.size(); ++k4)
		      {
			curves_[k4]->writeStandardHeader(out_file);
			curves_[k4]->write(out_file);
		      }

		    MESSAGE("Possible multiple connectivity in composite curve");
		  }
		if (pt2[0].dist(pt3[0]) <= toptol_.neighbour)
		  {
		    if (pt2[0].dist(pt3[0]) > toptol_.gap)
		      continuity_[perm_[ki]] = JOINT_GAP;
		    else if (pt2[1].angle(pt3[1]) > toptol_.bend)
		      continuity_[perm_[ki]] = JOINT_G0;
		    else if (pt2[1].angle(pt3[1]) > toptol_.kink)
		      continuity_[perm_[ki]] = JOINT_KINK;
		    else
		      continuity_[perm_[ki]] = JOINT_G1;

		    start_parameters_[perm_[kj]] = startpar;
		    startpar += (curves_[perm_[kj]]->endparam() - 
				 curves_[perm_[kj]]->startparam());

		    std::swap(perm_[ki+1],perm_[kj]);
		    break;  // The next curve is found
		  }

		if (pt2[0].dist(pt4[0]) <= toptol_.neighbour)
		  {
		    turned_[perm_[kj]] = true;
		    if (pt2[0].dist(pt4[0]) > toptol_.gap)
		      continuity_[perm_[ki]] = JOINT_GAP;
		    else if (pt2[1].angle(pt4[1]) > toptol_.bend)
		      continuity_[perm_[ki]] = JOINT_G0;
		    else if (pt2[1].angle(pt4[1]) > toptol_.kink)
		      continuity_[perm_[ki]] = JOINT_KINK;
		    else
		      continuity_[perm_[ki]] = JOINT_G1;

		    start_parameters_[perm_[kj]] = startpar;
		    startpar += (curves_[perm_[kj]]->endparam() - 
			 curves_[perm_[kj]]->startparam());
		    std::swap(perm_[ki+1],perm_[kj]);
		    break;
		  }
	      }

	    if (kj == curves_.size())
	      {
		// The end of a chain is found.
		// Check for a closed loop
		vector<Point> pt1(2), pt2(2);
		curves_[perm_[last]]->point(pt1, (turned_[perm_[last]]) ?
					    curves_[perm_[last]]->endparam() :
					    curves_[perm_[last]]->startparam(), 1);
		curves_[perm_[ki]]->point(pt2, (turned_[perm_[ki]]) ?
					  curves_[perm_[ki]]->startparam() :
					  curves_[perm_[ki]]->endparam(), 1);
		if (pt1[0].dist(pt2[0]) > toptol_.neighbour)
		  continuity_[perm_[ki]] = JOINT_DISC;
		else if (pt1[0].dist(pt2[0]) > toptol_.gap)
		  continuity_[perm_[ki]] = JOINT_GAP;
		else if (pt1[1].angle(pt2[1]) > toptol_.bend)
		  continuity_[perm_[ki]] = JOINT_G0;
		else if (pt1[1].angle(pt2[1]) > toptol_.kink)
		  continuity_[perm_[ki]] = JOINT_KINK;
		    else
		      continuity_[perm_[ki]] = JOINT_G1;

		// Look for a new start point
		last = ki+1;
		break;
	      }
	  }
      }
  }

}  // namespace Go
