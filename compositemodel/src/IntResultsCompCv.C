#include "GoTools/compositemodel/IntResultsCompCv.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/tesselator/LineStrip.h"
#include "GoTools/tesselator/CurveTesselator.h"

using std::vector;

namespace Go
{
  //===========================================================================
  // Constructor
  IntResultsCompCv::IntResultsCompCv(CompositeCurve* compcv, const ftLine& line)
  //===========================================================================
  : IntResultsModel(CompositeCurve_Line)
  {
    addLineInfo(line);
  }

  //===========================================================================
  // Constructor
  IntResultsCompCv::IntResultsCompCv(CompositeCurve* compcv, const ftPlane& plane)
  //===========================================================================
  : IntResultsModel(CompositeCurve_Plane)
  {
    addPlaneInfo(plane);
  }

  //===========================================================================
  // Destructor
  IntResultsCompCv::~IntResultsCompCv()
  //===========================================================================
  {
  }

  //===========================================================================
  void IntResultsCompCv::addIntPt(shared_ptr<ParamCurve> cv, double* parval)
  //===========================================================================
  {
    if (numpar_ == 1)
      int_pts_1cv_.push_back(PointOnCurve(cv, *parval));
  }


  //===========================================================================
  void IntResultsCompCv::addIntCv(shared_ptr<ParamCurve> cv, double* startpar,
			     double *endpar)
  //===========================================================================
  {
    if (numpar_ == 1)
      int_seg_1cv_.push_back(std::make_pair(PointOnCurve(cv, *startpar),
				       PointOnCurve(cv, *endpar)));
  }

  //===========================================================================
  void 
  IntResultsCompCv::getIntersectionPoints(std::vector<PointOnCurve>& int_points) const
  //===========================================================================
  {
    int_points = int_pts_1cv_;
  }

  //===========================================================================
  void 
  IntResultsCompCv::getIntersectionCurves(std::vector<std::pair<PointOnCurve, PointOnCurve> >& int_crvs) const
  //===========================================================================
  {
    int_crvs = int_seg_1cv_;
  }

  //===========================================================================
  void IntResultsCompCv::tesselate(std::vector<shared_ptr<LineStrip> >& meshes,
				   PointCloud3D& points) const
  //===========================================================================
  {
    int res = 100;
    tesselate(res, meshes, points);
  }

  //===========================================================================
  void IntResultsCompCv::tesselate(int resolution,
				   std::vector<shared_ptr<LineStrip> >& meshes,
				   PointCloud3D& points) const
  //===========================================================================
  {
    if (hasIntCurves())
      {
	for (size_t ki=0; ki<int_seg_1cv_.size(); ++ki)
	  {
	    double t1 = int_seg_1cv_[ki].first.getPar();
	    double t2 = int_seg_1cv_[ki].second.getPar();
	    shared_ptr<ParamCurve> cv = int_seg_1cv_[ki].first.getCurve();
	    shared_ptr<ParamCurve> sub_cv = 
	      shared_ptr<ParamCurve>(cv->subCurve(std::min(t1,t2),
						  std::max(t1,t2)));
	    
	    CurveTesselator tesselator(*sub_cv.get());
	    tesselator.changeRes(resolution);
	    shared_ptr<LineStrip> mesh = tesselator.getMesh();
	    meshes.push_back(mesh);
	  }
      }

    vector<double> coords;
    for (size_t ki=0; ki<int_pts_1cv_.size(); ++ki)
      {
	Point pt = int_pts_1cv_[ki].getPos();
	coords.insert(coords.end(), pt.begin(), pt.end());
      }

    points = PointCloud3D(&coords[0], (int)int_pts_1cv_.size());
  }

  //===========================================================================
  void IntResultsCompCv::tesselate(double density,
				   std::vector<shared_ptr<LineStrip> >& meshes,
				   PointCloud3D& points) const
  //===========================================================================
  {
    int min_nmb = 5;
    int max_nmb = (int)(1000000.0/(int)int_seg_1cv_.size());
    if (hasIntCurves())
      {
	for (size_t ki=0; ki<int_seg_1cv_.size(); ++ki)
	  {
	    double t1 = int_seg_1cv_[ki].first.getPar();
	    double t2 = int_seg_1cv_[ki].second.getPar();
	    shared_ptr<ParamCurve> cv = int_seg_1cv_[ki].first.getCurve();
	    shared_ptr<ParamCurve> sub_cv = 
	      shared_ptr<ParamCurve>(cv->subCurve(std::min(t1,t2),
						  std::max(t1,t2)));
	    double len = sub_cv->estimatedCurveLength();
	    int res = (int)(len/density);
	    res = std::max(min_nmb, std::min(res, max_nmb));
	    
	    CurveTesselator tesselator(*sub_cv.get());
	    tesselator.changeRes(res);
	    shared_ptr<LineStrip> mesh = tesselator.getMesh();
	    meshes.push_back(mesh);
	  }
      }

    vector<double> coords;
    for (size_t ki=0; ki<int_pts_1cv_.size(); ++ki)
      {
	Point pt = int_pts_1cv_[ki].getPos();
	coords.insert(coords.end(), pt.begin(), pt.end());
      }

    points = PointCloud3D(&coords[0], (int)int_pts_1cv_.size());
  }

} // namespace Go
