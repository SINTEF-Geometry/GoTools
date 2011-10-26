#include "GoTools/compositemodel/IntResultsSfModel.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/tesselator/LineStrip.h"

using namespace std;
using std::shared_ptr;

namespace Go
{
  //===========================================================================
  // Constructor
  IntResultsSfModel::IntResultsSfModel(SurfaceModel* sfmodel, const ftLine& line)
  //===========================================================================
    : IntResultsModel(SurfaceModel_Line), sfmodel1_(sfmodel)
  {
    addLineInfo(line);
  }

  //===========================================================================
  // Constructor
  IntResultsSfModel::IntResultsSfModel(SurfaceModel* sfmodel, const ftPlane& plane)
  //===========================================================================
    : IntResultsModel(SurfaceModel_Plane), sfmodel1_(sfmodel)
  {
    addPlaneInfo(plane);
  }

  //===========================================================================
  // Destructor
  IntResultsSfModel::~IntResultsSfModel()
  //===========================================================================
  {
  }

  //===========================================================================
  void IntResultsSfModel::addIntPts(vector<ftPoint>& intpts)
  //===========================================================================
  {
    int_pts_.insert(int_pts_.end(), intpts.begin(), intpts.end());
  }


  //===========================================================================
  void IntResultsSfModel::addIntCvs(ftCurve& cvs)
  //===========================================================================
  {
    intcvs_ = cvs;
  }

  //===========================================================================
  void IntResultsSfModel::tesselate(std::vector<std::shared_ptr<LineStrip> >& meshes,
				    PointCloud3D& points) const
  //===========================================================================
  {
    if (hasIntCurves())
      intcvs_.tesselate(meshes);
    
    vector<double> coords;
    for (size_t ki=0; ki<int_pts_.size(); ++ki)
      {
	const Point pt = int_pts_[ki].position();
	coords.insert(coords.end(), pt.begin(), pt.end());
      }

    points = PointCloud3D(&coords[0], (int)int_pts_.size());
  }

  //===========================================================================
  void IntResultsSfModel::tesselate(int resolution, 
				    std::vector<std::shared_ptr<LineStrip> >& meshes,
				    PointCloud3D& points) const
  //===========================================================================
  {
    if (hasIntCurves())
      intcvs_.tesselate(resolution, meshes);
    
    vector<double> coords;
    for (size_t ki=0; ki<int_pts_.size(); ++ki)
      {
	const Point pt = int_pts_[ki].position();
	coords.insert(coords.end(), pt.begin(), pt.end());
      }

    points = PointCloud3D(&coords[0], (int)int_pts_.size());
  }

  //===========================================================================
  void IntResultsSfModel::tesselate(double density, 
				    std::vector<std::shared_ptr<LineStrip> >& meshes,
				    PointCloud3D& points) const
  //===========================================================================
  {
    if (hasIntCurves())
      intcvs_.tesselate(density, meshes);
    
    vector<double> coords;
    for (size_t ki=0; ki<int_pts_.size(); ++ki)
      {
	const Point pt = int_pts_[ki].position();
	coords.insert(coords.end(), pt.begin(), pt.end());
      }

    points = PointCloud3D(&coords[0], (int)int_pts_.size());
  }

} // namespace Go
