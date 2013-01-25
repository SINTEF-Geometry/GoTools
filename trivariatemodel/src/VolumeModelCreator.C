#include "GoTools/trivariatemodel/VolumeModelCreator.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/Body.h"
#include "GoTools/compositemodel/RegularizeFace.h"
#include "GoTools/compositemodel/RegularizeFaceSet.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/Array.h"
#include "GoTools/utils/MatrixXD.h"
#include <fstream>

//#define DEBUG

using std::vector;
using namespace Go;

//===========================================================================
bool VolumeModelCreator::createRotationalModel(shared_ptr<SurfaceModel>& sfmodel,
					       shared_ptr<VolumeModel>& volmodel)
//===========================================================================
{
  // Check if the surface model represents a solid
  int nmb_bd = sfmodel->nmbBoundaries();
  if (nmb_bd > 0)
    return false;  // Not a solid

  // Create solid
  shared_ptr<Body> body(new Body(sfmodel));

  // Check if the surface model represents a rotational object
  Point centre, axis, vec;
  double angle;
  double min_ang;
  bool rotational = sfmodel->isAxisRotational(centre, axis, vec, angle,
					      min_ang);
  if (!rotational)
    return false;

  double eps = sfmodel->getTolerances().gap;  // Tolerance
  
  // Compute a section of the model. Check first if the position should be changed
  Point dir;
  //if (angle < 2.0*M_PI-eps)
  if (min_ang < 2.0*M_PI-eps)
    {
      // Move section to avoid coincident intersections
      Array<double,3> tmp_vec(vec[0], vec[1], vec[2]);
      MatrixXD<double, 3> mat;
      // mat.setToRotation(0.5*angle, axis[0], 
      // 			axis[1], axis[2]);  // Rotate the 
      mat.setToRotation(0.25*min_ang, axis[0], 
			axis[1], axis[2]);  // Rotate the 
      // start vector the angle 0.5*angle around axis
      Array<double,3> tmp_vec2 = mat*tmp_vec;
      dir = Point(tmp_vec2[0], tmp_vec2[1], tmp_vec2[2]);
    }
  else
    dir = vec;

  // Create a planar B-spline surface from the rotational axis to a point on
  // the model boundary
  // First compute intersections between the bounding box of the model and the
  // computed axes
  BoundingBox box = sfmodel->boundingBox();
  vector<Point> result;
  vector<Point> res = box.lineIntersect(centre, axis);
  if (res.size() > 0)
    result.insert(result.end(), res.begin(), res.end());
  res.clear();
  res = box.lineIntersect(centre, -1.0*axis);
  if (res.size() > 0)
    result.insert(result.end(), res.begin(), res.end());
  if (result.size() != 2)
    return false;  // Unexpected result
  res.clear();
  Point mid = 0.5*(result[0]+result[1]);  // Point internal in the box
  res = box.lineIntersect(mid, dir);
  if (res.size() > 0)
    result.insert(result.end(), res.begin(), res.end());
  
  // Remove duplicates
  for (size_t ki=0; ki<result.size(); ++ki)
    for (size_t kj=ki+1; kj<result.size(); ++kj)
      {
	if (result[ki].dist(result[kj]) < eps)
	  {
	    result.erase(result.begin()+kj);
	    kj--;
	  }
      }

  if (result.size() != 3)
    return false;  // Unexpected result
    
  // Create planar B-spline surface interpolating these three points
  Point vec1 = result[1] - result[0];
  shared_ptr<SplineCurve> cv1(new SplineCurve(result[0]-0.1*vec1,
					      result[1]+0.1*vec1));

  Point vec2 = result[2] - mid;
  shared_ptr<SplineCurve> cv2(new SplineCurve(mid, result[2]+0.1*vec2));

  SweepSurfaceCreator sweep;
  shared_ptr<SplineSurface> surf(sweep.linearSweptSurface(*cv1, *cv2, mid));

  // Make surface model
  vector<shared_ptr<ParamSurface> > sfs;
  sfs.push_back(surf);

  double neighbour = sfmodel->getTolerances().neighbour;
  double kink = sfmodel->getTolerances().kink;
  double bend = sfmodel->getTolerances().bend;
  shared_ptr<SurfaceModel> model2(new SurfaceModel(eps, eps, neighbour, 
						   kink, 10.0*kink, sfs));

  // Intersect the two models to get the part of the plane internal to the
  // initial model
  vector<shared_ptr<SurfaceModel> > submodels = 
    sfmodel->splitSurfaceModels(model2);
 
  // Check which part of the plane is internal. The plane is given by the third and fourth
  // sub model. First fetch a point in the third submodel
  double upar, vpar;
  if (submodels[2]->nmbEntities() != 1 || 
      submodels[3]->nmbEntities() != 1)
    return false;  // Something went wrong
  Point pt = submodels[2]->getSurface(0)->getInternalPoint(upar, vpar);

  // Check if this point is inside the solid
  bool inside = body->isInside(pt);
  int idx = (inside) ? 2 : 3;
  
  // Fetch face
  shared_ptr<ftSurface> face = submodels[idx]->getFace(0);

#ifdef DEBUG
  std::ofstream of1("section_sf.g2");
  face->surface()->writeStandardHeader(of1);
  face->surface()->write(of1);
#endif

  // Regularize face
  RegularizeFace reg(face, eps, neighbour, bend, true);
  Point centre2;  // Dummy
  reg.setAxis(centre2, axis);
  vector<shared_ptr<ftSurface> > reg_faces = reg.getRegularFaces();
  if (reg_faces.size() == 0)
    return false;

#ifdef DEBUG
  std::ofstream of2("section_sfs2.g2");
  for (size_t kr=0; kr<reg_faces.size(); ++kr)
    {
      reg_faces[kr]->surface()->writeStandardHeader(of2);
      reg_faces[kr]->surface()->write(of2);
    }
#endif

  // Create spline surfaces
  shared_ptr<SurfaceModel> model3 =
    shared_ptr<SurfaceModel>(new SurfaceModel(eps, eps, neighbour,
					      kink, bend, reg_faces,
					      true));

  model3->replaceRegularSurfaces();  
#ifdef DEBUG
  int nmb_sfs0 = model3->nmbEntities();
  std::ofstream of3("section_sfs3.g2");
  for (int kh=0; kh<nmb_sfs0; ++kh)
    {
      shared_ptr<ParamSurface> surf = model3->getSurface(kh);
      surf->writeStandardHeader(of3);
      surf->write(of3);
    }
 #endif


  // Ensure consistent spline spaces across surface boundaries
  if (!model3->allSplines())
    return false;  // Not all surfaces are represented as spline surfaces
  model3->makeCommonSplineSpaces();

  // Rotate all surfaces to create volume blocks
  if (angle < 2.0*M_PI-eps)
    {
      // Rotate section back to the initial position
    }
  
  int nmb_sfs = model3->nmbEntities();
  SweepVolumeCreator createvol;
  vector<shared_ptr<ftVolume> > volumes;
  for (int ki=0; ki<nmb_sfs; ++ki)
    {
      shared_ptr<ParamSurface> surf = model3->getSurface(ki);
      shared_ptr<SplineSurface> sf = dynamic_pointer_cast<SplineSurface, ParamSurface>(surf);
      if (!sf.get())
	return false;  // Should not happen

      shared_ptr<ParamVolume> vol =
	shared_ptr<ParamVolume>(createvol.rotationalSweptVolume(*sf, angle, mid, axis));
      shared_ptr<ftVolume> vol2(new ftVolume(vol, -1));
      volumes.push_back(vol2);
    }
#ifdef DEBUG
  std::ofstream of4("volume_blocks.g2");
  for (size_t kr=0; kr<volumes.size(); ++kr)
    {
      volumes[kr]->getVolume()->writeStandardHeader(of4);
      volumes[kr]->getVolume()->write(of4);
    }
#endif

  // Create volume model
  volmodel = shared_ptr<VolumeModel>(new VolumeModel(volumes, eps, neighbour, kink, bend));
      
  return true;
}

//===========================================================================
bool VolumeModelCreator::linearSweptModel(shared_ptr<SurfaceModel>& sfmodel,
					  shared_ptr<VolumeModel>& volmodel)
//===========================================================================
{
  // This function recognizes a planar surface set swept along a line
  // in a boundary represented model and creates the associated
  // volume model. 
  // Non-planar surface sets in a similar configuration is not recognized and
  // neither is linear sweep where the curve along which to sweep is not linear

  // Check if the surface model represents a solid
  int nmb_bd = sfmodel->nmbBoundaries();
  if (nmb_bd > 0)
    return false;  // Not a solid

  // Create solid
  shared_ptr<Body> body(new Body(sfmodel));

  // Check if the surface model represents a linear swept object with the
  // restrictions explained above
  Point pnt, axis;
  double len;
  bool linear_swept = sfmodel->isLinearSwept(pnt, axis, len);
  if (!linear_swept)
    return false;

  // Fetch the faces lying in the specified plane
  vector<shared_ptr<ftSurface> > faces = sfmodel->facesInPlane(pnt, axis);
  
#ifdef DEBUG
  std::ofstream of1("section_sf.g2");
  for (size_t kr=0; kr<faces.size(); ++kr)
    {
      faces[kr]->surface()->writeStandardHeader(of1);
      faces[kr]->surface()->write(of1);
    }
#endif

  // Regularize face set
  double eps = sfmodel->getTolerances().gap;  // Tolerance
  double kink = sfmodel->getTolerances().kink;  // Angular tolerance
  double bend = sfmodel->getTolerances().bend;  // G1 tolerance
  double neighbour = sfmodel->getTolerances().neighbour;  // Adjacency tolerance
  RegularizeFaceSet reg(faces, eps, neighbour, kink, bend);
  shared_ptr<SurfaceModel> model2 = reg.getRegularModel();
  if (model2->nmbEntities() == 0)
    return false;

#ifdef DEBUG
  std::ofstream of2("section_sfs2.g2");
  int nmb = model2->nmbEntities();
  for (int kh=0; kh<nmb; ++kh)
    {
      shared_ptr<ParamSurface> sf = model2->getSurface(kh);
      sf->writeStandardHeader(of2);
      sf->write(of2);
    }
#endif

  // Create spline surfaces

  model2->replaceRegularSurfaces();  
#ifdef DEBUG
  int nmb_sfs0 = model2->nmbEntities();
  std::ofstream of3("section_sfs3.g2");
  for (int kh=0; kh<nmb_sfs0; ++kh)
    {
      shared_ptr<ParamSurface> surf = model2->getSurface(kh);
      surf->writeStandardHeader(of3);
      surf->write(of3);
    }
 #endif


  // Ensure consistent spline spaces across surface boundaries
  if (!model2->allSplines())
    return false;  // Not all surfaces are represented as spline surfaces
  model2->makeCommonSplineSpaces();

  // Create curve along which to sweep
  Point pnt2 = pnt + len*axis;
  shared_ptr<SplineCurve> curve(new SplineCurve(pnt, pnt2));

  // Perform sweep
  int nmb_sfs = model2->nmbEntities();
  SweepVolumeCreator createvol;
  vector<shared_ptr<ftVolume> > volumes;
  for (int ki=0; ki<nmb_sfs; ++ki)
    {
      shared_ptr<ParamSurface> surf = model2->getSurface(ki);
      shared_ptr<SplineSurface> sf = dynamic_pointer_cast<SplineSurface, ParamSurface>(surf);
      if (!sf.get())
	return false;  // Should not happen

      shared_ptr<ParamVolume> vol =
	shared_ptr<ParamVolume>(createvol.linearSweptVolume(*sf, *curve, pnt));
      shared_ptr<ftVolume> vol2(new ftVolume(vol, -1));
      volumes.push_back(vol2);
    }
#ifdef DEBUG
  std::ofstream of4("volume_blocks.g2");
  for (size_t kr=0; kr<volumes.size(); ++kr)
    {
      volumes[kr]->getVolume()->writeStandardHeader(of4);
      volumes[kr]->getVolume()->write(of4);
    }
#endif

  // Create volume model
  volmodel = shared_ptr<VolumeModel>(new VolumeModel(volumes, eps, neighbour, 
						     kink, bend));
      
  return true;
}

