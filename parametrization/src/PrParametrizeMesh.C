#include "GoTools/parametrization/PrParametrizeMesh.h"
#include "GoTools/parametrization/PrPrmShpPres.h"
#include "GoTools/parametrization/PrDijkstra.h"
#include <algorithm>

#ifdef _WIN32
#define M_PI 3.14159265358979
#endif

using std::cerr;
using std::cout;
using std::endl;

//using namespace std;

//-----------------------------------------------------------------------------
PrParametrizeMesh::PrParametrizeMesh() :
  convTolerance(1.0e-5)
//-----------------------------------------------------------------------------
{
  pb_ = new PrParametrizeBdy();
  pi_ = new PrPrmShpPres();
  //distance_finder_ = new Dijkstra();

  pb_->setParamKind(PrCHORDLENGTHBDY);

  pi_->setStartVectorKind(PrBARYCENTRE);
  pi_->setBiCGTolerance(convTolerance);
}


//-----------------------------------------------------------------------------
PrParametrizeMesh::~PrParametrizeMesh()
//-----------------------------------------------------------------------------
{
    delete pi_;
    delete pb_;
    //    delete distance_finder_;
}

//-----------------------------------------------------------------------------
PrParamTriangulation* PrParametrizeMesh::parametrize() 
//-----------------------------------------------------------------------------
{
  int i;
  PrParamTriangulation* paramtri = new PrParamTriangulation ();
  paramtri->attach(mesh_, basemesh_);
  
  std::shared_ptr<PrSubTriangulation> sub_tri(new PrSubTriangulation);
  vector<int> corners;
  for (i=0; i<basemesh_->findNumFaces(); i++) {

    cerr << "Processing coarse Triangle No. " << i << " - ";

    makeSubTriangulationFromTriangle(i, *sub_tri, corners);
    parametrizeSubTriangulation(sub_tri, corners);
    paramtri->makeCorrespondences(i, *sub_tri);
  }
  
  return paramtri;
}

//-----------------------------------------------------------------------------
void 
PrParametrizeMesh::
parametrizeSubTriangulation(std::shared_ptr<PrSubTriangulation> sub_tri,
			    vector<int>& corners)
//-----------------------------------------------------------------------------
{
//  cerr << " - parametrizeSubTriangulation: parameterizing the boundary ";

  pb_->attach(sub_tri);
  
  // set the parameter values for the corners
  int localInd0 = sub_tri->getLocalIndex(corners[0]);
  sub_tri->setU (localInd0, 0.0);
  sub_tri->setV (localInd0, 0.0);

  int localInd1 = sub_tri->getLocalIndex(corners[1]);
  sub_tri->setU (localInd1, 1.0);
  sub_tri->setV (localInd1, 0.0);

  int localInd2 = sub_tri->getLocalIndex(corners[2]);
  sub_tri->setU (localInd2, 0.0);
  sub_tri->setV (localInd2, 1.0);

  // parameterize the boundary
  pb_->parametrizeSide(localInd0, localInd1);
  pb_->parametrizeSide(localInd1, localInd2);
  pb_->parametrizeSide(localInd2, localInd0);

//  cerr << "\n     ";
//  for (int i=0; i<sub_tri->findNumBdyNodes(); i++)
//    cerr << i << " -> (" << sub_tri->getU(i) << "," << sub_tri->getV(i) << "), ";
//  cerr << "\n - parametrizeSubTriangulation: parameterizing the interior ";

  // parameterize the interior
  pi_->attach(sub_tri);
  pi_->parametrize();

//  cerr << "\n     ";
//  for (i=sub_tri->findNumBdyNodes(); i<sub_tri->getNumNodes(); i++)
//    cerr << i << " -> (" << sub_tri->getU(i) << "," << sub_tri->getV(i) << "), ";
//  cerr << endl;
}


//-----------------------------------------------------------------------------
void 
PrParametrizeMesh::
makeSubTriangulationFromTriangle(int i, 
				 PrSubTriangulation& sub_tri,
				 vector<int>& corners) 
//-----------------------------------------------------------------------------
{
  corners.resize(3);

  const PrTriangle& tri = basemesh_->getPrTriangle(i);
  corners[0] = tri.n1();
  corners[1] = tri.n2();
  corners[2] = tri.n3();

//  cerr << " - makeSubTriangulationFromTriangle found corner node indices: ";
//  cerr << corners[0] << ", " << corners[1] << ", " << corners[2] << endl;

  makeSubTriangulation(*mesh_, corners, sub_tri);
}


//-----------------------------------------------------------------------------
bool 
PrParametrizeMesh::
makePolygon( PrTriangulation_OP& triangulation, 
	     const vector<int>&  nodes, 
	     vector<int>&        polygon) 
//-----------------------------------------------------------------------------
{
  // given a "triangulation" and a (ordered) set of "nodes",
  // this routine returns a "polygon", consisting of those
  // nodes (ordered) who form the shortest paths between
  // the given nodes in the triangulation
  //
  // i.e.: given a triangle (or arbitrary polygon) in the
  // coarse mesh, this algorithm finds the corresponding
  // triangle (polygon) in the original (fine) mesh

  polygon.empty();	

  if (nodes.size()<2)
    return false;

  vector<int> path;
  for (size_t i=0; i<nodes.size(); i++) {
    if(makePath(triangulation, nodes[i], nodes[(i+1)%nodes.size()], path)) {
      for (size_t j=0; j<path.size()-1; j++) {
        polygon.push_back(path[j]);
      }
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
bool 
PrParametrizeMesh::
makePath( PrTriangulation_OP& triangulation, 
	  int                 n1, 
	  int                 n2,
	  vector<int>&   path) 
//-----------------------------------------------------------------------------
{
  bool reverse = false;

  if (n2 < n1) {
    // always return the path from the node with the smaller index
    // to the one with the larger index!
    reverse = true;
    int c = n1;
    n1 = n2;
    n2 = c;
  }

  // Find the shortest path in triangulation connecting n1 and n2
  path.clear();	
  
  Dijkstra dist_finder;
  dist_finder.setGraph(&triangulation);
  dist_finder.initialize();
  dist_finder.setSource(n2);
  dist_finder.run(n1);

//   distance_finder_->setGraph(&triangulation);
//   distance_finder_->initialize();
//   distance_finder_->setSource(n2);
//   distance_finder_->run(n1);

  int curr = n1;
  path.push_back(curr);
  while (curr != n2) {
    //curr = distance_finder_->closestNeighbour(curr);
    curr = dist_finder.closestNeighbour(curr);  
    path.push_back(curr);
  }

  if (reverse) {
    // revert the order of entries in "path"
      int ps = (int)path.size();
    for (int i=0; i<ps/2; i++) {
      int c = path[i];
      path[i] = path[ps-i-1];
      path[ps-i-1] = c;
    }
  }

//  cerr << " - - makePath found the following path between " << n1 << " and ";
//  cerr << n2 << ":\n      ";
//  for (int i=0; i<path.size(); i++)
//    cerr << path[i] << ", ";
//  cerr << endl;

  return true;
}



//-----------------------------------------------------------------------------
bool 
PrParametrizeMesh::
makeSubTriangulation ( PrTriangulation_OP&  triangulation, 
		       const vector<int>&   nodes, 
		       PrSubTriangulation&  subtriangulation) 
//-----------------------------------------------------------------------------
{
  // Make a subtriangulation of triang along the polygon induced by the 
  // shortest paths between the nodes in nodes
  // Use the same nodeset ?

  vector<int> polygon;

  if (!makePolygon(triangulation, nodes, polygon))
    return false; 

  std::vector<int> unique_vect(polygon);
  std::sort(unique_vect.begin(), unique_vect.end());
  if (std::unique(unique_vect.begin(), unique_vect.end())!=unique_vect.end())
  {
    std::cerr << "polygon has more than one occourence of a vertex\n";
  }
//  cerr << " - makeSubTriangulation found the following path in the fine mesh: ";
//  for (int i=0; i<polygon.size(); i++)
//    cerr << polygon[i] << ", ";
//  cerr << endl;
  
  if (!makeSubTriangulationInsidePolygon(triangulation,
					 polygon,
					 subtriangulation))
    return false;
  
//  cerr << " - makeSubTriangulation: the SubTriangulation has been created!";
//  cerr << endl;

  return true;
}

//-----------------------------------------------------------------------------
bool 
PrParametrizeMesh::
makeSubTriangulationInsidePolygon (PrTriangulation_OP& triangulation, 
				   const vector<int>&  polygon, 
				   PrSubTriangulation& subtriangulation) 
//-----------------------------------------------------------------------------
{
  // Make a subtriangulation of triang inside polygon. Assume polygon to be a 
  // closed ring of neighbouring nodes

  set<int> int_set;
  getConnectedNodes(triangulation, polygon, int_set);

#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
  vector<int> interior;
  for (set<int>::const_iterator it = int_set.begin(); it != int_set.end(); ++it)
    interior.push_back(*it);
#else
  vector<int> interior(int_set.begin(), int_set.end());
#endif

//  cerr << " - makeSubTriangulationInsidePolygon found interior nodes:\n     ";
//  for (int i=0; i<interior.size(); i++)
//    cerr << interior[i] << ", ";
//  cerr << endl;

  subtriangulation.attach(&triangulation);
  subtriangulation.initialize(polygon, interior);

  return true;
}


//-----------------------------------------------------------------------------
bool 
PrParametrizeMesh::
getConnectedNodes( PrTriangulation_OP& triangulation, 
		   const vector<int>&  boundary, 
		   set<int>&           interior)
//-----------------------------------------------------------------------------
{
  // Assume the vertices in polygon is a closed oriented polygon on
  // triangulation. Returns the indices of the vertices interior to
  // polygon. If the orientation is reversed it returns the exterior
  // If the polygon is open and does not split triangulation all 
  // nodes will be found
  
  if (boundary.size()==0)
    return false;
  
  int curr, next;
  set<int> bdy_nodes(boundary.begin(), boundary.end());

  interior.clear();

  // set seed
  std::queue<int> queue;
  queue.push(boundary[0]);
  interior.insert(boundary[0]);

  // March out without crossing polygon
  vector<int> neighbours;

  while (!queue.empty()) {
    curr = queue.front(); queue.pop();
    triangulation.getNeighbours(curr, neighbours);

    if (bdy_nodes.find(curr) != bdy_nodes.end()) {
      // use only "interior" neighbours for "boundary" nodes:
      vector<int> new_nbrs;
      getInteriorNeighbours(curr, neighbours, boundary, new_nbrs);

      // transfer "new_nbrs" to "neighbours"
      neighbours.assign(new_nbrs.begin(), new_nbrs.end());
    }

    for (size_t i=0; i<neighbours.size(); i++) {
      next = neighbours[i];
      if (interior.find(next) == interior.end()) {
        interior.insert(next);
        queue.push(next);
      }
    }
  }

  // erase "boundary" nodes from the "interior" set
  for (size_t i=0; i<boundary.size(); i++)
    interior.erase(boundary[i]);
		   
  return (true);
}


//-----------------------------------------------------------------------------
int
PrParametrizeMesh::
getLeftNode (PrTriangulation_OP& triangulation,
	     int n1,
	     int n2)
//-----------------------------------------------------------------------------
{
  // find the third node of the triangle with nodes n1, n2 in
  // anti-clockwise order (i.e. the triangle to the "left")
  
  // abbreviations:

  #define Node triangulation.getPrNode
  #define Triangle triangulation.getPrTriangle

  int tr = Node(n1).tr();

  while (Triangle(tr).getAnticlockwiseNode(n1) != n2) {
    tr = Triangle(tr).getLeftTriangle(n1);
  }
  
  return (Triangle(tr).getClockwiseNode(n1));
}


//-----------------------------------------------------------------------------
void 
PrParametrizeMesh::
getInteriorNeighbours( int                v,
		       const vector<int>& neighbours,
		       const vector<int>& boundary,
		       vector<int>&       new_nbrs )
//-----------------------------------------------------------------------------
{
    int ns = (int)neighbours.size();
  int start = 0;
  while (!edgeInPolygon(v, neighbours[start], boundary))
    start++;
  bool inside = true;
  int next = start;
  do {
    if (edgeInPolygon(neighbours[next], v, boundary))
      inside = false;
    next = (next+1) % ns;
    if (edgeInPolygon(v, neighbours[next], boundary))
      inside = true;
    if (inside)
      new_nbrs.push_back(neighbours[next]);
  } while (next != start);
}

//-----------------------------------------------------------------------------
bool
PrParametrizeMesh::
edgeInPolygon( int n1, int n2, const vector<int>& polygon )
//-----------------------------------------------------------------------------
{
    int ps = (int)polygon.size();
  for (int i=0; i<ps-1; i++)
    if ((polygon[i] == n1) && (polygon[i+1] == n2))
      return true;
  if ((polygon[ps-1] == n1) && (polygon[0] == n2))
    return true;
  return false;
}



//-----------------------------------------------------------------------------
int
PrParametrizeMesh::
castRay( PrTriangulation_OP& triangulation, 
	 int vs, double angle, double dMax, const set<int>& T, 
	 vector<int>& t_path, vector<Vector3D>& v_path)
//-----------------------------------------------------------------------------
//
// casts a ray from the vertex "vs" into the direction of 
// "angle". angle = 0 is the direction to the ccw node in
// vs's leading triangle (remark: that is also the first
// node in the result of the "getNeighbours" routine).
// The ray is cast as long as its length does not exceed "dMax" 
// and no triangle in "T" is visited. The triangles visited are 
// stored in "t_path". If the path does not hit a triangle in "T", 
// the return value is -1, otherwise it is the triangle index
// "v_path" contains the vertices of the path
//
{
  double bs;
  int ts;
  convertAngle(triangulation, angle, vs, bs, ts);

  // initialize variables
  PrTriangle tri_curr = triangulation.getPrTriangle(ts);

  int v1 = vs;
  int v2 = tri_curr.getAnticlockwiseNode(v1);
  int v3 = tri_curr.getClockwiseNode(v1);

  int tn = tri_curr.getOppositeTriangle(v1);
  PrTriangle tri_next = triangulation.getPrTriangle(tn);

  int v4 = tri_next.getAnticlockwiseNode(v2);

  Vector3D p1 = triangulation.getPrNode(v1).point();
  Vector3D p2 = triangulation.getPrNode(v2).point();
  Vector3D p3 = triangulation.getPrNode(v3).point();
  Vector3D p4 = triangulation.getPrNode(v4).point();

  double a = 1;
  double b = bs;
  double c = 0.0;
  bool f = true;

  Vector3D pb = b*p2 + (1-b)*p3;
  double d = pb.dist(p1);

  // initialize paths
  t_path.clear();
  t_path.push_back(ts);
  
  v_path.clear();
  v_path.push_back(p1);
  v_path.push_back(pb);

  // if ts already is in T -> return 
  if (T.count(ts))
    return (ts);
  
  t_path.push_back(tn);

  // loop ray casting until maximal length is reached or triangle in T is hit
  while (!T.count(tn) && (d < dMax)) {

    // get next vertex
    f = getNextV(p1, p2, p3, p4, a, b, f, c);

    // update variables to go on
    if (f) {
      ts = tn; 
      tri_curr = tri_next;
      tn = tri_curr.getOppositeTriangle(v3);
      tri_next = triangulation.getPrTriangle(tn);

      v1 = v3; v3 = v4; v4 = tri_next.getAnticlockwiseNode(v2);
      p1 = p3; p3 = p4; p4 = triangulation.getPrNode(v4).point();

      a = 1-b; b = 1-c;
    } 
    else {
      ts = tn; 
      tri_curr = tri_next;
      tn = tri_curr.getOppositeTriangle(v2);
      tri_next = triangulation.getPrTriangle(tn);

      v1 = v2; v2 = v4; v4 = tri_next.getClockwiseNode(v3);
      p1 = p2; p2 = p4; p4 = triangulation.getPrNode(v4).point();

      a = b; b = c;
    }

    // update total path length
    Vector3D pc = b*p2 + (1-b)*p3;
    d += pc.dist(pb);
    pb.x() = pc.x(); pb.y() = pc.y(); pb.z() = pc.z(); 

    t_path.push_back(tn);
    v_path.push_back(pb);
  }
	       
  if (T.count(tn))
    return (tn);
  else
    return (-1);
}

//-----------------------------------------------------------------------------
bool
PrParametrizeMesh::
getNextV( const Vector3D& p1, const Vector3D& p2, 
	  const Vector3D& p3, const Vector3D& p4, 
	  double a, double b, bool f, double& c)
//-----------------------------------------------------------------------------

{
  // determine barycentric coordinates (r,s,t) of p4 w.r.t (p1,p2,p3)
  // after flattening

  Vector3D u = p1-p3;
  Vector3D v = p2-p3;
  Vector3D w = p4-p3;

  double cPvw = (v.cross(w)).length();
  double cPuv = (u.cross(v)).length();
  double dPvw = v*w;
  double dPuv = u*v;

  double r = -cPvw/cPuv;
  double s = (cPvw*dPuv + cPuv*dPvw)/v.length2()/cPuv;
  double t = 1-r-s;

  double q1,q2;
  if (f) {
    q1 = t; q2 = 1-b; 
  }
  else {
    q1 = 1-s; q2 = -b; 
  }

  double mu = a*(1-b)/(a*q1+r*q2);
  double nu = a*b/(a*(1-q1)-r*q2);

  #define epsi 0.000001
  if (fabs(mu) < epsi)
    mu = epsi;
  if (fabs(1.0-mu) < epsi)
    mu = 1.0-epsi;
  if (fabs(nu) < epsi)
    nu = epsi;
  if (fabs(1.0-nu) < epsi)
    nu = 1.0-epsi;

  if ((mu > 0.0) && (mu < 1.0)) {
    if ((nu > 0.0) && (nu < 1.0)) {
      cerr << "!!_OOPS_!! getNextV : mu and nu are BOTH within [0,1]: mu = ";
      cerr << mu << " / nu = " << nu << endl;

      if (fabs(0.5-mu) < fabs(0.5-nu)) {
	c = mu;
	return(true);
      }
      else {
	c = nu;
	return (false);
      }
    }

    c = mu;
    return (true);
  }

  if ((nu > 0.0) && (nu < 1.0)) {
    c = nu;
    return (false);
  }

  cerr << "!!_OOPS_!! getNextV : mu and nu are NOT within [0,1]: mu = ";
  cerr << mu << " / nu = " << nu << endl;
  return (false);
}


//-----------------------------------------------------------------------------
double
PrParametrizeMesh::
getIsoline(const PrTriangulation_OP& triangulation, int vs, int vd, 
	    vector<int>& t_path, vector<Vector3D>& v_path)
//-----------------------------------------------------------------------------
//
// runs the Dijkstra algorithm from source point "vs" to
// destination point "vd" and determines the iso-distance-line
// of vd.
// The vertices of that isoline are stored in "v_path", whereas
// the indices of the triangles visited by that line are
// stored in "t_path".
// Return value is the distance from "vs" to "vd".
//
{
  Dijkstra dist_finder;

  #define Node triangulation.getPrNode
  #define Triangle triangulation.getPrTriangle
  #define Distance dist_finder.getDistance
  //#define Distance distance_finder_->getDistance

  t_path.clear();
  v_path.clear();

  // run Dijkstra's algorithm
  dist_finder.setGraph(&triangulation);
  dist_finder.initialize();
  dist_finder.setSource(vs);
  dist_finder.run(vd);
//   distance_finder_->setGraph(&triangulation);
//   distance_finder_->initialize();
//   distance_finder_->setSource(vs);
//   distance_finder_->run(vd);

  double dist = Distance(vd);

  // traverse the isoline of "vd" ccw (seen from vs)
  v_path.push_back(Node(vd).point());

  vector<int> nbrs;
  triangulation.getNeighbours(vd,nbrs);
//  cout << "vd's distance: " << dist << endl;

  // find first triangle around vd
  int t_first = Node(vd).tr();
  size_t cnt = 0;
  while (((Distance(Triangle(t_first).getClockwiseNode(vd)) > dist) ||
	  (Distance(Triangle(t_first).getAnticlockwiseNode(vd)) < dist)) &&
	 (cnt <= nbrs.size())) {
    t_first = Triangle(t_first).getLeftTriangle(vd);
    cnt++;
  }
  if (cnt > nbrs.size()) {
    cerr << "could not find leading triangle in isoline determination!" <<endl;
    return (-1);
  }

  // find last triangle around vd
  int t_last = Node(vd).tr();
  cnt = 0;
  while (((Distance(Triangle(t_last).getAnticlockwiseNode(vd)) > dist) ||
	  (Distance(Triangle(t_last).getClockwiseNode(vd)) < dist)) &&
	 (cnt <= nbrs.size())) {
    t_last = Triangle(t_last).getLeftTriangle(vd);
    cnt++;
  }
  if (cnt > nbrs.size()) {
    cerr << "could not find closing triangle in isoline determination!" <<endl;
    return (-1);
  }

  // traverse isoline from "t_first" to "t_last"
  int direct = vd;
  int pMinus, pPlus, pNext;
  double dMinus, dPlus, dNext;
  int t_next;
  Vector3D vi;
  double w;

  while (t_first != t_last) {
    t_next = Triangle(t_first).getOppositeTriangle(direct);

    pMinus = Triangle(t_first).getClockwiseNode(direct);
    pPlus  = Triangle(t_first).getAnticlockwiseNode(direct);
    pNext  = Triangle(t_next).getClockwiseNode(pMinus);

    dMinus = Distance(pMinus);
    dPlus  = Distance(pPlus);
    dNext  = Distance(pNext);

    if ((dMinus > dist) || (dPlus < dist))
      cerr << "!_OOPS_! getIsoline: something strange happened!!!" << endl;

    // interpolate pMinus and pPlus:
    w = (dist-dMinus)/(dPlus-dMinus);
    vi = (1-w)*(Node(pMinus).point()) + w*(Node(pPlus).point());

    // push back values in the triangle and vertex list
    t_path.push_back(t_first);
    v_path.push_back(vi);
  
//    if (dNext == dist)
//      cerr << "here's trouble!!!" << endl;

    // carry on, carry on...
    t_first = t_next;
    if (dNext > dist)
      direct = pPlus;
    else
      direct = pMinus;
  }

  t_path.push_back(t_first);
  v_path.push_back(Node(vd).point());

  return (dist);
}

//-----------------------------------------------------------------------------
void
PrParametrizeMesh::
convertAngle(const PrTriangulation_OP& triangulation, 
	      double theta, int v, double& bc, int& t)
//-----------------------------------------------------------------------------
//
// converts the direction indicated by "angle" around the node
// "v" into a (double,int) pair. The "int" value is the index of
// the triangle around "v" in direction "angle", the "double"
// value is the barycentric coordinate of the intersection
// point of a ray in direction "angle" with the edge opposite
// to "v" in the triangle w.r.t to the vertices of that edge.
//
{
  #define Node triangulation.getPrNode
  #define Triangle triangulation.getPrTriangle

  size_t i;
  static int recentV = -1;
  static vector<int> triangles;
  static vector<double> angles;

  if (v != recentV) {
    // update "triangles" information
    int tr1 = Node(v).tr();
    triangles.clear();
    triangles.push_back(tr1);
    int tr = Triangle(tr1).getLeftTriangle(v);
    while ((tr > -1) && (tr != tr1)) {
      triangles.push_back(tr);
      tr = Triangle(tr).getLeftTriangle(v);
    }

    // update "angles" information
    Vector3D v0 = Node(v).point();
    double alpha = 0.0;
    angles.clear();
    angles.push_back(0.0);
    for (i=0; i<triangles.size(); i++) {
      Vector3D v1 = Node(Triangle(triangles[i]).getAnticlockwiseNode(v)).point()-v0;
      Vector3D v2 = Node(Triangle(triangles[i]).getClockwiseNode(v)).point()-v0;
      alpha += v1.angle(v2);
      angles.push_back(alpha);
    }
#ifdef MICROSOFT
    const double PI = 3.14159265358979323846264338327950288419716939937510;
    double scaleFactor = 2.0*PI/alpha;
#else
    double scaleFactor = 2.0*M_PI/alpha;
#endif
    for (i=0; i<angles.size(); i++)
      angles[i] *= scaleFactor;

    recentV = v;
  }

  // find triangle index
  int ti = 0;
  while (angles[ti] < theta)
    ti++;
  t = triangles[ti-1];

  // linearly interpolate to get barycentric coordinate
  // (of course, this is only approximately correct,
  // but good enough for our purposes!)
  bc = (angles[ti]-theta) / (angles[ti]-angles[ti-1]);

  #define epsi 0.000001
  if (fabs(bc) < epsi)
    bc = epsi;
  if (fabs(1.0-bc) < epsi)
    bc = 1.0-epsi;
}

