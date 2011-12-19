/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

/********************************************************************
 FILENAME    : PrParametrizeInt.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : March 97
 DESCRIPTION : Implementation of methods in the class PrParametrizeInt.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrParametrizeInt.h"


#include "GoTools/parametrization/PrBiCGStab.h"
#include "GoTools/parametrization/PrMatSparse.h"
#include "GoTools/parametrization/PrVec.h"

#include <fstream>

using namespace std;

// PRIVATE METHODS

//-----------------------------------------------------------------------------
int PrParametrizeInt::getNumIntNghrs(int i)
//-----------------------------------------------------------------------------
// How many interior neighbours does the i-th node have?
{
  g_->getNeighbours(i,neighbours_);

  int nInt = 0;
  for(size_t j=0; j< neighbours_.size(); j++)
  {
    if(!g_->isBoundary(neighbours_[j])) nInt++;
  }
  return nInt;
}


//-----------------------------------------------------------------------------
int PrParametrizeInt::getNumInt2Nghrs(int i, vector<int>& fixedPnts)
//-----------------------------------------------------------------------------
// How many interior neighbours in the 2-neighbourhood does the i-th node have?
// without counting fixed neighbours
{
  size_t j;
  g_->getNeighbours(i,neighbours_);

  vector<int> neighbours1;
  vector<int> neighbours2;

  for (j=0; j<neighbours_.size(); j++)
    if(!isFixed(neighbours_[j], fixedPnts))
      neighbours1.push_back(neighbours_[j]);

  for (j=0; j<neighbours_.size(); j++)
  {
    g_->getNeighbours(neighbours_[j], neighbours2);
    
    for (size_t k=0; k<neighbours2.size(); k++)
    {
      bool newVertex = true;

      if ((neighbours2[k] == i) || (isFixed(neighbours2[k], fixedPnts)))
	newVertex = false;
      else
	for (size_t l=0; l<neighbours1.size(); l++)
	  if (neighbours1[l] == neighbours2[k])
	  {
	    newVertex = false;
	    break;
	  }

      if(newVertex)
	neighbours1.push_back(neighbours2[k]);
    }
  }

  return (int)neighbours1.size();
}


//-----------------------------------------------------------------------------
void PrParametrizeInt::findBarycentre(double& ucentre, double& vcentre)
//-----------------------------------------------------------------------------
//   Compute barycentre of the parameter points of the boundary
//   of the planar graph g_.
{
  ucentre = 0.0;
  vcentre = 0.0;
  int numBdyNodes = 0;
  for(int i=0; i < g_->getNumNodes(); i++)
  {
    if(g_->isBoundary(i))
    {
      ucentre += g_->getU(i);
      vcentre += g_->getV(i);
      numBdyNodes++;
    }
  }

  double nbdouble = (double)numBdyNodes;
  ucentre /= nbdouble;
  vcentre /= nbdouble;
}


// PUBLIC METHODS

//-----------------------------------------------------------------------------
PrParametrizeInt::PrParametrizeInt()
//-----------------------------------------------------------------------------
{
  tolerance_ = 1.0e-6;
  startvectortype_ = PrBARYCENTRE;
}
//-----------------------------------------------------------------------------
PrParametrizeInt::~PrParametrizeInt()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
void PrParametrizeInt::attach(shared_ptr<PrOrganizedPoints> graph)
//-----------------------------------------------------------------------------
{
  g_ = graph;
}

//-----------------------------------------------------------------------------
bool PrParametrizeInt::parametrize()
//-----------------------------------------------------------------------------
//   Parametrize the interior nodes of the planar graph g_.
//   The boundary parameter points are assumed to be already set.
//   M.F. Dec. 98.
{
  // To find the internal (u,v) nodes, we solve a linear system
  // twice: once for u and once for v. The matrix A is the same
  // in both cases. Only the right-hand side changes.

  // There is one equation for each interior node,
  // so we solve Au=b1 and Av=b2 where A is ni * ni and ni is the
  // number of interior nodes.

  int n = g_->getNumNodes();
  int ni = n - g_->findNumBdyNodes();
  if(ni == 0) return true;

  // find number of non-zeros in the matrix A
  int i;
  int numNonZeros = 0;
  for(i=0; i<n; i++)
  {
    if(!g_->isBoundary(i))
       numNonZeros += (getNumIntNghrs(i) + 1);
  }

  PrMatSparse A (ni,ni,numNonZeros);

  PrVec b1(ni, 0.0);
  PrVec b2(ni, 0.0);
  PrVec uvec(ni);
  PrVec vvec(ni);

  // Make a permutation array for mapping global i to interior i.
  vector<int> permute(n);
  int j = 0;
  for(i=0; i<n; i++)
  {
    if(!g_->isBoundary(i))
    {
      permute[i] = j;
      j++;
    }
  }
  
  //debug

  //  s_o << "Here's permute...\n";
  //  for(i=1; i<=n; i++) s_o << permute(i) << endl;

   
  // Set start vector.

  if(startvectortype_ == PrBARYCENTRE)
  {
    // Set start vector to barycentre of boundary.
    double ucentre,vcentre;
    findBarycentre(ucentre,vcentre);
//     std::cout << "using barycentre" << std::endl;
//     std::cout << ucentre << " " << vcentre << std::endl;
    for(i = 0; i < ni; i++)
    {
      uvec(i) = ucentre;
      vvec(i) = vcentre;
    }
  }
  else if(startvectortype_ == PrFROMUV)
  {
    for(i=0; i<n; i++)
    {
      if(!g_->isBoundary(i))
      {
        uvec(permute[i]) = g_->getU(i);
        vvec(permute[i]) = g_->getV(i);
      }
    }
  }
  else return false;

  // DEBUG
//    for (i=0; i<n; i++)
//      {
//        if(g_->isBoundary(i))
//        {
//  	std::cout << "parameter value: " << g_->getU(i) << " " << g_->getV(i) << std::endl;
//        }
//      }
//    int i2;
//    vector<int> neighbours;
//    for (i=0; i<n; i++)
//      if (g_->isBoundary(i))
//        break;
//    std::cout << "parameter value: " << g_->getU(i) << " " << g_->getV(i) << std::endl;
//    g_->getNeighbours(i,neighbours); 
//    i2 = neighbours[0];
//    while (i2 != i)
//      {
//        std::cout << "parameter value: " << g_->getU(i2) << " " << g_->getV(i2) << std::endl; 
//        g_->getNeighbours(i2,neighbours); 
//        i2 = neighbours[0];
//      }
      

  // END DEBUG


  // Initialize right hand side to zero
  // @afr: Already done in construction of b1, b2.

  int offset = 0;

  for(i=0; i<n; i++)
  {
    if(!g_->isBoundary(i))
    {
      A.irow(permute[i]) = offset;
      A(offset) = 1.0;
      A.jcol(offset) = permute[i];
      offset++;

      // @afr: Seems this is set for the benefit of makeWeights() below.
      g_->getNeighbours(i,neighbours_);
      int degree = (int)neighbours_.size();
      makeWeights(i); // Ignoring return value.
         // set the weights_ array, indices 0,1,2,...,degree-1.

      for(j=0; j< degree; j++)
      {
        int k = neighbours_[j];
        if(g_->isBoundary(k))
        {
          b1(permute[i]) += weights_[j] * g_->getU(k);
          b2(permute[i]) += weights_[j] * g_->getV(k);
        }
        else
        {
          A(offset) = -weights_[j];
          //debug
          //s_o << oform("%2.20f\n",(*A)(offset));

          A.jcol(offset) = permute[k];
          offset++;
        }
      }
    }
  }

// USEFUL DEBUG!

//    A.print(std::cout);
//    b1.print(std::cout);
//    b2.print(std::cout);
//    uvec.print(std::cout);
//    vvec.print(std::cout);
//    /*for(i=1; i<=ni; i++) s_o << oform("%2.20f\n",(*b1)(i));
//    for(i=1; i<=ni; i++) s_o << oform("%2.20f\n",(*b2)(i)); */

// END OF USEFUL DEBUG

  PrBiCGStab solver;
  solver.setMaxIterations(ni);
  solver.setTolerance(tolerance_);
  solver.solve(A,uvec,b1);
//   std::cout << "Converge " << solver.converged() << std::endl;

#ifdef PRDEBUG
  double cpu_time = solver.getCPUTime();
  int noIts = solver.getItCount();
  bool converged = solver.converged();
  std::cout << "unknowns = " << ni << "  cpu_time = " << cpu_time
       << "  no_its = " << noIts << "  converged = " << converged << std::endl;
#endif
// END OF DEBUG

  solver.solve(A,vvec,b2);

#ifdef PRDEBUG
  cpu_time = solver.getCPUTime();
  noIts = solver.getItCount();
  converged = solver.converged();
  std::cout << "unknowns = " << ni << "  cpu_time = " << cpu_time
       << "  no_its = " << noIts << "  converged = " << converged << std::endl;
#endif
// END OF DEBUG

  //uvec->print(s_o,"solution1");
  //vvec->print(s_o,"solution2");
  // Copy the results to the u and v arrays.
  for(i=0; i<n; i++)
  {
    if(!g_->isBoundary(i))
    {
      g_->setU(i, uvec(permute[i]));
      g_->setV(i, vvec(permute[i]));
    }
  }
  // DEBUG
//    for (i=0; i<n; i++)
//      {
//        std::cout << i << ", parameter value: " << g_->getU(i) << " " << g_->getV(i) << std::endl;
//      }
  // END DEBUG


  return true;
}

//-----------------------------------------------------------------------------
bool PrParametrizeInt::isFixed(int k, vector<int>& fixedPnts)
//-----------------------------------------------------------------------------
{
  for(size_t i=0; i<fixedPnts.size(); i++)
    {
      if(k == fixedPnts[i]) return true;
    }
 return false;
}

//-----------------------------------------------------------------------------
bool PrParametrizeInt::parametrize3d(vector<int>& fixedPnts, 
				     vector<double>& uvw)
//-----------------------------------------------------------------------------
//   Parametrize the nodes of the 3D grapg g_ except those
//   with indices in "fixedPnts". The parameterization is done
//   in 3D and the result is returned as vector "uvw"
{
  // To find the internal (u,v,w) nodes, we solve a linear system
  // thrice: once for u, v and w. The matrix A is the same
  // in both cases. Only the right-hand side changes.

  // There is one equation for each interior node,
  // so we solve Au=b1, Av=b2 and Av=b3 where A is ni * ni and ni is the
  // number of interior nodes.

  int n = g_->getNumNodes();
  int ni = n - (int)fixedPnts.size();
  if(ni == 0) return true;

  //
  // setup the "full" matrix A (i.e. without boundary vertices)
  //
  //  remark: this will only work for closed meshes!!!
  //

  // find number of non-zeros in the matrix A
  int i;
  int numNonZeros = 0;
  for(i=0; i<n; i++)
  {
    g_->getNeighbours(i, neighbours_);
    numNonZeros += ((int)neighbours_.size() + 1);
  }

  PrMatSparse A (n, n, numNonZeros);
  int offset = 0;

  for(i=0; i<n; i++)
  {
    A.irow(i) = offset;
    A(offset) = 1.0;
    A.jcol(offset) = i;
    offset++;

    g_->getNeighbours(i,neighbours_);
    int degree = (int)neighbours_.size();
    makeWeights(i); // Ignoring return value.
    
    for(int j=0; j< degree; j++)
    {
      int k = neighbours_[j];

      A(offset) = weights_[j];

      A.jcol(offset) = k;
      offset++;
    }
  }

  cerr << "Matrix A upset!" << endl;
//  A.print(cerr);

  //
  // setup the matrix B
  //

  // setup the linear system

  // find number of non-zeros in the matrix B
  numNonZeros = 0;
  for(i=0; i<n; i++)
  {
    if(!isFixed(i,fixedPnts))
       numNonZeros += (getNumInt2Nghrs(i, fixedPnts) + 1);
  }

  PrMatSparse B (ni,ni,numNonZeros);

  PrVec b1(ni, 0.0);
  PrVec b2(ni, 0.0);
  PrVec b3(ni, 0.0);
  PrVec uvec(ni);
  PrVec vvec(ni);
  PrVec wvec(ni);

  // Make a permutation array for mapping global i to interior i.
  vector<int> permute(n);
  int j = 0;
  for(i=0; i<n; i++)
  {
    if(!isFixed(i,fixedPnts))
    {
      permute[i] = j;
      j++;
    }
  }
     
  // Set start vector to the original xyz-coordinates.
  for(i = 0; i < ni; i++)
  {
    Vector3D p=g_->get3dNode(i);
    uvec(i) = p.x();
    vvec(i) = p.y();
    wvec(i) = p.z();
  }

  // Initialize right hand side to zero
  // @afr: Already done in construction of b1, b2, b3.

  offset = 0;

  for(i=0; i<n; i++)
  {
    if(!isFixed(i,fixedPnts))
    {
      g_->getNeighbours(i,neighbours_);

      // set the diagonal entry

      B.irow(permute[i]) = offset;
      B(offset) = 1.0;

      for (size_t k=0; k<neighbours_.size(); k++)
      {
	int kk = neighbours_[k];
	// A^t * A
	B(offset) += A(kk,i) * A(kk,i);
	// A^2
	// B(offset) += A(kk,i) * A(i, kk);
      }
      
      B.jcol(offset) = permute[i];
      offset++;

      // handle neighbours in 1-neighbourhood

      int degree = (int)neighbours_.size();

      for(j=0; j<degree; j++)
      {
        int jj = neighbours_[j];

	// compute the b_ij

	double b_ij; 

	// A^t * A
	b_ij = - A(i, jj) - A(jj, i);
	// A^2
	// b_ij = - 2*A(i, jj);

	vector<int> commonNeighbours;
	g_->getCommonNeighbours(i, jj, commonNeighbours);

	for (size_t k=0; k<commonNeighbours.size(); k++)
	{
	  int kk = commonNeighbours[k];
	  // A^t * A
	  b_ij += A(kk, i) * A(kk,jj);
	  // A^2
	  // b_ij += A(i, kk) * A(kk,jj);
	}

        if(isFixed(jj,fixedPnts))
        {
          Vector3D p=g_->get3dNode(jj);
          b1(permute[i]) -= b_ij * p.x();
          b2(permute[i]) -= b_ij * p.y();
          b3(permute[i]) -= b_ij * p.z();
        }
        else
        {
	  B(offset) = b_ij;
          B.jcol(offset) = permute[jj];
          offset++;
        }
      }

      // handle neighbours in 2-neighbourhood

      g_->get2Neighbours(i,neighbours_);
      degree = (int)neighbours_.size();

      for(j=0; j<degree; j++)
      {
        int jj = neighbours_[j];

	// compute the b_ij

	double b_ij = 0.0; 

	vector<int> commonNeighbours;
	g_->getCommonNeighbours(i, jj, commonNeighbours);

	for (size_t k=0; k<commonNeighbours.size(); k++)
	{
	  int kk = commonNeighbours[k];
	  // A^t * A
	  b_ij += A(kk, i) * A(kk,jj);
	  // A^2
	  // b_ij += A(i, kk) * A(kk,jj);
	}

        if(isFixed(jj,fixedPnts))
        {
          Vector3D p=g_->get3dNode(jj);
          b1(permute[i]) -= b_ij * p.x();
          b2(permute[i]) -= b_ij * p.y();
          b3(permute[i]) -= b_ij * p.z();
        }
        else
        {
	  B(offset) = b_ij;
          B.jcol(offset) = permute[jj];
          offset++;
        }
      }
    }
  }

  cerr << "Matrix B upset!" << endl;
//  B.print(cerr);

  // solve the system thrice
  PrBiCGStab solver;
  solver.setMaxIterations(ni);
  solver.setTolerance(tolerance_);
  solver.solve(B,uvec,b1);
  solver.solve(B,vvec,b2);
  solver.solve(B,wvec,b3);

  // Copy the results to the u and v arrays.
  uvw.resize(3*n);

  for(i=0; i<n; i++)
  {
    if(!isFixed(i,fixedPnts))
    {
      uvw[3*i+0] = uvec(permute[i]);
      uvw[3*i+1] = vvec(permute[i]);
      uvw[3*i+2] = wvec(permute[i]);
    }
    else
    {
      Vector3D p=g_->get3dNode(i);
      uvw[3*i+0] = p.x();
      uvw[3*i+1] = p.y();
      uvw[3*i+2] = p.z();
    }
  }

  cerr << "System solved!" << endl;

  // project parameter values to a unit sphere around the barycentre of the
  // fixed points

  Vector3D barycentre = 0.25 * (g_->get3dNode(fixedPnts[0]) +
				  g_->get3dNode(fixedPnts[1]) +
				  g_->get3dNode(fixedPnts[2]) +
				  g_->get3dNode(fixedPnts[3]) );

  for(i=0; i<n; i++)
  {
    Vector3D vec = Vector3D(uvw[3*i+0], uvw[3*i+1], uvw[3*i+2]) - barycentre;
    vec.normalize();
    uvw[3*i+0] = vec.x();
    uvw[3*i+1] = vec.y();
    uvw[3*i+2] = vec.z();
  }
  
                      
  return true;
}


//-----------------------------------------------------------------------------
bool PrParametrizeInt::new_parametrize3d(vector<int>& fixedPnts, 
					 vector<double>& uvw)
//-----------------------------------------------------------------------------
//   Parametrize the nodes of the 3D grapg g_ except those
//   with indices in "fixedPnts". The parameterization is done
//   in 3D and the result is returned as vector "uvw"
{
  // To find the internal (u,v,w) nodes, we solve a linear system
  // thrice: once for u, v and w. The matrix A is the same
  // in both cases. Only the right-hand side changes.

  // There is one equation for each interior node,
  // so we solve Au=b1, Av=b2 and Av=b3 where A is ni * ni and ni is the
  // number of interior nodes.

  int n = g_->getNumNodes();

  //
  // setup the "full" matrix A (i.e. without boundary vertices)
  //
  //  remark: this will only work for closed meshes!!!
  //          (because "makeWeights" only works for interior nodes)
  //

  // find number of non-zeros in the matrix A
  int i;
  int numNonZeros = 0;
  for(i=0; i<n; i++)
  {
    g_->getNeighbours(i, neighbours_);
    numNonZeros += ((int)neighbours_.size() + 1);
  }

  PrMatSparse A (n, n, numNonZeros);
  int offset = 0;

  for(i=0; i<n; i++)
  {
    A.irow(i) = offset;
    A(offset) = 1.0;
//A(offset) = 0.96;
    A.jcol(offset) = i;
    offset++;

    g_->getNeighbours(i,neighbours_);
    makeWeights(i); // Ignoring return value.
    
    for(size_t j=0; j<neighbours_.size(); j++)
    {
      A(offset) = -weights_[j];
      A.jcol(offset) = neighbours_[j];
      offset++;
    }
  }

  PrMatSparse B (0,0,0);
  A.matProd(A,B);
  B.matProd(B,A);

  // replace the lines of the matrix that correspond to fixed
  // nodes by a transposed unit vector
  int fix_size = (int)fixedPnts.size();
  for (i=0; i<fix_size; i++)
    for(int j=A.irow(fixedPnts[i]); j<A.irow(fixedPnts[i]+1); j++)
      if (A.jcol(j) == i)
	A(j) = 1.0;
      else
	A(j) = 0.0;

  cerr << "Matrix A upset!" << endl;
  cerr << A.irow(n) << " of " << n*n << " entries are non-zero (";
  cerr << A.irow(n)/n/n*100 << " %)" << endl;
     
  // Set start vector to the original xyz-coordinates.

  PrVec uvec(n);
  PrVec vvec(n);
  PrVec wvec(n);

  for(i=0; i<n; ++i)
  {
    Vector3D p=g_->get3dNode(i);
    uvec(i) = p.x();
    vvec(i) = p.y();
    wvec(i) = p.z();
  }

  // Set right hand side to zero except for the fixed points, since
  // they are set to the original xyz-coordinate

  PrVec b1(n, 0.0);
  PrVec b2(n, 0.0);
  PrVec b3(n, 0.0);
  
  for (i=0; i<fix_size; ++i)
  {
    b1(fixedPnts[i]) = uvec(fixedPnts[i]);
    b2(fixedPnts[i]) = vvec(fixedPnts[i]);
    b3(fixedPnts[i]) = wvec(fixedPnts[i]);
  }

  // solve the system thrice
  PrBiCGStab solver;
  solver.setMaxIterations(n);
  solver.setTolerance(tolerance_);

  solver.solve(A,uvec,b1);

  double cpu_time = solver.getCPUTime();
  int noIts = solver.getItCount();
  bool converged = solver.converged();
//   cout << "unknowns = " << n << "  cpu_time = " << cpu_time
//        << "  no_its = " << noIts << "  converged = " << converged << endl;

  solver.solve(A,vvec,b2);

  cpu_time = solver.getCPUTime();
  noIts = solver.getItCount();
  converged = solver.converged();
//   cout << "unknowns = " << n << "  cpu_time = " << cpu_time
//        << "  no_its = " << noIts << "  converged = " << converged << endl;

  solver.solve(A,wvec,b3);

  cpu_time = solver.getCPUTime();
  noIts = solver.getItCount();
  converged = solver.converged();
//   cout << "unknowns = " << n << "  cpu_time = " << cpu_time
//        << "  no_its = " << noIts << "  converged = " << converged << endl;


//   cerr << "System solved!" << endl;

  // Copy the results to the u and v arrays.
  uvw.resize(3*n);

  for(i=0; i<n; i++)
  {
    uvw[3*i+0] = uvec(i);
    uvw[3*i+1] = vvec(i);
    uvw[3*i+2] = wvec(i);
  }

  return true;
}



//-----------------------------------------------------------------------------
void PrParametrizeInt::findFixedPntsFromXYZ(vector<int>& fixedPnts)
//-----------------------------------------------------------------------------
//   This is a simple routine which finds the indices
//   fixedPnts[0..3] of the four vertices of the graph
//   whose (x,y,z) points are the furthest in the directions
//   (-1/-1/1), (1/1/1), (-1/1/-1), and (1/-1/-1) in that order.
//   The indices can be used as fixed points for 
//   parametrising in 3D.
{
  int k;
  fixedPnts.resize(4);

  Vector3D dir[4];
  dir[0] = Vector3D(-1,-1,1);
  dir[1] = Vector3D(1,1,1);
  dir[2] = Vector3D(-1,1,-1);
  dir[3] = Vector3D(1,-1,-1);

  Vector3D p0 = g_->get3dNode(0);
  double dp[4];
  for(k=0; k<4; k++) {
    fixedPnts[k] = 0;
    dp[k] = p0 * dir[k];
  }

  for(int i=1; i< g_->getNumNodes(); i++)
  {
    Vector3D p = g_->get3dNode(i);
    for(k=0; k<4; k++)
      {
	double newdp = p * dir[k];
	if (newdp > dp[k]) {
	  fixedPnts[k] = i; 
	  dp[k] = newdp;
	}
      }
  }

  cout << "the following 4 points have been fixed: " << endl;
  for(k=0; k<4; k++) {
    Vector3D p = g_->get3dNode(fixedPnts[k]);
    cout << fixedPnts[k] << " : " << p.x() << "," << p.y() << "," << p.z() << endl;
  }
}


//-----------------------------------------------------------------------------
void PrParametrizeInt::computeWeights()
//-----------------------------------------------------------------------------
//   computes and stores the weights for each node
{
  cout << "Computing all weights...";

  for(int i=0; i<g_->getNumNodes(); i++) {
    g_->getNeighbours(i,neighbours_);
    allNeighbours_.push_back(neighbours_);
    makeWeights(i); // Ignoring return value.
    allWeights_.push_back(weights_);
  }

  cout << "done" << endl;
}


//-----------------------------------------------------------------------------
void PrParametrizeInt::smooth(int nmb, vector<int>& fixedPnts)
//-----------------------------------------------------------------------------
//   performs "nmb" Gauss-Seidel smoothing steps on the sphere
{
  for(int j=0; j<nmb; j++) {
    for(int i=0; i<g_->getNumNodes(); i++) {
      if (!isFixed(i, fixedPnts)) {
	Vector3D node_i (0.0, 0.0, 0.0);
	
	// average node i and scale it back to the unit sphere
	for(size_t k=0; k<allWeights_[i].size(); k++) {
	  node_i += allWeights_[i][k] * g_->get3dNode(allNeighbours_[i][k]);
	}
	node_i.normalize();
	
	g_->set3dNode(i, node_i);
      }
    }
  }
}

