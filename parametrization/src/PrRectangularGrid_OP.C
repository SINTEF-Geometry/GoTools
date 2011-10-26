/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/


/********************************************************************
 FILENAME    : PrRectangularGrid_OP.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Jan. 97
 DESCRIPTION : Implementation of methods in the class PrRectangularGrid_OP.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrRectangularGrid_OP.h"
#include "GoTools/parametrization/PrParamUtil.h"

// PRIVATE MEMBER FUNCTIONS

// PUBLIC MEMBER FUNCTIONS

//-----------------------------------------------------------------------------
PrRectangularGrid_OP::PrRectangularGrid_OP(int numcols, int numrows,
					   const double* grid)
//-----------------------------------------------------------------------------
//   Construct the PrRectangularGrid_OP from an array of nodes.
    : ncols_(numcols), nrows_(numrows)
{
  int tot = numcols * numrows;
  grid_.reserve(tot);
  uvgrid_.reserve(tot);
  for(int j=0; j< tot; j++)
  {
    grid_.push_back(Vector3D(grid[3*j], grid[3*j+1], grid[3*j+2]));
    uvgrid_.push_back(Vector2D(0.0,0.0));
  }
}

//-----------------------------------------------------------------------------
PrRectangularGrid_OP::PrRectangularGrid_OP(int numcols, int numrows,
					   const double* grid,
					   const double* uvgrid)
//-----------------------------------------------------------------------------
//   Construct the PrRectangularGrid_OP from an array of nodes.
    : ncols_(numcols), nrows_(numrows)
{
  int tot = numcols * numrows;
  grid_.reserve(tot);
  uvgrid_.reserve(tot);
  for(int j=0; j< tot; j++)
  {
    grid_.push_back(Vector3D(grid[3*j], grid[3*j+1], grid[3*j+2]));
    uvgrid_.push_back(Vector2D(uvgrid[2*j], uvgrid[2*j+1]));
  }
}

//-----------------------------------------------------------------------------
PrRectangularGrid_OP::~PrRectangularGrid_OP()
//-----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
PrRectangularGrid_OP::getNeighbours(int k, vector<int>& neighbours) const
//-----------------------------------------------------------------------------
//   Return the neighbours of the k-th node in:
//   1. any anticlockwise order if the k-th node is an interior node
//   2. the unique anticlockwise order if the k-th node is a boundary node.
{
  neighbours.clear();
  neighbours.reserve(4);
  int i,j;
  graphToGrid(k,i,j);

  if(i == 0 && j == 0)
  {
      neighbours.push_back(gridToGraph(i+1,j));
      neighbours.push_back(gridToGraph(i,j+1));
  }
  else if(i == ncols_-1 && j == 0)
  {
      neighbours.push_back(gridToGraph(i,j+1));
      neighbours.push_back(gridToGraph(i-1,j));
  }
  else if(i == ncols_-1 && j == nrows_-1)
  {
      neighbours.push_back(gridToGraph(i-1,j));
      neighbours.push_back(gridToGraph(i,j-1));
  }
  else if(i == 0 && j == nrows_-1)
  {
      neighbours.push_back(gridToGraph(i,j-1));
      neighbours.push_back(gridToGraph(i+1,j));
  }

  else if(j == 0)
  {
      neighbours.push_back(gridToGraph(i+1,j));
      neighbours.push_back(gridToGraph(i,j+1));
      neighbours.push_back(gridToGraph(i-1,j));
  }
  else if(i == ncols_-1)
  {
      neighbours.push_back(gridToGraph(i,j+1));
      neighbours.push_back(gridToGraph(i-1,j));
      neighbours.push_back(gridToGraph(i,j-1));
  }
  else if(j == nrows_-1)
  {
      neighbours.push_back(gridToGraph(i-1,j));
      neighbours.push_back(gridToGraph(i,j-1));
      neighbours.push_back(gridToGraph(i+1,j));
  }
  else if(i == 0)
  {
      neighbours.push_back(gridToGraph(i,j-1));
      neighbours.push_back(gridToGraph(i+1,j));
      neighbours.push_back(gridToGraph(i,j+1));
  }
  else // interior node
  {
    neighbours.push_back(gridToGraph(i+1,j));
    neighbours.push_back(gridToGraph(i,j+1));
    neighbours.push_back(gridToGraph(i-1,j));
    neighbours.push_back(gridToGraph(i,j-1));
  }
  return;
}

//----------------------------------------------------------------------------
bool PrRectangularGrid_OP::isBoundary(int k) const
//-----------------------------------------------------------------------------
//   Given a node and its leading edge,
//   return true if the node is a boundary node and false otherwise.
{
  int i,j;
  graphToGrid(k,i,j);
  if(i == 0 || j == 0 || i == ncols_-1 || j == nrows_-1) return true;
  return false;
}

//-----------------------------------------------------------------------------
void PrRectangularGrid_OP::setDim(int numcols, int numrows)
//-----------------------------------------------------------------------------
// Set the dimensions of the grid. 
{
  int tot = numcols * numrows;
  grid_.resize(tot);
  uvgrid_.resize(tot);
  ncols_ = numcols;
  nrows_ = numrows;
}

//----------------------------------------------------------------------------
void PrRectangularGrid_OP::setXYZVertices(const double *grid)
//-----------------------------------------------------------------------------
// Reset all the xyz points from an array
{
  int tot = ncols_ * nrows_;
  for(int j=0; j< tot; j++)
  {
    grid_[j] = Vector3D(grid[3*j], grid[3*j+1], grid[3*j+2]);
  }
}

//----------------------------------------------------------------------------
void PrRectangularGrid_OP::setUVVertices(const double *grid)
//-----------------------------------------------------------------------------
// Reset all the uv points from an array
{
  int tot = ncols_ * nrows_;
  for(int j=0; j< tot; j++)
  {
    uvgrid_[j] = Vector2D(grid[2*j], grid[2*j+1]);
  }
}

//----------------------------------------------------------------------------
void PrRectangularGrid_OP::getCorners(int& c1, int& c2, int& c3, int& c4) const
//-----------------------------------------------------------------------------
//   Get the indices of the four corner nodes in an anticlockwise
//   direction, starting with the bottom left hand corner.
{
  c1 = gridToGraph(0,0);
  c2 = gridToGraph(ncols_-1,0);
  c3 = gridToGraph(ncols_-1,nrows_-1);
  c4 = gridToGraph(0,nrows_-1);
}

//-----------------------------------------------------------------------------
bool PrRectangularGrid_OP::getUVTriangle(double u, 
					 double v,
					 int& i, 
					 int&j, 
					 bool& right,
					 double& tau0, 
					 double& tau1, 
					 double& tau2) const
//-----------------------------------------------------------------------------
{
  int ii,jj;

  for(jj=0; jj<(nrows_-1); jj++)
  {
    for(ii=0; ii<(ncols_-1); ii++)
    {
      if(getTriangle(u, v, ii, jj, right, tau0, tau1, tau2))
      {
        i = ii;
        j = jj;
        return true;
      }
    }
  }

  i = -1;
  j = -1;
  return false;
}

//-----------------------------------------------------------------------------
bool PrRectangularGrid_OP::getUVTriangle(double& u, double& v,
                             int ilast, int jlast,
                             int& i, int&j, bool& right,
                             double& tau0, double& tau1, double& tau2) const
//-----------------------------------------------------------------------------
{
  if(getTriangle(u,v,ilast,jlast,right,tau0,tau1,tau2)) {
      i = ilast;
      j = jlast;
      return true;
  }
  if (ilast < ncols_-2) {
      if(getTriangle(u,v,ilast+1,jlast,right,tau0,tau1,tau2)) {
	  i = ilast+1;
	  j = jlast;
	  return true;
      }
  }
  if (jlast < nrows_-2) {
      if(getTriangle(u,v,ilast,jlast+1,right,tau0,tau1,tau2)) {
	  i = ilast;
	  j = jlast+1;
	  return true;
      }
  }
  if (ilast > 0) {
      if(getTriangle(u,v,ilast-1,jlast,right,tau0,tau1,tau2)) {
	  i = ilast-1;
	  j = jlast;
	  return true;
      }
  }
  if (jlast > 0) {
      if(getTriangle(u,v,ilast,jlast-1,right,tau0,tau1,tau2)) {
	  i = ilast;
	  j = jlast-1;
	  return true;
      }
  }
  if (ilast == ncols_-2) {
      if(getTriangle(u,v,0,jlast,right,tau0,tau1,tau2)) {
	  i = 0;
	  j = jlast;
	  return true;
      }
  }

  return getUVTriangle(u,v,i,j,right,tau0,tau1,tau2);
}

//-----------------------------------------------------------------------------
bool PrRectangularGrid_OP::getTriangle(double u, double v,
				       int ii, int jj, bool& right,
				       double& tau0, double& tau1, double& tau2) const
//-----------------------------------------------------------------------------
{
  // @afr: This seems to assume that (ii, jj) is the lower left corner
  // of the gridcell from which we want a triangle. That is, ii must
  // be less than n_-1 and jj must be less than m_-1
  int k1,k2,k3,k4;
  k1 = gridToGraph(ii,jj);
  k2 = gridToGraph(ii+1,jj);
  k3 = gridToGraph(ii,jj+1);
  k4 = gridToGraph(ii+1,jj+1);
  baryCoords(u,v,getU(k1),getV(k1),getU(k2),getV(k2),getU(k4),getV(k4),
             tau0, tau1, tau2);
  if(tau0 >= -1.0e-10 && tau1 >= -1.0e-10 && tau2 >= -1.0e-10)
  {
    right = true;
    return true;
  }

  baryCoords(u,v,getU(k1),getV(k1),getU(k4),getV(k4),getU(k3),getV(k3),
             tau0, tau1, tau2);
  if(tau0 >= -1.0e-10 && tau1 >= -1.0e-10 && tau2 >= -1.0e-10)
  {
    right = false;
    return true;
  }

  return false;
}

//-----------------------------------------------------------------------------
void PrRectangularGrid_OP::print(std::ostream& os)
//-----------------------------------------------------------------------------
{
  os << ncols_ << ' ' << nrows_ << '\n';
  int i,j,k=0;
  for(i=0; i<nrows_; i++)
  {
    for(j=0; j<ncols_; j++)
    {
//	os << uvgrid_[k] << grid_[k];
      uvgrid_[k].write(os);
      grid_[k].write(os);
	k++;
    }
    os << std::endl;
  }
}

//-----------------------------------------------------------------------------
void PrRectangularGrid_OP::scan(std::istream& is)
//-----------------------------------------------------------------------------
{
  is >> ncols_ >> nrows_;
  int tot = ncols_*nrows_;
  grid_.resize(tot);
  uvgrid_.resize(tot);
  int i,j,k=0;
  for(i=0; i<nrows_; i++)
  {
    for(j=0; j<ncols_; j++)
    {
      uvgrid_[k].read(is);
      grid_[k].read(is);
//	is >> uvgrid_[k] >> grid_[k];
	k++;
    }
  }
}

//-----------------------------------------------------------------------------
void PrRectangularGrid_OP::printRawData(std::ostream& os)
//-----------------------------------------------------------------------------
{
  os << ncols_ << ' ' << nrows_ << '\n';
  int i,j,k=0;
  for(i=0; i<nrows_; i++)
  {
    for(j=0; j<ncols_; j++)
    {
//	os << grid_[k];
      grid_[k].write(os);
	k++;
    }
    os << std::endl;
  }
}

//-----------------------------------------------------------------------------
void PrRectangularGrid_OP::scanRawData(std::istream& is)
//-----------------------------------------------------------------------------
{
  is >> ncols_ >> nrows_;
  int tot = ncols_*nrows_;
  grid_.resize(tot);
  uvgrid_.resize(tot);
  int i,j,k=0;
  for(i=0; i<nrows_; i++)
  {
    for(j=0; j<ncols_; j++)
    {
      grid_[k].read(is);
//	is >> grid_[k];
	uvgrid_[k] = Vector2D(0.0,0.0);
	k++;
    }
  }
}

