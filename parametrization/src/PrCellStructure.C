/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1999 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/


/********************************************************************
 FILENAME    : PrCellStructure.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Sep. 99
 DESCRIPTION : Implementation of methods in the class PrCellStructure.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrCellStructure.h"
#include <queue>
using std::cout;
using std::endl;


class HeapNode {

public:
  int     idx_;
  double  key_;

  HeapNode( int idx, double key) { idx_ = idx; key_ = key;}
  ~HeapNode(){}

  bool operator >   ( const HeapNode& x) const { return (this->key_
>   x.key_);}
  bool operator >=  ( const HeapNode& x) const { return (this->key_
>=  x.key_);}
};
typedef std::priority_queue< HeapNode, vector<HeapNode>, std::greater<HeapNode> > HeapType;

// PRIVATE MEMBER FUNCTIONS

//-----------------------------------------------------------------------------
void PrCellStructure::makeCellStructure()
//-----------------------------------------------------------------------------
{
  // Find bounding box
  double maxarr[3];
  min_[0] = xyz_[0].x();
  maxarr[0] = xyz_[0].x();
  min_[1] = xyz_[0].y();
  maxarr[1] = xyz_[0].y();
  min_[2] = xyz_[0].z();
  maxarr[2] = xyz_[0].z();
  int i;
  int xyz_size = (int)xyz_.size();
  for(i=1; i< xyz_size; i++)
  {
    if(xyz_[i].x() < min_[0]) min_[0] = xyz_[i].x();
    else if(xyz_[i].x() > maxarr[0]) maxarr[0] = xyz_[i].x();

    if(xyz_[i].y() < min_[1]) min_[1] = xyz_[i].y();
    else if(xyz_[i].y() > maxarr[1]) maxarr[1] = xyz_[i].y();

    if(xyz_[i].z() < min_[2]) min_[2] = xyz_[i].z();
    else if(xyz_[i].z() > maxarr[2]) maxarr[2] = xyz_[i].z();
  }


  // Choose a good cell size;
  cell_size_ = std::max( (maxarr[0] - min_[0]) / (double)max_no_cells_,
                   (maxarr[1] - min_[1]) / (double)max_no_cells_ );
  cell_size_ = std::max( cell_size_,
                   (maxarr[2] - min_[2]) / (double)max_no_cells_ );
  

  // Decide how many cells are needed in each direction
  ncells_[0] = (int)((maxarr[0] - min_[0]) / cell_size_) + 1;
  ncells_[1] = (int)((maxarr[1] - min_[1]) / cell_size_) + 1;
  ncells_[2] = (int)((maxarr[2] - min_[2]) / cell_size_) + 1;


  //allocate space for cells
  ind_.resize(ncells_[0]*ncells_[1]*ncells_[2]);

  // Now start filling in the cell information
  int ii,i1,j1,k1;
  for(i=0; i< xyz_size; i++)
  {
    // Find which cell the point lies in
    whichCell(xyz_[i],i1,j1,k1);
    //cout << oform("(i1,j1,k1) = (%d,%d,%d)\n",i1,j1,k1);
    ii = getI(i1,j1,k1);
    //cout << oform("ii = %d\n",ii);
    ind_[ii].push_back(i);
  }

  size_t min_num = xyz_.size();
  int ind_size = (int)ind_.size();
  size_t max_num = 0;
  for (i=0; i<ind_size; i++ )
  {
    if(ind_[i].size() > max_num) max_num = ind_[i].size();
    else if(ind_[i].size() < min_num) min_num = ind_[i].size();
  }


}

// PUBLIC MEMBER FUNCTIONS

//-----------------------------------------------------------------------------
PrCellStructure::PrCellStructure(int n, double* xyz_points, int max_no_cells)
//-----------------------------------------------------------------------------
{
  xyz_.resize(n);
  int j;
  for(j=0; j< n; j++)
  {
    xyz_[j].x() = xyz_points[3*j];
    xyz_[j].y() = xyz_points[3*j+1];
    xyz_[j].z() = xyz_points[3*j+2];
  }
  max_no_cells_ = max_no_cells;
  makeCellStructure();
}

//-----------------------------------------------------------------------------
void PrCellStructure::attach(int n, const double* xyz_points)
//-----------------------------------------------------------------------------
{
  xyz_.resize(n);
  int j;
  for(j=0; j< n; j++)
  {
    xyz_[j].x() = xyz_points[3*j];
    xyz_[j].y() = xyz_points[3*j+1];
    xyz_[j].z() = xyz_points[3*j+2];
  }
  makeCellStructure();
}

//----------------------------------------------------------------------------
void
PrCellStructure::getBall(const Vector3D& p, double radius,
                               vector<int>& neighbours, int notP) const
//-----------------------------------------------------------------------------
//   Return all points within the ball of radius radius around
//   the point p. Don't include p itself if notP = 1.
{ 
  neighbours.clear();
  double r2 = radius * radius, dist2;
  // find which cell p lies in:
  int i1,j1,k1;
  whichCell(p,i1,j1,k1);
  int imin[3],imax[3];
  int diff = (int)(radius/cell_size_) + 1;
  imin[0] = std::max(i1-diff,0);
  imax[0] = std::min(i1+diff,ncells_[0]-1);
  imin[1] = std::max(j1-diff,0);
  imax[1] = std::min(j1+diff,ncells_[1]-1);
  imin[2] = std::max(k1-diff,0);
  imax[2] = std::min(k1+diff,ncells_[2]-1);

  int iq,ii;
  for(int k=imin[2]; k<=imax[2]; k++)
    for(int j=imin[1]; j<=imax[1]; j++)
      for(int i=imin[0]; i<=imax[0]; i++)
      {
        ii = getI(i,j,k);
        for(size_t ell=0; ell<ind_[ii].size(); ell++)
        {
          iq = ind_[ii][ell];
          dist2 = xyz_[iq].dist2(p);
          if(dist2 <= r2)
          {
            if(notP == 0 || dist2 > 0) neighbours.push_back(iq);
          }
        }
      }
}


//----------------------------------------------------------------------------
void
PrCellStructure::getKNearest(const Vector3D& p, int k,
                               vector<int>& neighbours, int notP) const
//-----------------------------------------------------------------------------
//   Return k nearest points to the point p.
//   Don't include p itself if notP = 1.
{
  neighbours.clear();
  if(k > int(xyz_.size()) - 1) return; // max k is xyz_.size() - 1

  // The following method is not optimal.
  // One should grow out in layers of cubes rather than in
  // larger and larger spheres, but this'll do as a first try.

  double r = cell_size_;

  vector<int> current_nghrs;
  getBall(p,r,current_nghrs,notP);
 
  while(int(current_nghrs.size()) < k)
  {
    r *= 2.0;
    getBall(p,r,current_nghrs,notP);
  }

  // Now we've got too many points in current_nghrs.
  // We now sort the points in distance from p and keep the
  // nearest k.


  int n = (int)current_nghrs.size();
  double dist2;
//  PrHeap heap(n);
//  std::priority_queue< HeapNode , std::vector<HeapNode>, std::greater<HeapNode> > heap(n);
  HeapType heap;
  int i;
  for(i = 0; i < n; i++)
  {
    dist2 = xyz_[current_nghrs[i]].dist2(p);
    heap.push(HeapNode(i,dist2));
  }

  int index;
  for(i = 0; i < k; i++)
  {
    index = heap.top().idx_;
    heap.pop();
//    heap.pop(key,index);
        // numbers are popped in increasing order
    //s_o << "  " << key << "  " << index << endl;
    neighbours.push_back(current_nghrs[index]);
  }
}


//----------------------------------------------------------------------------
void
PrCellStructure::getIJK(int ii, int& i, int& j, int& k) const
//-----------------------------------------------------------------------------
{
  k = ii / (ncells_[0] * ncells_[1]);
  int jj = ii - k * (ncells_[0] * ncells_[1]);
  j = jj / ncells_[0];
  i = jj - j * ncells_[0];
}

//----------------------------------------------------------------------------
void
PrCellStructure::whichCell(const Vector3D& xyz, int& i, int& j, int& k) const
//-----------------------------------------------------------------------------
{
  i = (int)((xyz.x() - min_[0]) / cell_size_);
  j = (int)((xyz.y() - min_[1]) / cell_size_);
  k = (int)((xyz.z() - min_[2]) / cell_size_);
}

//-----------------------------------------------------------------------------
void PrCellStructure::print(ostream& os)
//-----------------------------------------------------------------------------
{
  os << xyz_.size() << endl;
  size_t i;
  for(i=0; i<xyz_.size(); i++)
  {
    os << xyz_[i].x() << " " << xyz_[i].y() << " " << xyz_[i].z() << "\n";
  }
}

//-----------------------------------------------------------------------------
void PrCellStructure::scan(istream& is)
//-----------------------------------------------------------------------------
{
  int numpnts;
  is >> numpnts;
  xyz_.resize(numpnts);
  int i;
  for(i=0; i<numpnts; i++)
  {
    is >> xyz_[i].x() >> xyz_[i].y() >> xyz_[i].z();
  }
  makeCellStructure();
}

