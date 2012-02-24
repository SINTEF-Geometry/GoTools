#ifndef _FACECONNECTIVITY_H_
#define _FACECONNECTIVITY_H_

#include <vector>

using std::vector;
using std::pair;

namespace Go
{

/// A structure storing the connectivity information between
/// two adjacent faces
template <class edgeType>
struct FaceConnectivity
{
public:
  /// Edge of first face along the common boundary
  edgeType *e1_;
  /// Edge of second face along the common boundary
  edgeType *e2_;
  /// The status is:
  /// 0 : edges join smoothly. \f$G^1\f$.
  /// 1 : edges join, but normals are slightly discontinous. A kink.
  /// 2 : edges join, but the normals are discontinous. \f$G^0\f$.
  /// 3 : edges almost join. A gap.
  /// 4 : edges are totally discontinous.
  /// The minimal tpTopologicalInfo has a one-element status vector
  /// and a two-element parameters vector
  vector<int> status_;
  /// Parameter intervals limiting the areas of the found state of continuity 
  vector< pair<double, double> > parameters_;

  /// Constructor
  FaceConnectivity(edgeType* e1, edgeType *e2)
  {
    e1_ = e1;
    e2_ = e2;
  }

  /// Reset edge info
  void setEdges(edgeType* e1, edgeType *e2)
  {
    e1_ = e1;
    e2_ = e2;
  }

  /// The highest continuity between the two faces
  int BestStatus() const
  {
    int s = 4;
    for (size_t i = 0; i < status_.size(); ++i)
      if (status_[i] < s)
	s = status_[i];
    return s;
  }
  /// The lowest continuity between the two faces
  int WorstStatus() const
  {
    int s = 0;
    for (size_t i = 0; i < status_.size(); ++i)
      if (status_[i] > s) s = status_[i];
    return s;
  }
};

} // namespace Go
#endif  //_FACECONNECTIVITY_H_
