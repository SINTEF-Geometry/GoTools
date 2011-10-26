#ifndef _FACECONNECTIVITY_H_
#define _FACECONNECTIVITY_H_

#include <vector>

using std::vector;
using std::pair;

template <class edgeType>
struct FaceConnectivity
{
public:

  // The status is:
  // 0 : edges join smoothly. G1.
  // 1 : edges join, but normals are slightly discontinous. A kink.
  // 2 : edges join, but the normals are discontinous. G0.
  // 3 : edges almost join. A gap.
  // 4 : edges are totally discontinous.
  // The minimal tpTopologicalInfo has a one-element status vector
  // and a two-element parameters vector
  edgeType *e1_;
  edgeType *e2_;
  vector<int> status_;
  vector< pair<double, double> > parameters_;

  // Constructor
  FaceConnectivity(edgeType* e1, edgeType *e2)
  {
    e1_ = e1;
    e2_ = e2;
  }

  // Reset edge info
  void setEdges(edgeType* e1, edgeType *e2)
  {
    e1_ = e1;
    e2_ = e2;
  }

  int BestStatus() const
  {
    int s = 4;
    for (size_t i = 0; i < status_.size(); ++i)
      if (status_[i] < s)
	s = status_[i];
    return s;
  }
  int WorstStatus() const
  {
    int s = 0;
    for (size_t i = 0; i < status_.size(); ++i)
      if (status_[i] > s) s = status_[i];
    return s;
  }
};

#endif  //_FACECONNECTIVITY_H_
