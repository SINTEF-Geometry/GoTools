#ifndef _FTTANGPRIORITY_H
#define _FTTANGPRIORITY_H

namespace Go
{

  /// Whether the group of objects (surfaces) are a master, a
  /// slave or not specified. Related to tangent plane continuity
  /// between groups of surfaces
  enum ftTangPriority
  {
    ftNoType = 0,
    ftMaster,
    ftSlave
  };

} // namespace Go

#endif //  _FTTANGPRIORITY_H
