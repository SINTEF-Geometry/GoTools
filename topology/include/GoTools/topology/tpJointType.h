#ifndef _TPJOINTTYPE_H
#define _TPJOINTTYPE_H

namespace Go
{

enum tpJointType
{
  ///< \f$G^1\f$ continous (angle diff. less than kink tolerance)
  JOINT_G1   = 0, 
  ///< Continous (angle diff. between kink and bend tolerance)
  JOINT_KINK = 1,  
  ///< Continous (angle diff. greater than bend tolerance)
  JOINT_G0   = 2,  
  ///< Discontinous (distance between gap and neighbour tol.)
  JOINT_GAP  = 3, 
  ///< Discontinous (distance greater than neighbour tolerance) 
  JOINT_DISC = 4, 
  ///< No joint, this was the last segment
  JOINT_NONE = 5   
};

} // namespace GO

#endif //  _TPJOINTTYPE_H
