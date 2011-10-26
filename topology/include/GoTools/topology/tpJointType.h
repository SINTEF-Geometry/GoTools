#ifndef _TPJOINTTYPE_H
#define _TPJOINTTYPE_H

enum tpJointType
{
  JOINT_G1   = 0,  ///< G1 continous (angle diff. less than kink tolerance)
  JOINT_KINK = 1,  ///< Continous (angle diff. between kink and bend tolerance)
  JOINT_G0   = 2,  ///< Continous (angle diff. greater than bend tolerance)
  JOINT_GAP  = 3,  ///< Discontinous (distance between gap and neighbour tol.)
  JOINT_DISC = 4,  ///< Discontinous (distance greater than neighbour tolerance)
  JOINT_NONE = 5   ///< No joint, this was the last segment
};

#endif //  _TPJOINTTYPE_H
