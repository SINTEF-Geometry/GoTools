#ifndef _TPTOLERANCES_H_
#define  _TPTOLERANCES_H_

namespace Go
{

/** tpTolerances -  Tolerances used in adjacency analysis for faces
 */
struct tpTolerances
{
public:
  /// Tolerance for when two faces are assumed to be C0 continous
    double gap;
  /// Tolerance for when two faces are assumed to be neighbours
    double neighbour;
  /// Tolerance for when two adjacent faces are assumed to be C1 continous
    double kink;
  /// Tolerance for when two adjacent faces are assumed to have an 
  /// intentially smooth connection
    double bend;
    tpTolerances(double g, double n, double k, double b)
	: gap(g), neighbour(n), kink(k), bend(b)
    {}
    tpTolerances(const tpTolerances& tol)
	: gap(tol.gap), neighbour(tol.neighbour), kink(tol.kink),
          bend(tol.bend)
    {}
};

} // namespace Go

#endif //  _TPTOLERANCES_H_
