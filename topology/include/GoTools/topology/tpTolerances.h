#ifndef _TPTOLERANCES_H_
#define  _TPTOLERANCES_H_
/** tpTolerances -  Short description.
 * Detailed description.
 */
struct tpTolerances
{
public:
    double gap;
    double neighbour;
    double kink;
    double bend;
    tpTolerances(double g, double n, double k, double b)
	: gap(g), neighbour(n), kink(k), bend(b)
    {}
    tpTolerances(const tpTolerances& tol)
	: gap(tol.gap), neighbour(tol.neighbour), kink(tol.kink),
          bend(tol.bend)
    {}
};

#endif //  _TPTOLERANCES_H_
