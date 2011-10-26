#ifndef _SMOOTHSURFSET_H_
#define _SMOOTHSURFSET_H_


//   -----------------------------------------------------------------------
//      Interface file for class SmoothSurfSET
//   -----------------------------------------------------------------------
//
//       Used to modify a tensor product B-spline surface with respect
//       to conditions on smoothness, editing constraints and boundary
//       conditions.
//
//       Implementation of the member functions are given in the
//       following files:
//
//          1. SmoothSurfSet.C
//
//   -----------------------------------------------------------------------
//    Written by: Vibeke Skytt                            04.2002.
//   -----------------------------------------------------------------------

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/creators/ConstraintDefinitions.h"
#include <vector>

namespace Go
{

    /// This class modifies a set of tensor product B-spline surfaces with respect
    /// to conditions on smoothness, editing constraints and boundary conditions.
class SmoothSurfSet
{


public:
    /// Default constructor. Initializes class variable
    SmoothSurfSet();

    /// Constructor. Initializes class variable.
    /// \param copy_coefs whether to copy coefs in input surface or work directly on them.
    SmoothSurfSet(bool copy_coefs);

    /// Destructor.
    virtual
    ~SmoothSurfSet();

    /// Initializes data given by an intermediate surface.
    /// \param insf input set of surfaces to be smoothed.
    /// \param coef_known whether the coefs are free to be altered.
    ///                   0: not known, 1: known
    ///                   >= kpointer_ = 3: the coefficients indexed by ki &
    ///                                     coef_known[ki] - kpointer_ should be equal.
    ///                   Size of coef_known equal to that of insf, size of coef_known[ki]
    ///                   equals the number of control points in insf[ki].
    /// \param num_side_constraints number of side constraints to the minimization problem
    /// \param has_normal_cond indicates if normal conditions will be given
    void attach(std::vector<std::shared_ptr<SplineSurface> >& insf, 
		std::vector<std::vector<int> >& coef_known,
		int num_side_constraints = 0, 
		int has_normal_cond = 0);

    // @@@ VSK. Can it be relevant to use different weights for 
    // different surfaces, and how to define the weights in that case?
    /// Compute the smoothing part of the equation system.
    /// All weights should be positive, and their sum not exceed 1.0.
    /// \param weight1 smoothing weight w.r.t. the 1st derivative.
    /// \param weight2 smoothing weight w.r.t. the 2nd derivative.
    /// \param weight3 smoothing weight w.r.t. the 3rd derivative.
    virtual
    void setOptimize(double weight1, double weight2, double weight3);

    /// Compute matrices for least squares approximation.
    /// First index of all vectors corresponds to the indexing of the attached surfaces.
    /// \param pnts points on surfaces to be approximated.
    ///             Stored as (for 3D): (x0, y0, z0, x1, y1, z1, ...)
    /// \param param_pnts the corresponding 2-dimensional parametric points.
    ///                   Stored as: (u0, v0, u1, v1, ...)
    /// \param pnt_weights each of the pnts is assigned a weight lying in
    ///                    the unit interval. 1.0 if all pnts are equally important.
    /// \param weight the contribution of the approximation of the pnts in the system.
    ///               weight should lie in the unit interval.
    void setLeastSquares(std::vector<std::vector<double> >& pnts,
			 std::vector<std::vector<double> >& param_pnts,
			 std::vector<std::vector<double> >& pnt_weights,
			 double weight);

    /// Compute matrices for approximation of normal directions.			 
    /// \param pnts normals in sf to be approximated.
    ///             Stored as (for 3D): (x0, y0, z0, x1, y1, z1, ...)
    /// \param param_pnts the corresponding 2-dimensional parametric points.
    ///                   Stored as: (u0, v0, u1, v1, ...)
    /// \param pnt_weights each of the pnts is assigned a weight lying in
    ///                    the unit interval. 1.0 if all pnts are equally important.
    /// \param weight the contribution of the approximation of the normals in the system.
    ///               weight should lie in the unit interval.
    /// \return status value: 0 = OK, 1 = warning (system not prepared for normal conditions).
    int setNormalCond(std::vector<std::vector<double> >& pnts,
		      std::vector<std::vector<double> >& param_pnts,
		      std::vector<std::vector<double> >& pnt_weights,
		      double weight);

    /// Compute the contribution to the equation system from the approximation
    /// of an original surface, i.e. the contribution of the original coefficients.
    /// \param weight the relative contribution of the original coefs.
    ///               Should lie in the unit interval.
    void approxOrig(double weight);

    /// Add linear side constraints to the system equation.
    /// \param constraints the linear side constraints between surface coefficients.
    void setSideConstraints(std::vector<sideConstraintSet>& constraints);

    // We may have side constraints which are not suitable for exact equality as
    // spline solution space may not be large enough. We therefore allow using least
    // squares to minimize the error.
    // This applies in particular to constraint involving higher order derivatives.
    // @@ Currently function uses newmat. As matrix is sparce we may do it another way.
    ///Use least squares to minimize the error given by the constraints.
    void setApproxSideConstraints(std::vector<sideConstraintSet>& constraints,
				  double weight);

    /// Solve equation system, and produce output surfaces.
    /// If failing to solve the routine may throw an exception.
    /// \param surfaces the output surface.
    /// \return 0 = OK, negative = failed solving system.
    int equationSolve(std::vector<std::shared_ptr<SplineSurface> >& surfaces);


protected:
    int kdim_;             // 3 if normal conditions are included, otherwise 1.
    int idim_;             // Dimension of geometry space. 
    int idim1_;            // Dimension of projective space.
    int ider_;             // Maximum derivative involved in the computations. 

    int kncond_;           // Number of unknown coefficents + number of constraints.
    int knconstraint_;     // Number of side constraints.

    std::vector<std::vector<int>::iterator> coef_known_;
    // Array indicating the status of coefficients, i.e.
    // free, fixed, not involved, equal to a given
    // coefficients.
    std::vector<std::vector<int> > pivot_; // Array giving the position of the
    // free coefficients in the equation system.

    // Parameters defining the spline space.
    std::vector<std::shared_ptr<SplineSurface> > srfs_;   // Pointer to input surfaces.

    const int copy_coefs_; // Whether we are to copy the coefs from input surfaces
                           // (or work directly on coefficients in input surfaces).

    // Parameters used to define the specific input spline surfaces.
    std::vector<std::vector<double> > coef_array_;

    //    std::vector<double>::iterator  scoef;   // Pointer to surface coefficients.      

    // Storage of the equation system.
    std::vector<double> gmat_;       // Matrix at left side of equation system.
    std::vector<double> gright_;     // Right side of equation system.      

    /// Given the value of non-zero B-spline functions, compute the value
    /// of the corresponding surface basis function (i.e. the products of
    /// the u and v basis functions).
    /// \param sb1 the basis values in the first parameter direction (u).
    /// \param sb2 the basis values in the second parameter direction (v).
    /// \param kk1 order in the first (u) direction.
    /// \param kk2 order in the second (v) direction.
    /// \param kleft1 index of the first knot interval in u-dir.
    /// \param kleft2 index of the first knot interval in v-dir.
    /// \param ider the number of derivatives to compute.
    /// \param sbasis the computed basis values in sf.
    ///               size = order_u()*order_v()*(ider + 1)*(ider + 1).
    ///               The space must be allocated on the outside.
    virtual
    void getBasis(double *sb1, double *sb2, int kk1, int kk2, 
		  int kleft1, int kleft2, int ider, double *sbasis);


private:

    /// Struct for storing integral information of the surface.
    typedef struct integralInfo
    {
	// Parameters used in integration
	std::vector<double> vec1, vec2;
	double ***integral1;  // Array used to store integrals of inner product
	// of derivatives of B-splines in 1. par. dir.
	// The 1st index runs over the derivatives, the 2nd & 3rd in u- and v-dir.
	double ***integral2;  // Array used to store integrals of inner product
	// of derivatives of B-splines in 2. par. dir.
	// The 1st index runs over the derivatives, the 2nd & 3rd in u- and v-dir.
	bool integralset; // Whether integral1 & integral2 have been computed.
	int der; // The number of derivatives to compute.
      
	/// Constructor.
	integralInfo() 
	{ integral1 = 0; integral2 = 0; integralset = false; der = -1; }

	/// Destructor.
	~integralInfo()
	{ erase(); }

	/// Resize/set the struct variables based in input agrguments.
	/// \param ider the new number of derivatives to compute.
	/// \param in1 number of coefficients in the u direction.
	/// \param in2 number of coefficients in the v direction.
	void resize(int ider, int in1, int in2)
	{
	    int ki, kj;
	    vec1.resize((ider+1)*in1*in1);
	    vec2.resize((ider+1)*in2*in2);
	    std::fill(vec1.begin(), vec1.end(), 0.0);
	    std::fill(vec2.begin(), vec2.end(), 0.0);

	    integral1 = new double**[ider+1];
	    integral2 = new double**[ider+1];

	    for (ki=0; ki<=ider; ki++)
		{
		    integral1[ki] = new double*[in1];
		    integral2[ki] = new double*[in2];

		    for (kj=0; kj<in1; kj++)
			integral1[ki][kj] = &vec1[(ki*in1+kj)*in1];

		    for (kj=0; kj<in2; kj++)
			integral2[ki][kj] = &vec2[(ki*in2+kj)*in2];
		}

	    integralset = 0;
	    der = ider;
	}

	/// Free the memory of the arrays in the struct.
	void erase()
	{
	    int ki;
	    for (ki=0; ki<=der; ki++)
		{
		    delete [] integral1[ki];
		    delete [] integral2[ki];
		}
	    delete [] integral1;
	    delete [] integral2;
	    integralset = false;
	}
    } integralInfo;

    std::vector<integralInfo> surf_integral_; // Integral calculations for the surfaces.

    /// Remove known coefs from the constraints. If there are no unknown coefs left
    /// in a constraint it is removed.
    /// \param constraints the new linear side constraints.
    /// \param surfs
    void removeKnownCoefs(std::vector<sideConstraintSet>& constraints) const;

};

}

#endif //
