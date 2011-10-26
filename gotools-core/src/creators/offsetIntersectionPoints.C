
#include "GoTools/creators/SmoothTransition.h"

#include <math.h>
#include "GoTools/geometry/LineCloud.h"

#include <vector>
#include <fstream>

using std::vector;
using std::shared_ptr;
using std::max;
using std::min;
using namespace Go;


void SmoothTransition::
offsetIntersectionPoints(std::vector<Point>& ep, std::vector<Point>& eq,
			 std::vector<Point>& eoffp, std::vector<Point>& eoffq,
			 Point& eparp, Point& eparq,
			 std::vector<Point>& espine, std::vector<Point>& egeobb1,
			 std::vector<Point>& egeobb2, std::vector<Point>& ecrtan1,
			 std::vector<Point>& ecrtan2, std::vector<Point>& egeop,
			 std::vector<Point>& egeoq, std::vector<double>& curv_radis)
    /*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the position, unit tangent, curvature vector
*              and radius of curvature at a point at the intersection
*              curve between two offset surfaces. The first and second
*              derivatives of the point in the two surfaces and their
*              offset surfaces are given as input. Corresponding data on
*              the boundary curves and cross derivative curves of the
*              surfaces is also computed.
*
*
* INPUT      : ep      - 0-2 order derivatives of first surface.
*                        The sequence is position, first derivative in first
*                        parameter direction, first derivative in second
*                        parameter direction, (2,0) derivative, (1,1)
*                        derivative, (0,2) derivative and normal. (21 numbers)
*                        Compatible with output of s1421
*              eq      - 0-2 order derivatives of second surface. Same
*                        sequence as ep.
*              eoffp   - 0-2 order derivatives of offset surface of first
*                        surface. Same sequence as ep.
*              eoffq   - 0-2 order derivatives of offset surface of second
*                        surface. Same sequence as ep.
*              eparp   - Parameter pair in first surface of point
*              eparq   - Parameter pair in second surface of point
*
* OUTPUT     : 
*              jstat  - status messages  
*                         = 1      : Curvature radius infinit
*                         = 0      : ok, curvature radius
*                         < 0      : error
*              espine - 3-D geometry description of the intersection . The
*                       array contains: position, unit tangent, curvature
*                       and radius of curvature. (A total of 10 numbers)
*                       A radius of curvature =-1, indicates that the radius
*                       of curvature is infinit.
*              egeobb1 - 3-D geometry description of the intersection in the
*                       first input surface. Same sequence as espine.
*              egeobb2 - 3-D geometry description of the intersection in the
*                       second input surface. Same sequence as espine.
*              ecrtan1 - 3-D geometry description of the cross-derivative at
*                       the intersection in the first input surface. Same
*                       sequence as espine, but curvature vector and radius
*                       is not computed. Curvature radius is given as -1.0.
*              ecrtan2 - 3-D geometry description of the cross-derivative at
*                       the intersection in the second input surface. Same
*                       sequence as espine, but curvature vector and radius
*                       is not computed. Curvature radius is given as -1.0.
*              egeop  - Description of the intersection in the parameter plane
*                       of the first surface. The array contains: position,
*                       unit tangent, curvature and radius of curvature.
*                       (A total of 7 numbers)
*              egeoq  - Description of the intersection in the parameter plane
*                       of the second surface. The array contains: position,
*                       unit tangent, curvature and radius of curvature.
*                       (A total of 7 numbers)
*
* METHOD     : First the most lineary independent selection of 3 vectors
*              from the derivative vectors are found. Then equation systems
*              are made to determine the derivatives of the parameter values
*              with respect to the parameter value not present in the
*              selection of the three vectors. Then the double derivatives
*              of the parameter direction are found. This information is used
*              for expressing the 3-D tangent, curvature and the radius
*              of curvature of the intersection point in question.
*              Corresponding values in input surfaces, tangent surfaces
*              and both parameter planes are also found. 
*
*
* REFERENCES : 
*-
* CALLS      : s6norm,s6scpr,sqrt,fabs,s6length,s6crss
*              
*
* WRITTEN BY : Tor Dokken, SI, Oslo , Norway, 30 May 1988
* REWISED BY : Vibeke Skytt, SI, Oslo, Norway, Nov 1989
*
* PORTED BY: Sverre Briseid, Sintef, Oslo, Norway, 18 Dec 2002
*********************************************************************
*/
{
#ifdef CREATORS_DEBUG
    std::ofstream of2("data/debug.g2");
    for (int ki = 1; ki < 3; ++ki) {
	vector<double> line_pts1;
	line_pts1.insert(line_pts1.end(), ep[0].begin(), ep[0].end());
	Point to_pt = ep[0] + ep[ki];
	line_pts1.insert(line_pts1.end(), to_pt.begin(), to_pt.end());
	line_pts1.insert(line_pts1.end(), eoffp[0].begin(), eoffp[0].end());
	to_pt = eoffp[0] + eoffp[ki];
	line_pts1.insert(line_pts1.end(), to_pt.begin(), to_pt.end());
	LineCloud line_cloud1(&line_pts1[0], 2);
	line_cloud1.writeStandardHeader(of2);
	line_cloud1.write(of2);
	vector<double> line_pts2;
	line_pts2.insert(line_pts2.end(), eq[0].begin(), eq[0].end());
	to_pt = eq[0] + eq[ki];
	line_pts2.insert(line_pts2.end(), to_pt.begin(), to_pt.end());
	line_pts2.insert(line_pts2.end(), eoffq[0].begin(), eoffq[0].end());
	to_pt = eoffq[0] + eoffq[ki];
	line_pts2.insert(line_pts2.end(), to_pt.begin(), to_pt.end());
	LineCloud line_cloud2(&line_pts2[0], 2);
	line_cloud2.writeStandardHeader(of2);
	line_cloud2.write(of2);
    }
#endif // CREATORS_DEBUG

    Point snpu; //double snpu[3];         /* Normalized version of eoffpu                */
    Point snpv; //double snpv[3];         /* Normalized version of eoffpv                */
    Point spn; //double spn[3];          /* Vector snpu x snpv                          */
    Point snqs; //double snqs[3];         /* Normalized version of eoffqs                */
    Point snqt; //double snqt[3];         /* Normalized version of eoffqt                */
    Point sqn; //double sqn[3];          /* Vector snqs x snqt                          */
    Point sright; //double sright[3];       /* Right hand side when finding s"             */
    Point sdc; //double sdc[3];          /* Derivative of curve by w                    */
    Point sddc; //double sddc[3];         /* Second derivative of curve by w             */
    Point sddsp; //double sddsp[3];        /* Second derivative of spine curve.           */
    vector<Point>::iterator sqs;            /* Pointer to first row of matrix, offset      */
    vector<Point>::iterator sqt;            /* Pointer to second row of matrix, offset     */
    vector<Point>::iterator spu;            /* Pointer to third row of matrix, offset      */
    vector<Point>::iterator spw;            /* Pointer to fourth row of matrix, offset     */
    vector<Point>::iterator spuu;           /* Pointer to renamed (2,0)-derivative, offset */
    vector<Point>::iterator spuw;           /* Pointer to renamed (1,1)-derivative, offset */
    vector<Point>::iterator spww;           /* Pointer to renamed (0,2)-derivative, offset */
    vector<Point>::iterator sqss;           /* Pointer to renamed (2,0)-derivative, offset */
    vector<Point>::iterator sqst;           /* Pointer to renamed (1,1)-derivative, offset */
    vector<Point>::iterator sqtt;           /* Pointer to renamed (0,2)-derivative, offset */
    vector<Point>::iterator s2qs;           /* Pointer to renamed (1,0)-der, 1. surface    */
    vector<Point>::iterator s2qt;           /* Pointer to renamed (0,1)-der, 1. surface    */
    vector<Point>::iterator s2pu;           /* Pointer to renamed (1,0)-der, 2. surface    */
    vector<Point>::iterator s2pw;           /* Pointer to renamed (0,1)-der, 2. surface    */
    vector<Point>::iterator s2puu;          /* Pointer to renamed (2,0)-derivative         */
    vector<Point>::iterator s2puw;          /* Pointer to renamed (1,1)-derivative         */
    vector<Point>::iterator s2pww;          /* Pointer to renamed (0,2)-derivative         */
    vector<Point>::iterator s2qss;          /* Pointer to renamed (2,0)-derivative         */
    vector<Point>::iterator s2qst;          /* Pointer to renamed (1,1)-derivative         */
    vector<Point>::iterator s2qtt;          /* Pointer to renamed (0,2)-derivative         */

    vector<Point>::iterator sgeobb1;        /* Pointer to blend boundary curve of 1. surf  */
    vector<Point>::iterator sgeobb2;        /* Pointer to blend boundary curve of 2. surf  */
    vector<Point>::iterator scrtan1;        /* Pointer to cross tangent curve of 1. surf   */
    vector<Point>::iterator scrtan2;        /* Pointer to cross tangent curve of 2. surf   */
//     vector<Point> sgeobb1(3);        /* Pointer to blend boundary curve of 1. surf  */
//     vector<Point> sgeobb2(3);        /* Pointer to blend boundary curve of 2. surf  */
//     vector<Point> scrtan1(3);        /* Pointer to cross tangent curve of 1. surf   */
//     vector<Point> scrtan2(3);        /* Pointer to cross tangent curve of 2. surf   */

    double tt;              /* Value of det(snpu,snpv,snqs)                */
    double ts;              /* Value of det(snpu,snpv,snqt)                */
    double tv;              /* Value of det(snqs,snqt,snpu)                */
    double tu;              /* Value of det(snqs,snqt,snpv)                */
    double tlpu;            /* Length of eoffpu                            */
    double tlpv;            /* Length of eoffpv                            */
    double tlqs;            /* Length of eoffqs                            */
    double tlqt;            /* Length of eoffqt                            */
    double tmax1;           /* Variable used for maximal value             */
    double tmax2;           /* Variable used for maximal value             */
    double tdum;            /* Dummy variable                              */
    double tdom;            /* The denominator in an equation              */
    double tds;             /* ds/dw                                       */
    double tdt;             /* dt/dw                                       */
    double tdu;             /* du/dw                                       */
    double tddu;            /* ddu/dwdw                                    */
    double tdds;            /* dds/dwdw                                    */
    double tddt;            /* ddt/dwdw                                    */
    double twds;            /* ds/dw after renaming variable second time   */
    double twdt;            /* dt/dw after renaming variable second time   */
    double twdds;           /* dds/dwdw after renaming variable second time*/
    double twddt;           /* ddt/dwdw after renaming variable second time*/
    double twdu;            /* du/dw after renaming variable second time   */
    double twdv;            /* dv/dw after renaming variable second time   */
    double twddu;           /* ddu/dwdw after renaming variable second time*/
    double twddv;           /* ddv/dwdw after renaming variable second time*/

    /* Variables used to compute the cross-derivatives.  */

    Point sdiff; //double sdiff[3];        /* Difference vector between spine and BB curve*/
    Point sdiff2; //double sdiff2[3];       /* Difference vector between derivatives of
    //			   spine and BB curve.                         */
    Point scrtan; //double scrtan[3];       /* Not normalized cross tangent vector.        */
    Point svec1, svec2; //double svec1[3],svec2[3]; /* Vectors used to compute derivative of
    //			     cross tangent curve.                      */
    double tleng,tl3;       /* Length of cross tangent vector, length**3   */
    double tdot;            /* Scalar product of scrtan and svec1.         */
    double tfac;            /* Sign of cross product for cross-tangents.   */


    if (curv_radis.size() != 7)
	curv_radis.resize(7);
 

    /* Get input values into output. */
    egeop[0] = eparp;
    egeoq[0] = eparq;
    for (int ki = 1; ki < 3; ++ki) {
	egeop[ki].setValue(0.0, 0.0);
	egeoq[ki].setValue(0.0, 0.0);
    }
    curv_radis[5] = 0.0;
    curv_radis[6] = 0.0;

    egeobb1[0] = ep[0];
    egeobb2[0] = eq[0];

    /* Make position of intersection */

    espine[0] = 0.5*(eoffp[0] + eoffq[0]);

    espine[1].setValue(0.0, 0.0, 0.0);
    espine[2].setValue(0.0, 0.0, 0.0);
    curv_radis[0] = 0.0;

    /* Normalize derivative vectors */
    tlpu = eoffp[1].length();
    snpu = eoffp[1];
    snpu.normalize();
    tlpv = eoffp[2].length();
    snpv = eoffp[2];
    snpv.normalize();
    tlqs = eoffq[1].length();
    snqs = eoffq[1];
    snqs.normalize();
    tlqt = eoffq[2].length();
    snqt = eoffq[2];
    snqt.normalize();

    /* Make normal vector for both derivative pairs */
    spn = snpu%snpv;
    sqn = snqs%snqt;

    /* Make four scalar product to decide which of the 3 vectors snpu, snpv, 
     *  snqs, snqt spans the 3-D space best. (Have the biggest determinant)
     *  Remember (axb)c = det(a,b,c). The naming convention is that the
     *  name of the variable not present on the left hand side is used for
     *  the naming of the determinants. The determinants tt, ts, tu and tv tells
     *  which direction is most linearly dependent on the other directions  */

    tt = fabs(spn*snqs);
    ts = fabs(spn*snqt);
    tv = fabs(sqn*snpu);
    tu = fabs(sqn*snpv);

    /* We want to use the parameter direction names s, t and u on the left
     *  hand side of the equation system, thus we want to express all derivatives
     *  of the curve as functions of a new parameter  value w, which is chosen
     *  to be the parameter direction with partial first derivative most
     *  lineary dependent of the other parameter directions */

    tmax1 = max(tt,ts);
    tmax2 = max(tv,tu);

    if (tmax1 > tmax2) {

	/*  The s or t variable should not be used on the left hand side of
	 *   the equation system                                                */

	if (ts>tt) {

	    /*      The s variable should not be used on the left hand side of the 
	     *       equation system
	     *       The renaming of variables is as follows s->w, t->u, u->s, v->t */

	    spu  = eoffq.begin() + 2;//+6;
	    spw  = eoffq.begin() + 1;//+3;
	    spuu = eoffq.begin() + 5;//+15;
	    spuw = eoffq.begin() + 4;//+12;
	    spww = eoffq.begin() + 3;//+9;
	    sqs  = eoffp.begin() + 1;//+3;
	    sqt  = eoffp.begin() + 2;//+6;
	    sqss = eoffp.begin() + 3;//+9;
	    sqst = eoffp.begin() + 4;//+12;
	    sqtt = eoffp.begin() + 5;//+15;

	    s2pu  = eq.begin() + 2;//+6;
	    s2pw  = eq.begin() + 1;//+3;
	    s2puu = eq.begin() + 5;//+15;
	    s2puw = eq.begin() + 4;//+12;
	    s2pww = eq.begin() + 3;//+9;
	    s2qs  = ep.begin() + 1;//+3;
	    s2qt  = ep.begin() + 2;//+6;
	    s2qss = ep.begin() + 3;//+9;
	    s2qst = ep.begin() + 4;//+12;
	    s2qtt = ep.begin() + 5;//+15;

	    /* Set pointers to output arrays.  */

	    sgeobb1 = egeobb2.begin();
	    sgeobb2 = egeobb1.begin();
	    scrtan1 = ecrtan2.begin();
	    scrtan2 = ecrtan1.begin();
	    tfac = 1.0;
	} else {

	    /*      The t variable should not be used on the left hand side of the 
	     *       equation system
	     *       The renaming of variables is as follows s->u, t->w, u->s, v->t */

	    spu  = eoffq.begin() + 1;//+3;
	    spw  = eoffq.begin() + 2;//+6;
	    spuu = eoffq.begin() + 3;//+9;
	    spuw = eoffq.begin() + 4;//+12;
	    spww = eoffq.begin() + 5;//+15;
	    sqs  = eoffp.begin() + 1;//+3;
	    sqt  = eoffp.begin() + 2;//+6;
	    sqss = eoffp.begin() + 3;//+9;
	    sqst = eoffp.begin() + 4;//+12;
	    sqtt = eoffp.begin() + 5;//+15;

	    s2pu  = eq.begin() + 1;//+3;
	    s2pw  = eq.begin() + 2;//+6;
	    s2puu = eq.begin() + 3;//+9;
	    s2puw = eq.begin() + 4;//+12;
	    s2pww = eq.begin() + 5;//+15;
	    s2qs  = ep.begin() + 1;//+3;
	    s2qt  = ep.begin() + 2;//+6;
	    s2qss = ep.begin() + 3;//+9;
	    s2qst = ep.begin() + 4;//+12;
	    s2qtt = ep.begin() + 5;//+15;	

	    /* Set pointers to output arrays.  */

	    sgeobb1 = egeobb2.begin();
	    sgeobb2 = egeobb1.begin();
	    scrtan1 = ecrtan2.begin();
	    scrtan2 = ecrtan1.begin();
	    tfac = 1.0;
	}
    } else {
	/*  The u or v variable should not be used on the left hand side of
	 *   the equation system                                                */
	if (tu>tv)
	    {

		/*      The u variable should not be used on the left hand side of the
		 *       equation system.
		 *       The renaming of variables is as follows s->s, t->t, u->w, v->u */

		spu  = eoffp.begin() + 2;//+6;
		spw  = eoffp.begin() + 1;//+3;
		spuu = eoffp.begin() + 5;//+15;
		spuw = eoffp.begin() + 4;//+12;
		spww = eoffp.begin() + 3;//+9;
		sqs  = eoffq.begin() + 1;//+3;
		sqt  = eoffq.begin() + 2;//+6;
		sqss = eoffq.begin() + 3;//+9;
		sqst = eoffq.begin() + 4;//+12;
		sqtt = eoffq.begin() + 5;//15;
	
		s2pu  = ep.begin() + 2;//+6;
		s2pw  = ep.begin() + 1;//+3;
		s2puu = ep.begin() + 5;//+15;
		s2puw = ep.begin() + 4;//+12;
		s2pww = ep.begin() + 3;//+9;
		s2qs  = eq.begin() + 1;//+3;
		s2qt  = eq.begin() + 2;//+6;
		s2qss = eq.begin() + 3;//+9;
		s2qst = eq.begin() + 4;//+12;
		s2qtt = eq.begin() + 5;//+15;

		/* Set pointers to output arrays.  */

		sgeobb1 = egeobb1.begin();
		sgeobb2 = egeobb2.begin();
		scrtan1 = ecrtan1.begin();
		scrtan2 = ecrtan2.begin();
		tfac = 1.0;	
	    } else {
		/*      The v variable should not be used on the left hand side of the
		 *       equation system.
		 *       The renaming of variables is as follows s->s, t->t, u->u, v->w */

		spu  = eoffp.begin() + 1;//+3;
		spw  = eoffp.begin() + 2;//+6;
		spuu = eoffp.begin() + 3;//+9;
		spuw = eoffp.begin() + 4;//+12;
		spww = eoffp.begin() + 5;//+15;
		sqs  = eoffq.begin() + 1;//+3;
		sqt  = eoffq.begin() + 2;//+6;
		sqss = eoffq.begin() + 3;//+9;
		sqst = eoffq.begin() + 4;//+12;
		sqtt = eoffq.begin() + 5;//+15;

		s2pu  = ep.begin() + 1;//+3;
		s2pw  = ep.begin() + 2;//+6;
		s2puu = ep.begin() + 3;//+9;
		s2puw = ep.begin() + 4;//+12;
		s2pww = ep.begin() + 5;//+15;
		s2qs  = eq.begin() + 1;//+3;
		s2qt  = eq.begin() + 2;//+6;
		s2qss = eq.begin() + 3;//+9;
		s2qst = eq.begin() + 4;//+12;
		s2qtt = eq.begin() + 5;//+15;	

		/* Set pointers to output arrays.  */

		sgeobb1 = egeobb1.begin();
		sgeobb2 = egeobb2.begin();
		scrtan1 = ecrtan1.begin();
		scrtan2 = ecrtan2.begin();
		tfac = 1.0;
	    }
    }

    /* Now we can solve the equation systems for finding
     *  ds/dw, dt/dw and du/dw and afterwards for
     *  dds/(dwdw), ddt/(dwdw) and ddu/(dwdw), using Cramers Rule.
     *
     *  This equation system is derived in the following way:
     *  Our problem is defined as OP(u,w) - OQ(s,t) = 0. By taking the
     *  derivative of this equation with repsect to w, we get:
     *
     *  dOP(u,w) du   dOP(u,w)   dOQ(s,t) ds   dOQ(s,t) dt
     *  -------- -- + -------- - -------- -- - -------- -- = 0
     *  du       dw   dw         ds       dw   dt       dw    
     *
     *  By using a simplified notation this can be written:
     *
     *  OP u' + OP  - OQ s' - OQ t' = 0
     *    u       w     s       t
     *
     *  We can thus set up the equation system:
     *
     *                 s'
     *  (OQ   Q -OP ) (t') = OP
     *     s   t   u   u'      w
     * 
     *
     *
     *  By making one futher derivative we get an equation systen for s",t" and
     *  u".
     *
     *                   s"          2                          2   
     *  (OQ  OQ  - OP ) (t") = OP  u'  + 2OP  u' + OP   - OQ  s'  - 
     *     s   t     u   u"      uu         uw       ww     ss      
     *
     *                                           2
     *                         2OQ  s't' - OQ  t'
     *                            st         tt
     *
     *
     */

    /* Prepare normal vectors for determinants */
    spn = (*spu)%(*spw);
    sqn = (*sqs)%(*sqt);

    tdom = -(sqn*(*spu));

    //     if (DEQUAL(tdom,(double)0.0)) goto war101;
    if (tdom == 0.0)
	THROW("Not advisable to divide by 0.");

    /*  Lineary dependent vectors on left hand side if tdom = 0.0 */

    tds = -(spn*(*sqt)/tdom);
    tdt = spn*(*sqs)/tdom;
    tdu = sqn*(*spw)/tdom;

    sright = (*spuu)*tdu + 2.0*(*spuw)*tdu + (*spww) -
	((*sqss)*tds + (*sqst)*tdt)*tds - ((*sqtt)*tdt + (*sqst)*tds)*tdt;

    /* Calculate second derivatives of parameter direction with respect to
     *  the w-direction */

    tddu = sqn*sright/tdom;

    /* Use sqn for temporary storage of cross products */

    sqn = sright%(*sqt);
    tdds = -(sqn*(*spu)/tdom);
    sqn = (*sqs)%sright;
    tddt = -(sqn*(*spu)/tdom);

    /* We will now express the intersection curve locally as a function
     *  of the w-parameter.
     *
     *  s(w) = OP(u(w),w)
     *
     *  This gives the derivative
     *                         
     *   
     *  s' = OP u' + OP                                                   
     *         u       w
     *
     *  And the second derivative
     *
     *            2
     *  s" = OP  u'  + 2OP  u' + OP   + OP u"
     *         uu         uw       ww     u
     *
     *  The curvature vector is defined as the derivative of the unit tangent
     *  vector with respect to the arc length a(w):
     *
     *         d         d    c'(w)    dw   d    c'(w)      da
     *  k(a) = -- T(a) = -- ---------- -- = -- ---------- / --
     *         da        dw sqrt(c'c') da   dw sqrt(c'c')   dw
     *
     *
     *         d       c'(w)                c"        c' (c'c'')
     *         -- ----------------- =   ---------- - ------------- 
     *         dw sqrt(c'(w) c'(w))     sqrt(c'c')   sqrt(c'c')**3
     *
     *
     *
     *         da
     *         -- = sqrt(c'c')
     *         dw 
     */


    /* Compute position, unit tangent, curvature vector and radius of
       curvature at the spine curve, i.e. the intersection curve between
       the offset surfaces.                                               */
 
    sdc = (*spu)*tdu + (*spw);
    sddsp = ((*spuu)*tdu + 2.0*(*spuw))*tdu + (*spu)*tddu + (*spww);

    /* To simplify futher calculations we want to normalize the tangent vector
     *  and correspondingly divide the second derivative by the tangent length
     */

    tlpu = sdc.length();
    espine[1] = sdc;
    espine[1].normalize();

    // debug
    vector<double> line_pts1;
    line_pts1.insert(line_pts1.end(), espine[0].begin(), espine[0].end());
    Point to_pt = espine[0] + espine[1];
    line_pts1.insert(line_pts1.end(), to_pt.begin(), to_pt.end());
    LineCloud line_cloud1(&line_pts1[0], 1);
    std::ofstream of3("data/debug.g2");
    line_cloud1.writeStandardHeader(of3);
    line_cloud1.write(of3);
     // end debug

    //     if (DEQUAL(tlpu,(double)0.0)) goto war101;
    if (tlpu == 0.0)
	THROW("avoid dividing by 0.");

    sddsp *= 1/tlpu;

    /* Make curvature vector */
 
    tdum = sddsp*espine[1];

    espine[2] = (sddsp - espine[1]*tdum)/tlpu;

    /* Make 3-D radius of curvature */

    tdum = espine[2].length();

    //     if (DNEQUAL(tdum,(double)0.0))
    if (tdum != 0.0) {
	// 	    espine[9] = (double)1.0/tdum;
	curv_radis[0] = 1.0/tdum;
    }else {
	// 	    espine[9] = (double)-1.0;
	curv_radis[0] = -1.0;
	MESSAGE("Curvature radius is -1.0. This should not happen!");
    }

    /* The blend boundary curve of the first surface can be expressed
       locally as

       c1(w) = P(u(w),w).
      
       Differentiate and compute the position, unit tangent, curvature 
       and radius of curvature of this curve.                          */


    sdc = (*s2pu)*tdu + (*s2pw);
    sddc = ((*s2puu)*tdu + 2.0*(*s2puw))*tdu + (*s2pu)*tddu + (*s2pww);

    /* To simplify futher calculations we want to normalize the tangent vector
     *  and correspondingly divide the second derivative by the tangent length
     */

    tlpu = sdc.length();
    sgeobb1[1] = sdc;
    sgeobb1[1].normalize();

    //     if (DEQUAL(tlpu,(double)0.0)) goto war101;
    ALWAYS_ERROR_IF(tlpu == 0.0, "Trying to divide by 0.");

    sddc *= 1/tlpu;

    /* Make curvature vector */
 
    tdum = sddc*sgeobb1[1];
    sgeobb1[2] = (sddc - sgeobb1[1]*tdum)/tlpu;

    /* Make 3-D radius of curvature */

    tdum = sgeobb1[2].length();

    //     if (DNEQUAL(tdum,(double)0.0))
    if (tdum != 0.0) {
	// 	    sgeobb1[9] = (double)1.0/tdum;
	curv_radis[1] = 1/tdum;
    }else {
	// 	    sgeobb1[9] = (double)-1.0;
	curv_radis[1] = -1.0;
	MESSAGE("Curvature radius is -1.0. This should not happen!");
    }

    /* Compute cross-derivative tangent at the BB-curve in the first input
       surface. Let

       t(w) = (c1(w)-s(w)) x s'(w)

       Normalize t(w) to get the cross derivative tangent.                */

    /* Compute directional vector of cross-derivative tangent. */

    sdiff = sgeobb1[0] - espine[0];
    scrtan = sdiff%espine[1];

    /* Find the right sign of the cross-product. */
    sdiff2 = sgeobb2[0] - sgeobb1[0];
    if (scrtan*sdiff2 < 0.0)
	tfac = -1.0*tfac;
    else
	tfac =  1.0*tfac;

    /* Get length of directional vector.  */

    tleng = scrtan.length();
    tl3 = tleng*tleng*tleng;

    /* Help variables used to compute the derivative of the cross
       tangent curve.                                                */

    sdiff2 = sgeobb1[1] - espine[1];
    svec1 = sdiff%sddsp;
    svec2 = sdiff2%espine[1];

    svec1 += svec2;

    tdot = scrtan*svec1;

    /* Compute position and derivative of cross tangent, curvature vector
       is not calculated.                                                  */

    scrtan1[0] = tfac*scrtan/tleng;
    scrtan1[1] = tfac*(svec1/tleng - scrtan*tdot/tl3);
    scrtan1[2] *= 0.0; //*scrtan1_normal = 0.0;


    /* The radius of curvature is not calculated. */

    //     scrtan1[9] = -1.0;
    curv_radis[3] = 0.0;

    /* The blend boundary curve of the second surface can be expressed
       locally as

       c2(w) = Q(s(w),t(w)).
      
       Differentiate and compute the position, unit tangent, curvature 
       and radius of curvature of this curve.                          */


    sdc = (*s2qs)*tds + (*s2qt)*tdt;
    sddc = ((*s2qss)*tds + 2.0*(*s2qst)*tdt)*tds + (*s2qs)*tdds + (*s2qtt)*tdt*tdt + (*s2qt)*tddt;

    /* To simplify futher calculations we want to normalize the tangent vector
     *  and correspondingly divide the second derivative by the tangent length
     */

    tlqs = sdc.length();
    sgeobb2[1] = sdc;
    sgeobb2[1].normalize();

    //     if (DEQUAL(tlqs,(double)0.0)) goto war101;
    if (tlqs == 0.0)
	THROW("Trying to divide by 0.");

    sddc *= 1/tlqs;

    /* Make curvature vector */
 
    tdum = sddc*sgeobb2[1];

    sgeobb2[2] = (sddc - sgeobb2[1]*tdum)/tlqs;

    /* Make 3-D radius of curvature */

    tdum = sgeobb2[2].length();

    //     if (DNEQUAL(tdum,(double)0.0))
    if (tdum != 0.0) {
	// 	    sgeobb2[9] = (double)1.0/tdum;
	// 	    *sgeobb2_normal = 1.0/tdum;
	curv_radis[2] = 1.0/tdum;
    } else {
	// 	    sgeobb2[9] = (double)-1.0;
	// 	    *sgeobb2_normal = -1.0;
	curv_radis[2] = -1.0;
	MESSAGE("Curvature radius is -1.0. This should not happen!");
    }

    /* Compute cross-derivative tangent at the BB-curve in the second input
       surface. Let

       t(w) = (c2(w) - s(w)) x s'(w)

       Normalize t(w) to get the cross derivative tangent.                */

    /* Compute directional vector of cross-derivative tangent. */

    /* Turn tangent direction. */
    tfac = -1.0*tfac;
    sdiff = sgeobb2[0] - espine[0];
    //     s6crss(sdiff,espine+3,scrtan);
    scrtan = sdiff%espine[1];

    /* Get length of directional vector.  */
    tleng = scrtan.length();
    tl3 = tleng*tleng*tleng;

    /* Help variables used to compute the derivative of the cross tangent. */

    sdiff2 = sgeobb2[1] - espine[1];
    svec1 = sdiff%sddsp;
    svec2 = sdiff2%espine[1];
    svec1 += svec2;

    tdot = scrtan*svec1;

    /* Compute position and derivative of cross tangent, curvature vector
       is not calculated.                                                  */

    scrtan2[0] = tfac*scrtan/tleng;
    scrtan2[1] = tfac*(svec1/tleng - scrtan*tdot/tl3);
    scrtan2[2].setValue(0.0, 0.0, 0.0);

    /* The radius of curvature is not calculated. */

    //     scrtan2[9] = -1.0;
    curv_radis[4] = 0.0;
    //     *scrtan2_normal = -1.0;

    
    /* TO CALCULATE UNIT TANGENT CURVATURE AND RADIUS OF CURVATURE OF THE
     *  INTERSECTION POINT IN THE PARAMETER PLANES OF THE TWO SURFACES, WE
     *  NOW WANT TO CALCULATE THE TRUE VALUES OF ds/dw, dt/dw, dds/dw,
     *  ddt/dwdw, du/dw, dv/dw, ddu/dwdw, ddv/dwdw, where w is the parameter
     *  direction we have chosen the other directions to be expressed in.
     *  THUS UNDO the changing of parameter directions */

    if (tmax1 > tmax2) {

	/* First and second row of the surface were originally interchanged
	 *  Thus change sequence of these back again */ 

	if (ts>tt) {
	    /*       We used the following renaming of variables:
	     *       s->w, t->u, u->s, v->t, now express the behavior in the parameter
	     *       plane with the original ordering */

	    twds  = 1.0;
	    twdt  = tdu;
	    twdds = 0.0;
	    twddt = tddu;
	    twdu  = tds;
	    twdv  = tdt;
	    twddu = tdds;
	    twddv = tddt;

	} else {

	    /*      We used the following renaming of variables:
	     *       s->u, t->w, u->s, v->t, now express the behavior in the parameter
	     *       plane with the original ordering */

	    /*      The renaming of variables is as follows s->v, t->u, u->s, v->t */

	    twds  = tdu;
	    twdt  = 1.0;
	    twdds = tddu;
	    twddt = 0.0;
	    twdu  = tds;
	    twdv  = tdt;
	    twddu = tdds;
	    twddv = tddt;

	}
    } else {
	/*  Keep the sequence of surfaces */
	if (tu>tv) {
	    /*      We used the following renaming of variables:
	     *       s->s, t->t, u->w, v->u, now express the behavior in the parameter
	     *       plane with the original ordering */

	    twds  = tds;
	    twdt  = tdt;
	    twdds = tdds;
	    twddt = tddt;
	    twdu  = 1.0;
	    twdv  = tdu;
	    twddu = 0.0;
	    twddv = tddu;

	} else {
	    /*      We used the following renaming of variables:
	     *       s->s, t->t, u->u, v->w, now express the behavior in the parameter
	     *       plane with the original ordering */

	    twds  = tds;
	    twdt  = tdt;
	    twdds = tdds;
	    twddt = tddt;
	    twdu  = tdu;
	    twdv  = 1.0;
	    twddu = tddu;
	    twddv = 0.0;

	}
    }

    /* Now the variable twds, twdt, twdu, twdv contains derivatives of the
     *  parameter directions with respect to the w-variable. Correspondingly
     *  the second derivatives with respect to w are contained in twdds, twddt,
     *  twddu and twddv.
     *
     *  THE UNIT TANGENT, CURVATURE VECTOR AND RADIUS OF CURVATURE CAN NOW
     *  BE CALCULATED IN BOTH PARAMETER PLANES                                */

    /* Make description of intersection curve in parameter plane of first patch
     */


    tdum = sqrt(twdu*twdu + twdv*twdv);
    //     if (DEQUAL(tdum,(double)0.0))
    if (tdum == 0.0) {
	egeop[1][0] = 0.0;
	egeop[1][1] = 0.0;
	egeop[2][0] = 0.0;
	egeop[2][1] = 0.0;
	// 	    egeop[6] = (double)0.0;
	curv_radis[5] = 0.0;
    } else {

	/* Make unit tangent        */

	egeop[1][0] = twdu/tdum;
	egeop[1][1] = twdv/tdum;
                
	/* Make curvature vector    */

	tdom = egeop[1][0]*twddu + egeop[1][1]*twddv;
	egeop[2][0] = (twddu/tdum - egeop[1][0]*tdom/tdum)/tdum;
	egeop[2][1] = (twddv/tdum - egeop[1][1]*tdom/tdum)/tdum;
    }

    /* Make radius of curvature in parameter plane 1 */
    tdum = egeop[2].length();

    //     if (DNEQUAL(tdum,(double)0.0))
    if (tdum != 0.0) {
	// 	    egeop[6] = 1.0/tdum;
	curv_radis[5] = 1.0/tdum;
    } else {
	// 	    egeop[6] = -1.0;
	curv_radis[5] = -1.0;
    }

    /* Make description of intersection curve in parameter plane of second patch
     */
    tdum = sqrt(twds*twds + twdt*twdt);
    //     if (DEQUAL(tdum,(double)0.0))
    if (tdum == 0.0) {
	egeoq[1][0] = 0.0;
	egeoq[1][1] = 0.0;
	egeoq[2][0] = 0.0;
	egeoq[2][1] = 0.0;
	// 	    egeoq[6] = (double)0.0;
	curv_radis[6] = 0.0;
    } else {

	/* Make unit tangent        */
	egeoq[1][0] = twds/tdum;
	egeoq[1][1] = twdt/tdum;

	/* Make curvature vector    */
	tdom = egeoq[1][0]*twdds + egeoq[1][1]*twddt;
	egeoq[2][0] = (twdds/tdum - egeoq[1][0]*tdom/tdum)/tdum;
	egeoq[2][1] = (twddt/tdum - egeoq[1][1]*tdom/tdum)/tdum;
    }

    /* Make radius of curvature in parameter plane 2 */
    tdum = egeoq[2].length();

    //     if (DNEQUAL(tdum,(double)0.0))
    if (tdum != 0.0) {
	// 	    egeoq[6] = 1.0/tdum;
	curv_radis[6] = 1.0/tdum;
    } else {
	// 	    egeoq[6] = -1.0;
	curv_radis[6] = -1.0;
    }
}
