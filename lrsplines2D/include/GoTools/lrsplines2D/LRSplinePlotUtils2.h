#ifndef LRSPLINEPLOT_UTILS2_H
#define LRSPLINEPLOT_UTILS2_H

#include "GoTools/lrsplines2D/Mesh2D.h"
//#include "BSplineFunction.h"

namespace Go
{

// NB: Plot functions below require locale to be properly set.
//     To do this, include <locale> in the main program 
//     compilation unit, and execute:
//     setlocale(LC_ALL, "en_US.UTF.8");
void plot_mesh(const Mesh2D& m, int thick_threshold = 2); // plots mesh to standard output
void plot_largest_tensorgrid(const Mesh2D& m); // extract/plot largest possible tensor grid within m.
void plot_supports_at_corner(const Mesh2D& m, int xpos, int ypos, int xdeg, int ydeg);
void plot_support(const Mesh2D& m, const int* const kvec1, const int* const kvec2, int len1, int len2);
void plot_all_supports(const Mesh2D& m, int x_deg, int y_deg);
void plot_rect_domain(const Mesh2D& m, int xmin, int ymin, int xmax, int ymax);
void plot_history(const Mesh2D& m); 
  //void plot_bspline_function(const Mesh2D& m, const BSplineFunction& b);

// Identify the largest tensorgrid found within a general Mesh (when run twice, once in each direction)
std::vector<int> orig_tensorgrid_knotpositions(const Mesh2D& m, Direction2D d);

}; // end namespace Go

#endif
