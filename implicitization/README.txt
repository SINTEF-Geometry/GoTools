Here is a description of my favourite apps. They are located on
...../gotools/implicitization/app.

test_implicit:

Takes a g2-file for a curve and an int (= degree) as input. The
program makes the implicit function of the given degree, and spits out
numbers that can be redirected to a file and gnusplotted (yes,
splotted!). For best results in Gnuplot, use 'set contour both', 'set
grid' and 'splot '<filename>' w l'.

test_find_self_intersections:

Takes a g2-file for a curve as input. The program does everything and
spits out gnuplottable numbers - the coordinates of the self-
intersections/cusps in the parameter plane. To plot in Gnuplot
together with the curve (see below), use 'plot '<curvefile>' w l,
'<selfintfile>' w p ps 3' (don't ask!).

test_implicit_for_surface:

Same as test_implicit for surfaces. The output is the implicit
function on the slice z=0.5. This z-value is hardcoded, but it should
be easy to change. Check out the q(x,y,0.5)=0 contour lines - they
should coincide with the intersection of the surface with the z=0.5
plane.

test_find_self_intersect_surface:

Same as test_self_intersections for surfaces. The output is a set of
points that lie on the self-intersection curves. In Gnuplot, use 'plot
'<filename>' for just the set of points, and 'plot '<filename>' w l'
for a plot with 3D-coinciding points connected two-by-two.


In addition there are two handy apps in ...../gotools/geometry/app:

gnuplot_curve:

Takes a g2-file as input and spits out gnuplottable numbers.

gnuplot_surface:

Ditto for surfaces.
