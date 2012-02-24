/**
\page viewlib GoTools viewlib

goview is an application program in the module viewlib. It is a 
utility to support visualization of curves, surface, point clouds and
line clouds.

goview can read curves and surfaces from an IGES file or from the GoTools
internal file format, \beginlink \link streamable_doc g2 \endlink, 
and visualize those. Currently, goview is not able
to visualize volumes. In the volume case, 
it is recommended to pick the boundary 
surfaces corresponding to the volume and draw them.

The program uses Qt for representing the GUI and openGl for graphics. Curves and
surfaces are tessellated in the submodule tessellate in the module 
gotools-core.
According to their type, curves and surfaces are approximated by triangles or
line segments that are convenient for visualization using openGl.

The model in the viewer is manipulated using the mouse keys. The left one 
rotates the model, the middle one is for zooming and the right one for 
translation of the model.

goview has got a graphical user interface. It has a graphical window, a window
containing an object list and a number of pull down
menus:
\arg \c file commands
to read and write geometry and to close the current session and make the 
viewer ready to read a new geometry file. 
\arg \c  view options
to choose shaded or wire frame mode for visualization. A highlight modus may be
toggled and some focusing facilities exist. 
\arg \c  select Selection of entities can be done
using this menu, in the object list or by
using the control key in combination with the left mouse key. 
\arg \c  group grouping of objects. This menu is not really used
\arg \c  object this menu
offers the possibility to alter the resolution of curves and surfaces and to
enable/disable selected entities from the view. 
Some commands have a key pad short cut.

goview is a utility and not a product. This implies unfortunately that the
help functionality is not implemented. Visualization of trimmed surfaces can
have some flaws, but they are normally repaired or minimized by increasing
the resolution. 
*/
