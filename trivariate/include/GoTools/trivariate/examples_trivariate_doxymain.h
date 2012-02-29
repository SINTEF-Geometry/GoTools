/**
\example coons_patch_volume_gen coons_patch_volume_gen.C
\verbatim
\endverbatim
This program demonstrates the use of the static function 
\em createCoonsPatch
in namespace \em CoonsPatchVolumeGen'
The function can create a new \em SplineVolume representing the coons patch of
six SplineSurfaces, the six faces of the volume.

\example createCoonsVolume createCoonsVolume.C
\verbatim
\endverbatim
The program creates a Coons volume and performs smoothing of this volume
keeping the boundary surfaces fixed.

\example linear_swept_volume linear_swept_volume.C
\verbatim
\endverbatim
This program demonstrates the use of the static function \em linearSweptVolume
in the class \em SweepVolumeCreator.
The function can generate a SplineVolume by sweeping a surface along a curve
or sweeping a curve along a surface. 
A sweeping point on the curve or the surface must be specified.
If the point lies on the the surface, the surface will be swept along the
curve. If the point lies on the the curve, the curve will be swept along the
surface. The curve and the surface must be such that it doesn't lead to
self-intersection.

\example loft_volume_creator loft_volume_creator.C
\verbatim
\endverbatim
This program demonstrates the use of a static function \em loftVolume
in namespace \em LoftVolumeCreator.
The function use lofting to create a new \em SplineVolume based on a set of
surfaces.  The surfaces are not changed during the lofting process.
The surfaces must lie in the same space.

\example rotational_swept_volume rotational_swept_volume.C
\verbatim
\endverbatim
This program demonstrates the use of the static function
\em rotationalSweptVolume in the class \em SweepVolumeCreator.
The function can generate a SplineVolume by rotating a surface around an axis.
The surface must be such that it doesn't lead to self-intersection.

 */
