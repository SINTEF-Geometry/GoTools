#ifndef _LINESTRIP_H
#define _LINESTRIP_H

#include "GoTools/tesselator/GeneralMesh.h"
#include <vector>


namespace Go
{

/** LineStrip: Structure for storing values for a line strip, i.e.
result of curve tesselation
 */

class GO_API LineStrip : public GeneralMesh
{
public:
  /// Constructor given size of mesh
    LineStrip(int n = 200);
    /// Destrcutor
    virtual ~LineStrip();

    /// Change mesh size
    void resize(int n);

    /// Number of nodes
    virtual int numVertices()
    { return (int)vert_.size()/3; }

    virtual double* vertexArray() { return &vert_[0]; }
    virtual double* paramArray() { return &param_[0]; }
     virtual int atBoundary(int idx);
     /// Indices for each line strip
    unsigned int* stripArray() { return &strip_[0]; }
    /// Indices for triangles. Not used
    virtual unsigned int* triangleIndexArray() { return &strip_[0]; }

    /// Casting. Return object as line strip.
    virtual LineStrip* asLineStrip();

private:
    std::vector<double> vert_;
    std::vector<double> param_;
    std::vector<unsigned int> strip_;
};



} // namespace Go





#endif // end of #ifdef _LINESTRIP_H_
