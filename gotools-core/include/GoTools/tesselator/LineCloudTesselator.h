//===========================================================================
//                                                                           
// File: LineCloudTesselator.h                                             
//                                                                           
// Created: Tue Oct 29 09:48:52 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _LINECLOUDTESSELATOR_H
#define _LINECLOUDTESSELATOR_H


#include "GoTools/tesselator/Tesselator.h"
#include "GoTools/geometry/LineCloud.h"

namespace Go
{

/** Tesselate line cloud. Preparation for visualization
 */
class GO_API LineCloudTesselator : public Tesselator
{
public:
  /// Constructor
    LineCloudTesselator(const Go::LineCloud& lc)
	: orig_cloud_(lc), render_cloud_(lc), scale_(1.0)
    {}
  /// Destructor
    virtual ~LineCloudTesselator();
  
    /// Perform tesselation
    virtual void tesselate();

    /// Fetch result
    const Go::LineCloud& getRenderCloud()
    {
	return render_cloud_;
    }

    /// The difference between the given line cload and the tesselated version
    /// is the scale. Set scaling factor.
    void setScale(double scale)
    {
	if (scale != scale_) {
	    scale = scale_;
	    tesselate();
	}
    }
    /// Fetch scaling factor.
    double getScale()
    {
	return scale_;
    }

private:
    const Go::LineCloud& orig_cloud_;
    Go::LineCloud render_cloud_;
    double scale_;
};

} // namespace Go



#endif // _LINECLOUDTESSELATOR_H

