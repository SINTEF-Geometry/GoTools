#include "GoTools/tesselator/LineStrip.h"

namespace Go
{


//===========================================================================
    LineStrip::LineStrip(int n)
    {
	resize(n);
    }


    LineStrip::~LineStrip()
    {
    }

    LineStrip* LineStrip::asLineStrip()
	{
	    return this;
	}

    void LineStrip::resize(int n)
    {
	vert_.resize(n*3);
	strip_.resize(n);
	param_.resize(n);
	for (int i = 0; i < n; ++i) {
	    strip_[i] = i;
	}
    }

    int LineStrip::atBoundary(int idx) 
    { 
	int n = numVertices();
	return (idx==0 || idx==n-1) ? 1 : 0; 
    }


} // namespace Go

