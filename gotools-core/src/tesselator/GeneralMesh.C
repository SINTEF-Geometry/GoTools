#include "GoTools/tesselator/GeneralMesh.h"

namespace Go
{


    GeneralMesh::~GeneralMesh()
    {
    }

    RegularMesh* GeneralMesh::asRegularMesh()
	{
	    return 0;
	}

    LineStrip* GeneralMesh::asLineStrip()
	{
	    return 0;
	}

    GenericTriMesh* GeneralMesh::asGenericTriMesh()
	{
	    return 0;
	}

} // namespace Go

