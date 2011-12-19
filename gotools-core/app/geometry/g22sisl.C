
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ClassType.h"
#include "GoTools/geometry/sisl_file_io.h"

#include <fstream>
#include <stdlib.h>
#include <stdio.h>


#include "GoTools/geometry/SISL_code.h"


int main(int argc, char** argv)
{
  if (argc != 3)
    {
      std::cout << "Usage; infile oufile" << std::endl;
      return -1;
    }

  std::ifstream infile(argv[1]);
  FILE* to_file = fopen(argv[2], "w" );

  int type, dummy, nmb_col_elem;
  std::vector<shared_ptr<Go::GeomObject> > objs;
  while (infile >> type)
    {
      infile >> dummy;
      infile >> dummy;
      infile >> nmb_col_elem;
      for (int ki = 0; ki < nmb_col_elem; ++ki)
	  infile >> dummy;

    if (type == 200)
      {
	shared_ptr<Go::SplineSurface> surf(new Go::SplineSurface());
	surf->read(infile);
	objs.push_back(surf);
      }
    else if (type == 100)
      {
	shared_ptr<Go::SplineCurve> crv(new Go::SplineCurve());
	crv->read(infile);
	objs.push_back(crv);
      }
    else
      {
	std::cout << "Expecting spline curves & surfaces only!" << std::endl;
      }
    }

  // We write to file the number of objects.
  fprintf(to_file,"$ Number of objects in file\n");
  fprintf(to_file,"%d\n\n",int(objs.size()));

  for (size_t ki = 0; ki < objs.size(); ++ki)
    {
      /* open the file with the name to_file */
      if (objs[ki]->instanceType() == Go::Class_SplineSurface)
	{
	  SISLSurf* srf=GoSurf2SISL
	    (*(dynamic_pointer_cast<Go::SplineSurface, Go::GeomObject>(objs[ki])));
	  surface_to_file(to_file, srf);
	}
      else if (objs[ki]->instanceType() == Go::Class_SplineCurve)
	{
	  SISLCurve* cv=Curve2SISL
	    (*(dynamic_pointer_cast<Go::SplineCurve, Go::GeomObject>(objs[ki])));
	  curve_to_file(to_file, cv);
	}
      fprintf(to_file,"\n");
    }

  fclose(to_file);
}
