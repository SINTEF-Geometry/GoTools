
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/SISL_code.h"
#include "GoTools/geometry/sisl_file_io.h"

#include <fstream>
#include <stdlib.h>
#include <stdio.h>


#ifndef FALSE
#define FALSE 0
#endif




int main(int argc, char** argv)
{
  if (argc!=3)
    {
      std::cerr << "usage: <infile name> <outfile name>\n";
      return -1;
    }

  FILE* from_file = fopen(argv[1], "r" );
  if (!from_file)
    {
      std::cout << "File does not exist!" << std::endl;
      return -1;
    }
  std::ofstream outfile(argv[2]);

  int nmb_objs = 1;
  // char string[80];       /* character to read $ which identifies a comment */
  //read_non_comment(from_file,string);
  //sscanf(string,"%d",&nmb_objs);
  int kstat = 0;
  int ki;
  for (ki = 0; ki < nmb_objs; ++ki)
    {
      SISLObject* ob=NULL;
      file_to_obj(from_file, &ob, &kstat);
      if (kstat<0)
	{
	  std::cerr << "file_to_obj returned :" << kstat << std::endl;
	  return -1;
	}
      if (ob && ob->iobj==2)
	{
	  Go::SplineSurface* gs=Go::SISLSurf2Go(ob->s1);
	  gs->writeStandardHeader(outfile);
	  gs->write (outfile);
	}
      else if (ob && ob->iobj==1)
	{
	  Go::SplineCurve* gc=Go::SISLCurve2Go(ob->c1);
	  gc->writeStandardHeader(outfile);
	  gc->write (outfile);
	}
      else
	{
	  std::cerr << "Nothing to write" << std::endl;	  
	}
      if (ob) delete ob;
    }

  if (from_file) fclose(from_file);
}
