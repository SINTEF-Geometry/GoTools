//  #ifdef MICROSOFT
//  #include <io.h>
//  #else
//  #include <unistd.h>
//  #endif

//  #include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "GoTools/viewlib/raster.h"
//#include "endian.h"




#define CRIT_ERR(stmnt) \
  printf("In file %s, line %d: ", __FILE__, __LINE__), (stmnt), exit(0)






//
// 000202: Added this. Modifying the writing function from 'picgen.C'.
//

void skip_ppm_comments(FILE *f)
{
  char buf[1000], *s;
  
  // printf("[%c=%d]\n", ungetc(getc(f), f), ungetc(getc(f), f));
  //
  // 000301: Was it necessary to include the test on ASCII(10) because
  //         only 3 bytes were read for the magic tag? (Instead of 4?!)
  //
  while (((ungetc(getc(f), f))=='#') || ((ungetc(getc(f), f))==10))
    {
      s=fgets(buf, 1000, f);
      if (s!=buf)
	CRIT_ERR(puts("Something went wrong!"));
      // printf("Comment skipped: [%s]\n", buf);
    }
}

unsigned char *read_ppm_file(const char * const name,
			     int * const xs, int * const ys)
{
  FILE *f=fopen(name, "r");
  if (f==NULL)
    CRIT_ERR(printf("Couldn't open %s.\n", name));

  int ppm, pgm;
  {
    char buf[3];

    fgets(buf, 3, f);
    // printf("[%s]\n", buf);
    ppm=(strncmp(buf, "P6", 2)==0);
    pgm=(strncmp(buf, "P5", 2)==0);
    if ((!ppm) && (!pgm))
      CRIT_ERR(printf("Oops! This (%s) was neither a ppm-file nor "
		      "a pgm-file!\n", name));
  }

  //
  // 000301: Just realized that xv inserts a comment. There may be
  //         comments also in other places. This could cause tricky to
  //         find problems at a later time!
  //
  skip_ppm_comments(f);

  fscanf(f, "%d\n", xs);
  // printf("x=%d\n", *xs);
  fscanf(f, "%d\n", ys);
  // printf("y=%d\n", *ys);
  
  {
    int maxval;

    fscanf(f, "%d\n", &maxval);
    if (maxval!=255)
      CRIT_ERR(printf("Oops! Expected maxval=255, not %d!\n", maxval));
  }

  unsigned char *d;
  int size=0;
  
  if (ppm)
    size=*xs*(*ys)*3;
  else
    if (pgm)
      size=*xs*(*ys);

  d=new unsigned char[size];
  if (d==NULL)
    CRIT_ERR(printf("Couldn't allocate %d bytes for image data.\n", size));
  

  {
    int tmp;
    
    if ((tmp=(int)fread(d, 1, size, f)) != size)
      CRIT_ERR(printf("Couldn't read all the %d bytes, only %d!\n",
		      size, tmp));
  }
  
  fclose(f);
  return d;
}






//
// Writes out a ppm file.
//

void write_ppm_file(const char * const name,
		    const unsigned char * const d, const int xs, const int ys, bool rgb)
{
  FILE *f=fopen(name, "w");
  if (f==NULL)
    CRIT_ERR(printf("Couldn't open %s for writing.\n", name));

  if (rgb) {
      fputs("P6\n", f);
  } else {
      fputs("P5\n", f);
  }
  fprintf(f, "%d %d\n", xs, ys);
  fprintf(f, "%d\n", 255);

  if (rgb) {
      if ((int)fwrite(d, 1, xs*ys*3, f) != xs*ys*3) {
	  CRIT_ERR(printf("Couldn't write all %d bytes to %s.\n", xs*ys*3, name));
      }
  } else {
      if ((int)fwrite(d, 1, xs*ys, f) != xs*ys) {
	  CRIT_ERR(printf("Couldn't write all %d bytes to %s.\n", xs*ys, name));
      }
  }

  fclose(f);
}
