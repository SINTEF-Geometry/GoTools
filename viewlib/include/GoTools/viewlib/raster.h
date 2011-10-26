#ifndef RASTERUTILS_H_INCLUDED







//
// Added 000202.
//

/// Returns RGB data in a new'ed array.
unsigned char *read_ppm_file(const char * const name,
			     int * const xs, int * const ys);

/// Writes out a ppm file.
void write_ppm_file(const char * const name,
		    const unsigned char * const d, const int xs, const int ys,
		    bool rgb = true);






#define RASTERUTILS_H_INCLUDED
#endif
