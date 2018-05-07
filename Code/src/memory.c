#include <stdlib.h>
#include <shalw.h>

void alloc(void) {
<<<<<<< HEAD

    int buffer_size = band_size_x * (band_size_y + 2);
    
    hFil = (double *) calloc(2 * buffer_size, sizeof(double));
    uFil = (double *) calloc(2 * buffer_size, sizeof(double));
    vFil = (double *) calloc(2 * buffer_size, sizeof(double));
    hPhy = (double *) calloc(2 * buffer_size, sizeof(double));
    uPhy = (double *) calloc(2 * buffer_size, sizeof(double));
    vPhy = (double *) calloc(2 * buffer_size, sizeof(double));
=======
	hFil = (double *) calloc(2*size_x*size_y, sizeof(double));
  	uFil = (double *) calloc(2*size_x*size_y, sizeof(double));
  	vFil = (double *) calloc(2*size_x*size_y, sizeof(double));
  	hPhy = (double *) calloc(2*size_x*size_y, sizeof(double));
  	uPhy = (double *) calloc(2*size_x*size_y, sizeof(double));
  	vPhy = (double *) calloc(2*size_x*size_y, sizeof(double));
>>>>>>> 7af82cdf57ff94e6a3794d1e0ea690a812094785
}

void dealloc(void) {
    free(hFil);
    free(uFil);
    free(vFil);
    free(hPhy);
    free(uPhy);
    free(vPhy);
}
