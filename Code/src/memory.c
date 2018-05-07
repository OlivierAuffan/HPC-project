#include <stdlib.h>
#include <shalw.h>

void alloc(void) {

    int buffer_size = band_size_x * (band_size_y + 2);
    
    hFil = (double *) calloc(2 * buffer_size, sizeof(double));
    uFil = (double *) calloc(2 * buffer_size, sizeof(double));
    vFil = (double *) calloc(2 * buffer_size, sizeof(double));
    hPhy = (double *) calloc(2 * buffer_size, sizeof(double));
    uPhy = (double *) calloc(2 * buffer_size, sizeof(double));
    vPhy = (double *) calloc(2 * buffer_size, sizeof(double));
}

void dealloc(void) {
    free(hFil);
    free(uFil);
    free(vFil);
    free(hPhy);
    free(uPhy);
    free(vPhy);
}
