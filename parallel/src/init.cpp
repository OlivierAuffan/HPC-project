#include <math.h>
#include <shalw.h>

void gauss_init(void) {
    double gmx, gmy, gsx, gsy;
    int sx, sy, ex, ey;
  
    gmx = size_x * dx / 2 ;
    gmy = size_y * dy / 2 ;
    gsx = 25000 ;
    gsy = 25000 ;

    sx = start_band_x;
    sy = start_band_y - 1;
    ex = end_band_x;
    ey = end_band_y + 1;

    for (int i = sx; i < ex;  i++) {
	for (int j = sy; j < ey; j++) {
	    HFIL(0, i - sx, j - sy) = height *
		(exp(- pow((i * dx - gmx) / gsx, 2) / 2.)) *
		(exp(- pow((j * dy - gmy) / gsy, 2) / 2.)) ;
	}
    }
}
