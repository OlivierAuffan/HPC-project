#include <math.h>
#include <shalw.h>

void gauss_init(void) {
  double gmx, gmy, gsx, gsy;

  gmx = size_x * dx / 2 ;
  gmy = size_y * dy / 2 ;
  gsx = 25000 ;
  gsy = 25000 ;

  int start_x = start_band_x;
  int start_y = start_band_y - 1; // one extra line en top
  int end_x   = end_band_x;
  int end_y   = end_band_y + 1; // one extra line en bottom

  printf("-------------------------- %d %d\n", size_x, size_y);
  
  for (int i = start_x; i < end_x;  i++) {
    for (int j = start_y; j < end_y; j++) {
	if (i < 0 || j < 0 || i > size_x || j > size_y)
	    continue;
	HFIL(0, i - start_x, j - start_y) = height *
	    (exp(- pow((i * dx - gmx) / gsx, 2) / 2.)) *
	    (exp(- pow((j * dy - gmy) / gsy, 2) / 2.)) ;
    }
  }
}
