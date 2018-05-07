#include <string>
extern double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
extern int size_x, size_y, nb_steps,w;
extern double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
extern bool file_export;
extern std::string export_path;

extern int id, p;
extern int band_size_x, band_size_y, band_size;
extern int start_band_x, start_band_y, end_band_x, end_band_y;

#define PRINT(message) \
    if (id == 0) \
	printf("%s", message);

#define HFIL(t, i, j) hFil[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define UFIL(t, i, j) uFil[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define VFIL(t, i, j) vFil[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define HPHY(t, i, j) hPhy[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define UPHY(t, i, j) uPhy[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define VPHY(t, i, j) vPhy[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
