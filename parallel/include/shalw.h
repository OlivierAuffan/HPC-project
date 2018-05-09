#include <string>
extern double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
extern int size_x, size_y, nb_steps,w;
extern double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
extern bool file_export;
extern std::string export_path;

extern int id, p;
extern int band_size_x, band_size_y, band_size;
extern int start_band_x, start_band_y, end_band_x, end_band_y;
extern int buffer_size;

#define PRINT(message) \
    if (id == 0) \
	printf("%s", message);

#define HFIL(t, i, j) hFil[ (i) +			\
			    (j) * size_y +		\
			    ((t)%2) * buffer_size ]
#define UFIL(t, i, j) uFil[ (i) +			\
			    (j) * size_y +		\
			    ((t)%2) * buffer_size ]
#define VFIL(t, i, j) vFil[ (i) +			\
			    (j) * size_y +		\
			    ((t)%2) * buffer_size ]
#define HPHY(t, i, j) hPhy[ (i) +			\
			    (j) * size_y +		\
			    ((t)%2) * buffer_size ]
#define UPHY(t, i, j) uPhy[ (i) +			\
			    (j) * size_y +		\
			    ((t)%2) * buffer_size ]
#define VPHY(t, i, j) vPhy[ (i) +			\
			    (j) * size_y +		\
			    ((t)%2) * buffer_size ]
