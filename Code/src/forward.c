#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>
#include <mpi.h>

double hFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //HPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return HPHY(t, i, j);
  return HPHY(t - 1, i, j) +
    alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

double uFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //UPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return UPHY(t, i, j);
  return UPHY(t - 1, i, j) +
    alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

double vFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //VPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return VPHY(t, i, j);
  return VPHY(t - 1, i, j) +
    alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
}

double hPhy_forward(int t, int i, int j) {
  double c, d;
  
  c = 0.;
  if (i > 0)
    c = UPHY(t - 1, i - 1, j);

  d = 0.;
  if (j < size_y - 1)
    d = VPHY(t - 1, i, j + 1);

  return HFIL(t - 1, i, j) -
    dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx +
		 (d - VPHY(t - 1, i, j)) / dy);
}

double uPhy_forward(int t, int i, int j) {
  double b, e, f, g;
  
  if (i == w*(id+1) - 1)
    return 0.;

  b = 0.;
  if (i < w - 1)
    b = HPHY(t - 1, i + 1, j);

  e = 0.;
  if (j < size_y - 1)
    e = VPHY(t - 1, i, j + 1);

  f = 0.;
  if (i < w - 1)
    f = VPHY(t - 1, i + 1, j);

  g = 0.;
  if (i < w - 1 && j < size_y - 1)
    g = VPHY(t - 1, i + 1, j + 1);

  return UFIL(t - 1, i, j) +
    dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) +
	  (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) -
	  (dissip * UFIL(t - 1, i, j)));
}

double vPhy_forward(int t, int i, int j) {
  double c, d, e, f;

  if (j == 0)
    return 0.;

  c = 0.;
  if (j > 0)
    c = HPHY(t - 1, i, j - 1);

  d = 0.;
  if (i > 0 && j > 0)
    d = UPHY(t - 1, i -1, j -1);

  e = 0.;
  if (i > 0)
    e = UPHY(t - 1, i - 1, j);

  f = 0.;
  if (j > 0)
    f = UPHY(t - 1, i, j - 1);

  return VFIL(t - 1, i, j) +
    dt * ((-grav / dy) * (HPHY(t - 1, i, j) - c) -
	  (pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) -
	  (dissip * VFIL(t - 1, i, j)));
}

void forward(void) {
  FILE *file = NULL;
  double svdt = 0.;
  int t = 0;
  
  if (file_export) {
    file = create_file();
    export_step(file, t);
  }
    
  for (t = 1; t < nb_steps; t++) {
    if (t == 1) {
      svdt = dt;
      dt = 0;
    }
    if (t == 2){
	dt = svdt / 2.;
    }

<<<<<<< HEAD
    if (t > 1) {
	if (id > 0)
	    {
		MPI_Recv(&HPHY(t - 1, 0, 0), band_size_x, MPI_DOUBLE,
			 id - 1, 0, MPI_COMM_WORLD, NULL);
		MPI_Recv(&UPHY(t - 1, 0, 0), band_size_x, MPI_DOUBLE,
			 id - 1, 0, MPI_COMM_WORLD, NULL);
		MPI_Send(&VPHY(t - 1, 0, 1), band_size_x, MPI_DOUBLE,
			 id - 1, 0, MPI_COMM_WORLD);
	    }
	if (id < p - 1)
	    {
		MPI_Send(&HPHY(t - 1, 0, band_size_y), band_size_x,
			 MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD);
		MPI_Send(&UPHY(t - 1, 0, band_size_y), band_size_x,
			 MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD);
		MPI_Recv(&VPHY(t - 1, 0, band_size_y + 1),
			 band_size_x, MPI_DOUBLE, id + 1, 0,
			 MPI_COMM_WORLD, NULL);
	    }
    }

    int start_x = 0;
    int start_y = 1; // skip first extra line
    int end_x   = band_size_x;
    int end_y   = band_size_y + 1; // skip last extra line
    
    for (int j = start_y; j < end_y; j++) {
      for (int i = start_x; i < end_x; i++) {
=======
    for (int j = 0; j < size_y; j++) {
      for (int i = 0; i < w; i++) {
>>>>>>> 7af82cdf57ff94e6a3794d1e0ea690a812094785
	HPHY(t, i, j) = hPhy_forward(t, i, j);
	UPHY(t, i, j) = uPhy_forward(t, i, j);
	VPHY(t, i, j) = vPhy_forward(t, i, j);
	HFIL(t, i, j) = hFil_forward(t, i, j);
	UFIL(t, i, j) = uFil_forward(t, i, j);
	VFIL(t, i, j) = vFil_forward(t, i, j);
      }
    }

    if (file_export) {
      export_step(file, t);
    }
    
    if (t == 2) {
      dt = svdt;
    }
  }

  if (file_export) {
    finalize_export(file);
  }
}
