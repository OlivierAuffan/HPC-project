#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "export.h"
#include "forward.h"
#include "shalw.h"

double hFil_forward(int t, int i, int j)
{
    // Phase d'initialisation du filtre
    // HPHY(t - 1, i, j) est encore nul
    if (t <= 2)
	return HPHY(t, i, j);
    return HPHY(t - 1, i, j) +
	alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

double uFil_forward(int t, int i, int j)
{
    // Phase d'initialisation du filtre
    // UPHY(t - 1, i, j) est encore nul
    if (t <= 2)
	return UPHY(t, i, j);
    return UPHY(t - 1, i, j) +
	alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

double vFil_forward(int t, int i, int j)
{
    // Phase d'initialisation du filtre
    // VPHY(t - 1, i, j) est encore nul
    if (t <= 2)
	return VPHY(t, i, j);
    return VPHY(t - 1, i, j) +
	alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
}

double hPhy_forward(int t, int i, int j)
{
    double c, d;

    c = 0.;
    if (i > 0)
	c = UPHY(t - 1, i - 1, j);

    d = 0.;
    if (j < size_y - 1)
	d = VPHY(t - 1, i, j + 1);

    return HFIL(t - 1, i, j) -
	dt * hmoy *
	((UPHY(t - 1, i, j) - c) / dx + (d - VPHY(t - 1, i, j)) / dy);
}

double uPhy_forward(int t, int i, int j)
{
    double b, e, f, g;

    if (i == size_x - 1)
	return 0.;

    b = 0.;
    if (i < size_x - 1)
	b = HPHY(t - 1, i + 1, j);

    e = 0.;
    if (j < size_y - 1)
	e = VPHY(t - 1, i, j + 1);

    f = 0.;
    if (i < size_x - 1)
	f = VPHY(t - 1, i + 1, j);

    g = 0.;
    if (i < size_x - 1 && j < size_y - 1)
	g = VPHY(t - 1, i + 1, j + 1);

    return UFIL(t - 1, i, j) +
	    dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) +
		  (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) -
		  (dissip * UFIL(t - 1, i, j)));
}

double vPhy_forward(int t, int i, int j)
{
    double c, d, e, f;

    if (j == 0)
	return 0.;

    c = 0.;
    if (j > 0)
	c = HPHY(t - 1, i, j - 1);

    d = 0.;
    if (i > 0 && j > 0)
	d = UPHY(t - 1, i - 1, j - 1);

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

void FORWARD(int t, int i, int j)
{
    HPHY(t, i, j) = hPhy_forward(t, i, j);
    HFIL(t, i, j) = hFil_forward(t, i, j);

    UPHY(t, i, j) = uPhy_forward(t, i, j);
    UFIL(t, i, j) = uFil_forward(t, i, j);

    VPHY(t, i, j) = vPhy_forward(t, i, j);
    VFIL(t, i, j) = vFil_forward(t, i, j);
}

// Cette fonction est longue et concentre les versions bandes/blocks
// sync/async pour factoriser le code et bien comprendre ce que ces modes
// impliquent sur une version standard. Ne pas hésiter à masquer des boucles
// ou des branchements dans un éditeur de code pour avoir un aperçu global
// de la fonction.
void forward(void)
{
    double svdt		 = 0.;
    int	t		 = 0;
    // double total_msg = 0, total_calc = 0, total_export = 0;

    // Special type MPI to exchange array columns
    MPI_Datatype column;
    MPI_Type_vector(band_size_y, 1, band_size_x + 2, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);

    // Requests for async mode
    MPI_Request r[3], s[3];
    for (int i = 0; i < 3; i++) {
	r[i] = MPI_REQUEST_NULL;
	s[i] = MPI_REQUEST_NULL;
    }

    if (file_export)
	create_file();
    
    // t = 0 is the initial state already in memory by gauss_init
    for (t = 1; t < nb_steps; t++) {// we talk about t iterations but actually we
	// calculate only t - 1 times
	if (t == 1) {
	    svdt = dt;
	    dt   = 0;
	}
	if (t == 2)
	    dt = svdt / 2.;

	if (file_export) {
	    // clock_t start_export = clock();
	    export_step(t - 1); // t - 1 is ready to export
	    // total_export += TIME(start_export, clock());
	}
	// recouvrement par le calcul en async

	// ECHANGE DE LIGNES
	// à l'instant t, la grille t-1 est calculée, on peut donc échanger les
	// lignes de t - 1
	// clock_t start_msg = clock();
	if (t > 1) {
	    if (async) {
		// Echange id-1 <=> id
		if (id > 0) {
		    MPI_Irecv(&HPHY(t - 1, 0, 0), band_size_x, MPI_DOUBLE,
			      id - 1, 0, MPI_COMM_WORLD, r);
		    MPI_Irecv(&UPHY(t - 1, 0, 0), band_size_x, MPI_DOUBLE,
			      id - 1, 0, MPI_COMM_WORLD, r + 1);
		    MPI_Isend(&VPHY(t - 1, 0, 1), band_size_x, MPI_DOUBLE,
			      id - 1, 0, MPI_COMM_WORLD, s);
		}
		// Echange id <=> id+1
		if (id < p - 1) {
		    MPI_Isend(&HPHY(t - 1, 0, band_size_y), band_size_x,
			      MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD, s + 1);
		    MPI_Isend(&UPHY(t - 1, 0, band_size_y), band_size_x,
			      MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD, s + 2);
		    MPI_Irecv(&VPHY(t - 1, 0, band_size_y + 1),
			      band_size_x, MPI_DOUBLE, id + 1, 0,
			      MPI_COMM_WORLD, r + 2);
		}
	    }
	    else {
		// Echange id-1 <=> id
		if (id > 0) {
		    MPI_Recv(&HPHY(t - 1, 0, 0), band_size_x, MPI_DOUBLE,
			     id - 1, 0, MPI_COMM_WORLD, NULL);
		    MPI_Recv(&UPHY(t - 1, 0, 0), band_size_x, MPI_DOUBLE,
			     id - 1, 0, MPI_COMM_WORLD, NULL);
		    MPI_Send(&VPHY(t - 1, 0, 1), band_size_x, MPI_DOUBLE,
			     id - 1, 0, MPI_COMM_WORLD);
		}
		
		// Echange id <=> id+1
		if (id < p - 1) {
		    MPI_Send(&HPHY(t - 1, 0, band_size_y), band_size_x,
			     MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD);
		    MPI_Send(&UPHY(t - 1, 0, band_size_y), band_size_x,
			     MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD);
		    MPI_Recv(&VPHY(t - 1, 0, band_size_y + 1),
			     band_size_x, MPI_DOUBLE, id + 1, 0,
			     MPI_COMM_WORLD, NULL);
		}
	    }
	}
	// total_msg += TIME(start_msg, clock());
	    
	// CALCULATIONS PREPARATION
	int start_x = 0;
	int start_y = 1; // skip first extra line
	int end_x   = band_size_x;
	int end_y   = band_size_y + 1; // skip last extra line
	    
	if (async) {// we have to reduce block by 1 line for each sides
	    start_y += 1; // no calculations for first line
	    end_y -= 1; // no calculations for last line
	}
	    
	// HERE ARE MOST CALCULATIONS for t
	// Peut facilement être parallélisé avec OpenMP
	// if async mode, messages are exchanged at the same time
	// clock_t start_calc = clock();
	#pragma omp parallel
	for (int y = start_y; y < end_y; y++)
	    for (int x = start_x; x < end_x; x++)
		FORWARD(t, x, y);
	// total_calc += TIME(start_calc, clock());
	    
	if (async) {
	    // Vérifier échange des bords t-1 avant de finir les calculs
	    // We need these receptions before finish calculations
	    // Should already be finished
	    // clock_t start_msg = clock();
	    MPI_Waitall(3, r, MPI_STATUSES_IGNORE); // Attente réception bords
	    // total_msg += TIME(start_msg, clock());
	    
	    // On fini les calculs des bords
	    start_x -= 1;
	    end_x += 1;
	    // clock_t start_calc = clock();
	    for (int x = start_x; x < end_x; x++) {
		FORWARD(t, x, 1);			 // first calculation line
		FORWARD(t, x, band_size_y); // last calculation line
	    }
	    // total_calc += TIME(start_calc, clock());
	    
	    // No need to wait before for these before calculations
	    // start_msg = clock();
	    MPI_Waitall(3, s, MPI_STATUSES_IGNORE); // Attente envoi bords
	    // total_msg += TIME(start_msg, clock());
	    // All messages have been exchanged : we can start new ones
	}
	    
	if (t == 2)
	    dt = svdt;
	
    
	if (file_export) {
	    // clock_t start_export = clock();
	    export_step(t - 1); // final iteration ready to export
	    finalize_export();
	    // total_export += TIME(start_export, clock());
	}
    
	// ID0_(printf("	Message exchange : %.2f\n", total_msg))
	//     ID0_(printf("	Calculations : %.2f\n", total_calc))
	// if (file_export) {
	//     ID0_(printf("	Export : %.2f\n", total_export))
	// 	}
	
	MPI_Barrier(MPI_COMM_WORLD);
    }
}
