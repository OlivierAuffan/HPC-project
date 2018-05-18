#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "export.h"
#include "forward.h"
#include "forward_utils.h"
#include "shalw.h"

void forward_simd(int t, int i, int j)
{
    HPHY(t, i, j) = hPhy_forward_simd(t, i, j);
    HFIL(t, i, j) = hFil_forward_simd(t, i, j);

    UPHY(t, i, j) = uPhy_forward_simd(t, i, j);
    UFIL(t, i, j) = uFil_forward_simd(t, i, j);

    VPHY(t, i, j) = vPhy_forward_simd(t, i, j);
    VFIL(t, i, j) = vFil_forward_simd(t, i, j);
}

void forward_(int t, int i, int j)
{
    HPHY(t, i, j) = hPhy_forward(t, i, j);
    HFIL(t, i, j) = hFil_forward(t, i, j);

    UPHY(t, i, j) = uPhy_forward(t, i, j);
    UFIL(t, i, j) = uFil_forward(t, i, j);

    VPHY(t, i, j) = vPhy_forward(t, i, j);
    VFIL(t, i, j) = vFil_forward(t, i, j);
}

void mpi_async(int t, MPI_Request r[], MPI_Request s[]) {

    if (id > 0) {
	MPI_Irecv(&HPHY(t - 1, 0, 0), band_size_x, MPI_DOUBLE,
		  id - 1, 0, MPI_COMM_WORLD, r);
	MPI_Irecv(&UPHY(t - 1, 0, 0), band_size_x, MPI_DOUBLE,
		  id - 1, 0, MPI_COMM_WORLD, r + 1);
	MPI_Isend(&VPHY(t - 1, 0, 1), band_size_x, MPI_DOUBLE,
		  id - 1, 0, MPI_COMM_WORLD, s);
    }
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

void mpi_(int t) {
    
    if (id > 0) {
	MPI_Recv(&HPHY(t - 1, 0, 0), band_size_x, MPI_DOUBLE,
		 id - 1, 0, MPI_COMM_WORLD, NULL);
	MPI_Recv(&UPHY(t - 1, 0, 0), band_size_x, MPI_DOUBLE,
		 id - 1, 0, MPI_COMM_WORLD, NULL);
	MPI_Send(&VPHY(t - 1, 0, 1), band_size_x, MPI_DOUBLE,
		 id - 1, 0, MPI_COMM_WORLD);
    }
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

void forward(void)
{
    double svdt;
    int t, start_x, start_y, end_x, end_y;
    MPI_Request r[3], s[3];
        
    svdt    = 0.;
    t       = 0;
   
    for (int i = 0; i < 3; i++) {
	r[i] = MPI_REQUEST_NULL;
	s[i] = MPI_REQUEST_NULL;
    }
    
    if (file_export) {
	create_file();
	export_step(t);
    }
    
    for (t = 1; t < nb_steps; t++) {
	if (t == 1) {
	    svdt = dt;
	    dt   = 0;
	}

	if (t == 2) dt = svdt / 2.;

	if (file_export) export_step(t);
	
	if (t > 1) {
	    if (async) mpi_async(t, r, s);
	    else mpi_(t);
	}
	
	start_x = 0;
	start_y = 1;
	end_x   = band_size_x;
	end_y   = band_size_y + 1;
	
	if (async) {
	    start_y += 1;
	    end_y -= 1;
	}

	// TODO :
	// mesurer temps normal, async, hybrid, ...
	// commencer le rapport

	if (t <= 2) {
	    if (hybrid) {
#pragma omp for
		for (int y = start_y; y < end_y; y++)
		    for (int x = start_x; x < end_x; x++)
			forward_(t, x, y);
	    } else {
		for (int y = start_y; y < end_y; y++)
		    for (int x = start_x; x < end_x; x++)
			forward_(t, x, y);   
	    }
	    
	} else {

	    if (hybrid) {
#pragma omp for
		for (int y = start_y; y < end_y; y++)
		    for (int x = start_x; x < end_x; x++)
			if (x == 0 || x == size_x - 1 ||
			    y == 0 || y == size_y - 1)
			    forward_(t, x, y);
			else
			    forward_simd(t, x, y);				
	    } else {
		for (int y = start_y; y < end_y; y++)
		    for (int x = start_x; x < end_x; x++)
			if (x == 0 || x == size_x - 1 ||
			    y == 0 || y == size_y - 1)
			    forward_(t, x, y);
			else
			    forward_simd(t, x, y);				
	    }
	}
	
	if (async) {
	    MPI_Waitall(3, r, MPI_STATUSES_IGNORE);

	    start_x -= 1;
	    end_x += 1;
	    for (int x = start_x; x < end_x; x++) {
		forward_(t, x, 1);
		forward_(t, x, band_size_y);
	    }

	    MPI_Waitall(3, s, MPI_STATUSES_IGNORE);
	}
	
	if (t == 2) dt = svdt;
    }
    
    if (file_export) {
	export_step(t);
	finalize_export();
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
}
