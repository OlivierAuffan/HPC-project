#include <mpi.h>
#include <stdio.h>

#include "shalw.h"

static MPI_File		fh;
static MPI_Request  request;

void create_file(void)
{
    request = 0;

    char fname[256];
    sprintf(fname, "%s/shalw_%dx%d_T%d.sav", export_path.c_str(), size_x,
	    size_y, nb_steps);

    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
		  MPI_INFO_NULL, &fh);
}

void export_step(int t)
{
    MPI_Offset disp =
	(t * size + start_band_y * band_size_x) * sizeof(double);
    MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native",
		      MPI_INFO_NULL);
    if (async) {
	if (request)
	    MPI_Wait(&request, NULL);
	MPI_File_iwrite(fh, (void*)&HFIL(t, 0, 1), band_size, MPI_DOUBLE,
			&request);
    }
    else
	MPI_File_write(fh, (void*)&HFIL(t, 0, 1), band_size, MPI_DOUBLE,
		       MPI_STATUS_IGNORE);
}

void finalize_export(void)
{
    if (request)
	MPI_Wait(&request, NULL);
    MPI_File_close(&fh);
}
