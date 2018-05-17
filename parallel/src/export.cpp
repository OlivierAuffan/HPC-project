#include <mpi.h>
#include <stdio.h>

#include "shalw.h"

static MPI_File file;
static MPI_Request req;

void create_file(void)
{
    req = 0;

    char fname[256];
    sprintf(fname, "%s/shalw_%dx%d_T%d.sav", export_path.c_str(), size_x,
	    size_y, nb_steps);

    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
		  MPI_INFO_NULL, &file);
}

void export_step(int t)
{
    MPI_Offset disp =
	(t * size + start_band_y * band_size_x) * sizeof(double);
    MPI_File_set_view(file, disp, MPI_DOUBLE, MPI_DOUBLE, "native",
		      MPI_INFO_NULL);
    if (async) {
	if (req) MPI_Wait(&req, NULL);
	
	MPI_File_iwrite(file, (void*)&HFIL(t, 0, 1), band_size, MPI_DOUBLE,
			&req);
    }
    else
	MPI_File_write(file, (void*)&HFIL(t, 0, 1), band_size, MPI_DOUBLE,
		       MPI_STATUS_IGNORE);
}

void finalize_export(void)
{
    if (req) MPI_Wait(&req, NULL);
    
    MPI_File_close(&file);
}
