#include <stdlib.h>
#include <shalw.h>
#include <parse_args.hpp>
#include <memory.h>
#include <init.h>
#include <forward.h>
#include <export.h>
#include <mpi.h>

double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
int size_x, size_y, nb_steps;
double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
bool file_export;
std::string export_path;

int id, p;
MPI_Status status;

int band_size_x, band_size_y, band_size;
int start_band_x, start_band_y, end_band_x, end_band_y;
int size, buffer_size;

int main(int argc, char **argv) {

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&id);

	parse_args(argc, argv);
	printf("Command line options parsed\n");

	if(size_y % p){
	    printf("le nombre de processus ne divise pas la taille\n");
	    exit(EXIT_FAILURE);
	}

	band_size_x = size_x;
	band_size_y = size_y / p;
	band_size = band_size_x * band_size_y;
	
	start_band_x = 0;
	start_band_y = id * band_size_y;
	end_band_x = band_size_x;
	end_band_y = (id + 1) * band_size_y;

	size = size_x * size_y;
	
	alloc();
	printf("%d Memory allocated\n", id);
 
	// /*elle sert a initialiser l'image*/ 
	gauss_init();
	printf("%d State initialised\n", id);
	
	forward();
	printf("%d State computed\n", id);
  
	dealloc();
	printf("%d Memory freed\n", id);

	MPI_Finalize();
	
	return EXIT_SUCCESS;
}
