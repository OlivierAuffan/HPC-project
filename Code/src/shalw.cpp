#include <stdlib.h>
#include <shalw.h>
#include <parse_args.hpp>
#include <memory.h>
#include <init.h>
#include <forward.h>
#include <export.h>

double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
int size_x, size_y, nb_steps;
double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
bool file_export;
std::string export_path;

int id,p,w;
MPI_Status status;

int main(int argc, char **argv) {

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&id);

	parse_args(argc, argv);
	printf("Command line options parsed\n");

	if(size_x%p){
		if(!id){ /*si processus de rang 0*/
			fprintf(stderr,"le nombre de processus ne divise pas la taille\n");
			exit(EXIT_FAILURE);
		}
	}
	w=size_x/p;
  
	alloc();
	printf("Memory allocated\n");
 
	/*elle sert a initialiser l'image*/ 
	gauss_init();
	printf("State initialised\n");

	forward();
	printf("State computed\n");
  
	dealloc();
	printf("Memory freed\n");
  
	return EXIT_SUCCESS;
}
