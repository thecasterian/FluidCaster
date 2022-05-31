#ifndef VIEWER_PCGNS_H
#define VIEWER_PCGNS_H

#include <mpi.h>
#include <petscsystypes.h>
#include "viewer.h"

PetscErrorCode FcViewerPCGNSOpen(MPI_Comm comm, const char *filename, FcViewer *viewer);

#endif
