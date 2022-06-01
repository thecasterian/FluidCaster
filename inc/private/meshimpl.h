#ifndef MESH_IMPL_H
#define MESH_IMPL_H

#include <petscdm.h>
#include "../mesh.h"

struct _p_FcMesh {
    /* MPI communicator. */
    MPI_Comm comm;
    /* Dimension. */
    PetscInt dim;
    /* Boundary types. */
    FcMeshBoundaryType bx, by, bz;
    /* Global number of elements. */
    PetscInt mx, my, mz;
    /* Number of processes. */
    PetscInt px, py, pz;
    /* Array containing the number of elements in each process */
    const PetscInt *lx, *ly, *lz;
    /* DMDA. */
    DM da;
    /* DMDA used by linear solvers. */
    DM dau, dav, daw, dap;
    /* DMStag. */
    DM stag;
};

#endif
