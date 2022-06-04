#ifndef MESH_IMPL_H
#define MESH_IMPL_H

#include <petscdm.h>
#include "../mesh.h"
#include "objectimpl.h"

struct _p_FcMesh {
    struct _p_FcObject obj;

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
    /* Domain range. */
    PetscReal xmin, xmax, ymin, ymax, zmin, zmax;
    /* Coordinate of element faces. */
    PetscReal *xf, *yf, *zf;
    /* Coordinate of element centers. */
    PetscReal *xc, *yc, *zc;
    /* DMDA. */
    DM da;
    /* DMDA used by linear solvers. */
    DM dau, dav, daw, dap;
    /* DMStag. */
    DM stag;
};

#endif
