#ifndef VIEWER_IMPL_H
#define VIEWER_IMPL_H

#include <mpi.h>
#include "../mesh.h"
#include "../solution.h"
#include "../viewer.h"
#include "objectimpl.h"

typedef struct _FcViewerOps *FcViewerOps;
struct _FcViewerOps {
    PetscErrorCode (*close)(FcViewer *);
    PetscErrorCode (*viewmesh)(FcViewer, FcMesh);
    PetscErrorCode (*viewsol)(FcViewer, FcSolution);
};

struct _p_FcViewer {
    /** Object. */
    struct _p_FcObject obj;
    /** Operations. */
    struct _FcViewerOps ops[1];
    /** Type-specific data. */
    void *data;

    /** Mesh. */
    FcMesh mesh;
    /** Solution. */
    FcSolution sol;
};

PetscErrorCode FcViewerCreate(MPI_Comm comm, FcViewerType type, FcViewerOps ops, void *data, FcViewer *viewer);

#endif
