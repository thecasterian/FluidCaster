#include "../inc/private/viewerimpl.h"

PetscErrorCode FcViewerCreate(MPI_Comm comm, FcViewerType type, FcViewerOps ops, void *data, FcViewer *viewer) {
    /* Allocate memory for the viewer. */
    PetscCall(PetscNew(viewer));
    FcObjectInit((FcObject)(*viewer), comm, "FcViewer");
    (*viewer)->obj.type = type;
    (*viewer)->ops[0] = *ops;
    (*viewer)->data = data;

    /* Initialize. */
    (*viewer)->mesh = NULL;
    (*viewer)->sol = NULL;

    return 0;
}

PetscErrorCode FcViewerClose(FcViewer *viewer) {
    PetscCall((*viewer)->ops->close(viewer));

    return 0;
}
