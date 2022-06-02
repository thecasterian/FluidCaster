#include "../inc/private/viewerimpl.h"

PetscErrorCode FcViewerCreate(MPI_Comm comm, FcViewerType type, void *data, FcViewer *viewer) {
    FcViewer v;

    /* Allocate memory for the viewer. */
    PetscCall(PetscNew(&v));

    /* Initialize. */
    FcObjectInit((FcObject)v, comm, "FcViewer");
    v->obj.type = type;
    v->data = data;
    v->mesh = NULL;
    v->sol = NULL;

    /* Get the reference. */
    PetscCall(FcObjectGetReference((FcObject)v, (FcObject *)viewer));

    return 0;
}

PetscErrorCode FcViewerDestroy(FcViewer *viewer) {
    PetscCall(FcObjectRestoreReference((FcObject *)viewer));

    return 0;
}
