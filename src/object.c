#include "../inc/object.h"
#include "../inc/private/objectimpl.h"

PetscErrorCode FcObjectGetComm(FcObject obj, MPI_Comm *comm) {
    FcObjectVerifyValidity(obj);

    *comm = obj->comm;

    return 0;
}

void FcObjectInit(FcObject obj, MPI_Comm comm, const char *class) {
    obj->comm = comm;
    obj->class = class;
    obj->type = NULL;
    obj->name = NULL;
    obj->magic = FC_OBJECT_MAGIC_NUMBER;
}

PetscErrorCode FcObjectGetReference(FcObject obj, FcObject *ref) {
    FcObjectVerifyValidity(obj);

    *ref = obj;
    obj->refcnt++;

    return 0;
}

PetscErrorCode FcObjectRestoreReference(FcObject *ref) {
    FcObjectVerifyValidity(*ref);

    (*ref)->refcnt--;
    if ((*ref)->refcnt == 0)
        PetscCall((*ref)->ops.destroy(ref));
    *ref = NULL;

    return 0;
}
