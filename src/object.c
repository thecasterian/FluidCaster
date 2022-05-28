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
