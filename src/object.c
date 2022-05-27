#include "../inc/object.h"
#include "../inc/private/objectimpl.h"

PetscErrorCode FcObjectGetComm(FcObject obj, MPI_Comm *comm) {
    FcObjectVerifyValidity(obj);

    *comm = obj->comm;

    return 0;
}

void FcObjectCreate_Private(FcObject obj, MPI_Comm comm, const char *class) {
    obj->comm = comm;
    obj->class = class;
    obj->name = "";
    obj->magic = FC_OBJECT_MAGIC_NUMBER;
}
