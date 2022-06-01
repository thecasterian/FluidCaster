#ifndef OBJECT_IMPL_H
#define OBJECT_IMPL_H

#include <mpi.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include "../object.h"

typedef struct _p_FcObjectOps {
    PetscErrorCode (*destroy)(FcObject *);
} *FcObjectOps;

struct _p_FcObject {
    /* Magic number. */
    PetscInt magic;

    /* Operations. */
    struct _p_FcObjectOps ops;
    /* Reference count. */
    PetscInt refcnt;

    /* MPI communicator. */
    MPI_Comm comm;
    /* Class name. */
    const char *class;
    /* Type name (if the subclass exists). */
    const char *type;
    /* Object name. */
    const char *name;
};

/* Magic number to verify the object validity. */
#define FC_OBJECT_MAGIC_NUMBER 0xFC42

#define FcObjectVerifyValidity(obj) \
    do { \
        if (!(obj) || (obj)->magic != FC_OBJECT_MAGIC_NUMBER) \
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_CORRUPT, "Invalid object"); \
    } while (0)

void FcObjectInit(FcObject obj, MPI_Comm comm, const char *class);
PetscErrorCode FcObjectGetReference(FcObject obj, FcObject *ref);
PetscErrorCode FcObjectRestoreReference(FcObject *ref);

#endif
