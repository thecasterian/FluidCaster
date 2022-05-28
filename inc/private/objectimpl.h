#ifndef OBJECT_IMPL_H
#define OBJECT_IMPL_H

#include <mpi.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include "../object.h"

struct _p_FcObject {
    /** MPI communicator. */
    MPI_Comm comm;
    /** Class name. */
    const char *class;
    /** Type name (if the subclass exists). */
    const char *type;
    /** Object name. */
    const char *name;
    /** Magic number. */
    PetscInt magic;
};

/** Magic number to verify the object validity. */
#define FC_OBJECT_MAGIC_NUMBER 0xFC42

/**
 * @brief Verifies the object validity.
 *
 * @param obj FluidCaster object.
 */
#define FcObjectVerifyValidity(obj) \
    do { \
        if (!obj || obj->magic != FC_OBJECT_MAGIC_NUMBER) \
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_CORRUPT, "Invalid object"); \
    } while (0)

/**
 * @brief Initializes an object.
 *
 * @param h Any FluidCaster object.
 * @param comm MPI communicator.
 * @param[in] class Class name.
 * @param[in] subclass Subclass name.
 */
void FcObjectInit(FcObject obj, MPI_Comm comm, const char *class);

#endif
