/**
 * @file object.h
 * @brief Common FluidCaster object interface.
 */

#ifndef OBJECT_H
#define OBJECT_H

#include <mpi.h>
#include <petscsystypes.h>

/**
 * @brief Any object in FluidCaster.
 */
typedef struct _p_FcObject *FcObject;

/**
 * @brief Gets the MPI communicator of an object.
 *
 * @param obj FluidCaster object.
 * @param[out] comm MPI comminicator.
 */
PetscErrorCode FcObjectGetComm(FcObject obj, MPI_Comm *comm);

#endif
