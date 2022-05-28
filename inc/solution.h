/**
 * @file solution.h
 * @brief Problem solutions.
 */

#ifndef SOLUTION_H
#define SOLUTION_H

#include <petscsystypes.h>
#include <petscvec.h>
#include "mesh.h"

/**
 * @brief Solution of a fluid dynamics problem.
 */
typedef struct _p_FcSolution *FcSolution;

/**
 * @brief Creates a solution.
 *
 * @param mesh Mesh to create the solution on.
 * @param[out] sol Resulting solution.
 */
PetscErrorCode FcSolutionCreate(FcMesh mesh, FcSolution *sol);

/**
 * @brief Destroys a solution.
 *
 * @param[in] sol Solution to destroy.
 */
PetscErrorCode FcSolutionDestroy(FcSolution *sol);

/**
 * @brief Gets the velocity and pressure vectors of a solution. Use NULL in place of an output parameter without
 * interest.
 *
 * @param sol Solution.
 * @param[out] u X-velocity.
 * @param[out] v Y-velocity.
 * @param[out] w Z-velocity.
 * @param[out] p Pressure vector.
 *
 * @warning User must not destroy the vectors.
 */
PetscErrorCode FcSolutionGetVelocityPressureVec(FcSolution sol, Vec *u, Vec *v, Vec *w, Vec *p);

#endif
