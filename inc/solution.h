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
 * @brief Gets velocity vectors of a solution.
 *
 * @param sol Solution.
 * @param[out] u X-velocity, or NULL.
 * @param[out] v Y-velocity, or NULL.
 * @param[out] w Z-velocity, or NULL.
 *
 * @warning User must not destroy the vectors.
 */
PetscErrorCode FcSolutionGetVelocityVec(FcSolution sol, Vec *u, Vec *v, Vec *w);

/**
 * @brief Gets the true pressure vector.
 *
 * @param sol Solution.
 * @param p_true True pressure.
 *
 * @warning User must not destroy the vector.
 */
PetscErrorCode FcSolutionGetTruePressureVec(FcSolution sol, Vec *p_true);

/**
 * @brief Views a solution.
 *
 * @param sol Solution.
 * @param viewer Viewer.
 */
PetscErrorCode FcSolutionView(FcSolution sol, FcViewer viewer);

#endif
