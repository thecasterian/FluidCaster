/**
 * @file ns.h
 * @brief Navier-Stokes equation solver interface.
 */

#ifndef NS_H
#define NS_H

#include <mpi.h>
#include <petscsystypes.h>
#include "solution.h"

/**
 * @brief Navier-Stokes solver that solves a fluid dynamics problem.
 */
typedef struct _p_FcNS *FcNS;

/**
 * @brief Creates a Navier-Stokes solver.
 *
 * @param mesh Mesh where the problem is defined.
 * @param sol Solution of the Navier-Stokes equation.
 * @param mat Fluid definition.
 * @param[out] ns Resulting Navier-Stokes solver.
 */
PetscErrorCode FcNSCreate(FcMesh mesh, FcSolution sol, FcMaterial mat, FcNS *ns);

/**
 * @brief Destroys a Navier-Stokes solver.
 *
 * @param[in] ns Navier-Stokes solver to destroy.
 */
PetscErrorCode FcNSDestroy(FcNS *ns);

/**
 * @brief Sets the current time of a solver.
 *
 * @param ns Navier-Stokes solver.
 * @param t Time.
 */
PetscErrorCode FcNSSetTime(FcNS ns, PetscReal t);

/**
 * @brief Sets the time step size.
 *
 * @param ns Navier-Stokes solver.
 * @param timestep Initial time step size.
 */
PetscErrorCode FcNSSetTimeStep(FcNS ns, PetscReal timestep);

/**
 * @brief Sets the maximum number of time steps.
 * @note If the solver is not transient, this function has no effect.
 *
 * @param ns Navier-Stokes solver.
 * @param maxsteps Maximum number of time steps to solve.
 */
PetscErrorCode FcNSSetMaxSteps(FcNS ns, PetscInt maxsteps);

/**
 * @brief Sets parameters in a Navier-Stokes solver from the options database.
 *
 * @param ns Navier-Stokes solver to set options for.
 */
PetscErrorCode FcNSSetFromOptions(FcNS ns);

/**
 * @brief Sets up the data structures inside a Navier-Stokes solver.
 *
 * @param ns Navier-Stokes solver to set up.
 */
PetscErrorCode FcNSSetUp(FcNS ns);

/**
 * @brief Gets the current time.
 *
 * @param ns Navier-Stokes solver.
 * @param[out] t Current time.
 */
PetscErrorCode FcNSGetTime(FcNS ns, PetscReal *t);

/**
 * @brief Solves the Navier-Stokes equation.
 *
 * @param ns Navier-Stokes solver.
 */
PetscErrorCode FcNSSolve(FcNS ns);

#endif
