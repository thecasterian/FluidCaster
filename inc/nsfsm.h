/**
 * @file nsfsm.h
 * @brief Fractional step method for a Navier-Stokes solver.
 */

#ifndef NS_FSM_H
#define NS_FSM_H

#include "ns.h"
#include "material.h"

/**
 * @brief Creates a Navier-Stokes solver with fractional step method (FSM).
 *
 * @param mesh Mesh where the problem is defined.
 * @param sol Solution of the Navier-Stokes equation.
 * @param mat Fluid definition.
 * @param[out] ns Resulting Navier-Stokes solver.
 */
PetscErrorCode FcNSFSMCreate(FcMesh mesh, FcSolution sol, FcMaterial mat, FcNS *ns);

#endif
