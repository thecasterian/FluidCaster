/**
 * @file material.h
 * @brief Materials used in the solver.
 */

#ifndef MATERIAL_H
#define MATERIAL_H

#include <petscsystypes.h>

/**
 * @brief Material definition.
 */
typedef struct _p_FcMaterial *FcMaterial;

/**
 * @brief Creates a meterial.
 *
 * @param rho Density.
 * @param mu Viscosity.
 * @param[out] mat Resulting material.
 */
PetscErrorCode FcMaterialCreate(PetscReal rho, PetscReal mu, FcMaterial *mat);

/**
 * @brief Destroys a material.
 *
 * @param[in] mat Material to destroy.
 */
PetscErrorCode FcMaterialDestroy(FcMaterial *mat);

#endif
