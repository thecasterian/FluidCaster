/**
 * @file material.h
 * @brief Materials.
 */

#ifndef MATERIAL_H
#define MATERIAL_H

#include <petscsystypes.h>

/**
 * @brief Material definition.
 */
typedef struct {
    /** Density. */
    PetscReal rho;
    /** Viscosity. */
    PetscReal mu;
} FcMaterial;

#endif
