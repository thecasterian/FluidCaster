#ifndef MATERIAL_IMPL_H
#define MATERIAL_IMPL_H

#include "../material.h"
#include "objectimpl.h"

struct _p_FcMaterial {
    /* Object. */
    struct _p_FcObject obj;

    /* Density. */
    PetscReal rho;
    /* Viscosity. */
    PetscReal mu;
};

#endif
