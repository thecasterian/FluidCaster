#ifndef NS_IMPL_H
#define NS_IMPL_H

#include <petscksp.h>
#include <petscsystypes.h>
#include "../material.h"
#include "../ns.h"
#include "objectimpl.h"

struct _p_FcNS {
    /* Object. */
    struct _p_FcObject obj;

    /* Mesh. */
    FcMesh mesh;
    /* Solution. */
    FcSolution sol;
    /* Material. */
    FcMaterial mat;

    /* Current time. */
    PetscReal t;
    /* Type step size. */
    PetscReal dt;
    /* Maximum number of time steps. */
    PetscInt maxsteps;

    /* KSP for solving the intermediate x-velocity. */
    KSP kspu;
    /* KSP for solving the intermediate y-velocity. */
    KSP kspv;
    /* KSP for solving the intermediate z-velocity. */
    KSP kspw;
    /* KSP for solving the pressure correction. */
    KSP kspp;
};

#endif
