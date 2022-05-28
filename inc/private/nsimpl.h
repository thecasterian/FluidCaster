#ifndef NS_IMPL_H
#define NS_IMPL_H

#include <petscsystypes.h>
#include "../material.h"
#include "../ns.h"
#include "objectimpl.h"

typedef struct _FcNSOps *FcNSOps;
struct _FcNSOps {
    PetscErrorCode (*destroy)(FcNS *);
    PetscErrorCode (*setfromoptions)(FcNS, PetscOptionItems *);
    PetscErrorCode (*setup)(FcNS);
    PetscErrorCode (*solve)(FcNS);
};

struct _p_FcNS {
    /** Object. */
    struct _p_FcObject obj;
    /** Operations. */
    struct _FcNSOps ops[1];
    /** Type-specific data. */
    void *data;

    /** Mesh. */
    FcMesh mesh;
    /** Solution. */
    FcSolution sol;
    /** Material. */
    FcMaterial mat;
    /** Is set up? */
    PetscBool setup;

    /** Maximum nuber of iterations in a time step. */
    PetscInt maxiters;

    /** Current time. */
    PetscReal t;
    /** Type step size. */
    PetscReal dt;
    /** Maximum number of time steps. */
    PetscInt maxsteps;
};

PetscErrorCode FcNSCreate(FcMesh mesh, FcSolution sol, FcMaterial mat, FcNSType type, FcNSOps ops, void *data, FcNS *ns);

#endif
