#include <petscsys.h>
#include "../inc/ns.h"
#include "../inc/private/meshimpl.h"
#include "../inc/private/nsimpl.h"
#include "../inc/private/solutionimpl.h"

PetscErrorCode FcNSCreate(FcMesh mesh, FcSolution sol, FcMaterial mat, FcNSType type, FcNSOps ops, void *data,
                          FcNS *ns) {
    /* Verify the mesh-solution relationship. */
    if (mesh != sol->mesh)
        SETERRQ(mesh->obj.comm, PETSC_ERR_ARG_INCOMP, "Mesh and solution must be associated.");

    /* Allocate memory for the Navier-Stokes solver. */
    PetscCall(PetscNew(ns));
    FcObjectInit((FcObject)(*ns), mesh->obj.comm, "FcNS");
    (*ns)->obj.type = type;
    (*ns)->ops[0] = *ops;
    (*ns)->data = data;

    /* Initialize. */
    (*ns)->mesh = mesh;
    (*ns)->sol = sol;
    (*ns)->mat = mat;
    (*ns)->setup = PETSC_FALSE;
    (*ns)->maxiters = 20;
    (*ns)->t = 0.0;
    (*ns)->dt = 1.0e-3;
    (*ns)->maxsteps = 100;

    return 0;
}

PetscErrorCode FcNSDestroy(FcNS *ns) {
    PetscCall((*ns)->ops->destroy(ns));

    return 0;
}

PetscErrorCode FcNSSetMaxIters(FcNS ns, PetscInt maxiters) {
    ns->maxiters = maxiters;

    return 0;
}

PetscErrorCode FcNSSetTime(FcNS ns, PetscReal t) {
    ns->t = t;

    return 0;
}

PetscErrorCode FcNSSetTimeStep(FcNS ns, PetscReal timestep) {
    ns->dt = timestep;

    return 0;
}

PetscErrorCode FcNSSetMaxSteps(FcNS ns, PetscInt maxsteps) {
    ns->maxsteps = maxsteps;

    return 0;
}

PetscErrorCode FcNSSetFromOptions(FcNS ns) {
    PetscOptionsBegin(ns->obj.comm, "fc_", "Navier-Stokes solver options.", "FcNS");

    /* Common options. */
    PetscCall(PetscOptionsInt("-ns_max_iters", "Maximum number of iterations in a time step", "FcNSSetMaxIters",
                              ns->maxiters, &ns->maxiters, NULL));
    PetscCall(PetscOptionsReal("-ns_time_step", "Time step size", "FcNSSetTimeStep", ns->dt, &ns->dt, NULL));
    PetscCall(PetscOptionsInt("-ns_max_steps", "Maximum number of time steps", "FcNSSetMaxSteps", ns->maxsteps,
                              &ns->maxsteps, NULL));

    /* Type-specific options. */
    PetscCall(ns->ops->setfromoptions(ns, PetscOptionsObject));

    PetscOptionsEnd();

    return 0;
}

PetscErrorCode FcNSSetUp(FcNS ns) {
    PetscCall(ns->ops->setup(ns));

    return 0;
}

PetscErrorCode FcNSGetTime(FcNS ns, PetscReal *t) {
    *t = ns->t;

    return 0;
}

PetscErrorCode FcNSSolve(FcNS ns) {
    PetscCall(ns->ops->solve(ns));

    return 0;
}
