#include <petscsys.h>
#include "../inc/private/meshimpl.h"
#include "../inc/private/nsimpl.h"
#include "../inc/private/solutionimpl.h"
#include "../inc/private/fsmimpl.h"

static PetscErrorCode FcObjectDestroy_NS(FcObject *obj);

PetscErrorCode FcNSCreate(FcMesh mesh, FcSolution sol, FcMaterial mat, FcNS *ns) {
    FcNS n;

    /* Verify the mesh-solution relationship. */
    if (mesh != sol->mesh)
        SETERRQ(mesh->obj.comm, PETSC_ERR_ARG_INCOMP, "Mesh and solution must be associated.");

    /* Create the new Navier-Stokes solver. */
    PetscCall(PetscNew(&n));

    /* Initialize. */
    FcObjectInit((FcObject)n, mesh->obj.comm, "FcNS");
    n->obj.ops.destroy = FcObjectDestroy_NS;
    PetscCall(FcObjectGetReference((FcObject)mesh, (FcObject *)&n->mesh));
    PetscCall(FcObjectGetReference((FcObject)sol, (FcObject *)&n->sol));
    PetscCall(FcObjectGetReference((FcObject)mat, (FcObject *)&n->mat));
    n->t = 0.0;
    n->dt = 1.0e-3;
    n->maxsteps = 0;

    /* Get the reference. */
    PetscCall(FcObjectGetReference((FcObject)n, (FcObject *)ns));

    return 0;
}

PetscErrorCode FcNSDestroy(FcNS *ns) {
    PetscCall(FcObjectRestoreReference((FcObject *)ns));

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
    PetscOptionsBegin(ns->obj.comm, "fc_ns_", "Navier-Stokes solver options.", "FcNS");

    PetscCall(PetscOptionsReal("-time_step", "Time step size", "FcNSSetTimeStep", ns->dt, &ns->dt, NULL));
    PetscCall(PetscOptionsInt("-max_steps", "Maximum number of time steps", "FcNSSetMaxSteps", ns->maxsteps,
                              &ns->maxsteps, NULL));
    PetscCall(PetscOptionsReal("-time", "Current time", "FcNSSetTime", ns->t, &ns->t, NULL));

    PetscOptionsEnd();

    return 0;
}

PetscErrorCode FcNSSetUp(FcNS ns) {
    PC pc;

    PetscCall(KSPCreate(ns->obj.comm, &ns->kspu));
    PetscCall(KSPSetDM(ns->kspu, ns->mesh->dau));
    if (ns->mesh->dim == 2) {
        PetscCall(KSPSetComputeRHS(ns->kspu, FSMComputeRHSUstar2d, ns));
        PetscCall(KSPSetComputeOperators(ns->kspu, FSMComputeOperatorsUVstar2d, ns));
    } else {
        PetscCall(KSPSetComputeRHS(ns->kspu, FSMComputeRHSUstar3d, ns));
        PetscCall(KSPSetComputeOperators(ns->kspu, FSMComputeOperatorsUVWstar3d, ns));
    }
    PetscCall(KSPGetPC(ns->kspu, &pc));
    PetscCall(PCSetType(pc, PCMG));
    PetscCall(KSPSetFromOptions(ns->kspu));

    PetscCall(KSPCreate(ns->obj.comm, &ns->kspv));
    PetscCall(KSPSetDM(ns->kspv, ns->mesh->dav));
    if (ns->mesh->dim == 2) {
        PetscCall(KSPSetComputeRHS(ns->kspv, FSMComputeRHSVstar2d, ns));
        PetscCall(KSPSetComputeOperators(ns->kspv, FSMComputeOperatorsUVstar2d, ns));
    } else {
        PetscCall(KSPSetComputeRHS(ns->kspv, FSMComputeRHSVstar3d, ns));
        PetscCall(KSPSetComputeOperators(ns->kspv, FSMComputeOperatorsUVWstar3d, ns));
    }
    PetscCall(KSPGetPC(ns->kspv, &pc));
    PetscCall(PCSetType(pc, PCMG));
    PetscCall(KSPSetFromOptions(ns->kspv));

    if (ns->mesh->dim == 3) {
        PetscCall(KSPCreate(ns->obj.comm, &ns->kspw));
        PetscCall(KSPSetDM(ns->kspw, ns->mesh->daw));
        PetscCall(KSPSetComputeRHS(ns->kspw, FSMComputeRHSWstar3d, ns));
        PetscCall(KSPSetComputeOperators(ns->kspw, FSMComputeOperatorsUVWstar3d, ns));
        PetscCall(KSPGetPC(ns->kspw, &pc));
        PetscCall(PCSetType(pc, PCMG));
        PetscCall(KSPSetFromOptions(ns->kspw));
    }

    PetscCall(KSPCreate(ns->obj.comm, &ns->kspp));
    PetscCall(KSPSetDM(ns->kspp, ns->mesh->dap));
    if (ns->mesh->dim == 2) {
        PetscCall(KSPSetComputeRHS(ns->kspp, FSMComputeRHSPprime2d, ns));
        PetscCall(KSPSetComputeOperators(ns->kspp, FSMComputeOperatorsPprime2d, ns));
    } else {
        PetscCall(KSPSetComputeRHS(ns->kspp, FSMComputeRHSPprime3d, ns));
        PetscCall(KSPSetComputeOperators(ns->kspp, FSMComputeOperatorsPprime3d, ns));
    }
    PetscCall(KSPGetPC(ns->kspp, &pc));
    PetscCall(PCSetType(pc, PCMG));
    PetscCall(KSPSetFromOptions(ns->kspp));

    return 0;
}

PetscErrorCode FcNSGetTime(FcNS ns, PetscReal *t) {
    *t = ns->t;

    return 0;
}

PetscErrorCode FcNSSolve(FcNS ns) {
    PetscInt i;

    PetscCall(FcNSSetUp(ns));

    for (i = 0; i < ns->maxsteps; i++) {
        PetscCall(FSMCalculateConvection2d(ns));
        PetscCall(FSMCalculateIntermediateVelocity2d(ns));
        PetscCall(FSMCalculatePressureCorrection2d(ns));
        PetscCall(FSMUpdate2d(ns));
    }

    return 0;
}

static PetscErrorCode FcObjectDestroy_NS(FcObject *obj) {
    FcNS ns = (FcNS)(*obj);

    PetscCall(KSPDestroy(&ns->kspu));
    PetscCall(KSPDestroy(&ns->kspv));
    if (ns->mesh->dim == 3)
        PetscCall(KSPDestroy(&ns->kspw));
    PetscCall(KSPDestroy(&ns->kspp));

    PetscCall(PetscFree(ns));

    return 0;
}
