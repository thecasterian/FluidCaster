#include <petscdm.h>
#include "../inc/private/meshimpl.h"
#include "../inc/private/solutionimpl.h"

static PetscErrorCode FcObjectDestroy_Solution(FcObject *obj);

PetscErrorCode FcSolutionCreate(FcMesh mesh, FcSolution *sol) {
    FcSolution s;

    /* Create the new solution. */
    PetscCall(PetscNew(&s));

    /* Initialize. */
    FcObjectInit((FcObject)s, mesh->obj.comm, "FcSolution");
    s->obj.ops.destroy = FcObjectDestroy_Solution;

    /* Set the mesh. */
    FcObjectGetReference((FcObject)mesh, (FcObject *)&s->mesh);

    /* Create vectors. */
    PetscCall(DMCreateLocalVector(mesh->da, &s->u));
    PetscCall(VecDuplicate(s->u, &s->v));
    if (mesh->dim == 3)
        PetscCall(VecDuplicate(s->u, &s->w));
    PetscCall(VecDuplicate(s->u, &s->p_true));

    PetscCall(VecDuplicate(s->u, &s->p));
    PetscCall(DMCreateLocalVector(mesh->stag, &s->UVW));
    PetscCall(VecDuplicate(s->u, &s->u_star));
    PetscCall(VecDuplicate(s->u, &s->v_star));
    if (mesh->dim == 3)
        PetscCall(VecDuplicate(s->u, &s->w_star));
    PetscCall(VecDuplicate(s->UVW, &s->UVW_star));
    PetscCall(VecDuplicate(s->u, &s->p_prime));
    PetscCall(VecDuplicate(s->u, &s->Nu));
    PetscCall(VecDuplicate(s->u, &s->Nv));
    if (mesh->dim == 3)
        PetscCall(VecDuplicate(s->u, &s->Nw));
    PetscCall(VecDuplicate(s->u, &s->p_prev));
    PetscCall(VecDuplicate(s->u, &s->Nu_prev));
    PetscCall(VecDuplicate(s->u, &s->Nv_prev));
    if (mesh->dim == 3)
        PetscCall(VecDuplicate(s->u, &s->Nw_prev));
    PetscCall(VecDuplicate(s->u, &s->u_tilde));
    PetscCall(VecDuplicate(s->u, &s->v_tilde));
    if (mesh->dim == 3)
        PetscCall(VecDuplicate(s->u, &s->w_tilde));

    /* Get the reference. */
    PetscCall(FcObjectGetReference((FcObject)s, (FcObject *)sol));

    return 0;
}

PetscErrorCode FcSolutionDestroy(FcSolution *sol) {
    PetscCall(FcObjectRestoreReference((FcObject *)sol));

    return 0;
}

PetscErrorCode FcSolutionGetVelocityVec(FcSolution sol, Vec *u, Vec *v, Vec *w) {
    if (u)
        *u = sol->u;
    if (v)
        *v = sol->v;
    if (sol->mesh->dim == 3 && w)
        *w = sol->w;

    return 0;
}

PetscErrorCode FcSolutionGetTruePressureVec(FcSolution sol, Vec *p_true) {
    if (p_true)
        *p_true = sol->p_true;

    return 0;
}

static PetscErrorCode FcObjectDestroy_Solution(FcObject *obj) {
    FcSolution sol = (FcSolution)(*obj);

    /* Destroy vectors. */
    PetscCall(VecDestroy(&sol->u));
    PetscCall(VecDestroy(&sol->v));
    if (sol->mesh->dim == 3)
        PetscCall(VecDestroy(&sol->w));
    PetscCall(VecDestroy(&sol->p_true));

    PetscCall(VecDestroy(&sol->p));
    PetscCall(VecDestroy(&sol->UVW));
    PetscCall(VecDestroy(&sol->u_star));
    PetscCall(VecDestroy(&sol->v_star));
    if (sol->mesh->dim == 3)
        PetscCall(VecDestroy(&sol->w_star));
    PetscCall(VecDestroy(&sol->UVW_star));
    PetscCall(VecDestroy(&sol->p_prime));
    PetscCall(VecDestroy(&sol->Nu));
    PetscCall(VecDestroy(&sol->Nv));
    if (sol->mesh->dim == 3)
        PetscCall(VecDestroy(&sol->Nw));
    PetscCall(VecDestroy(&sol->p_prev));
    PetscCall(VecDestroy(&sol->Nu_prev));
    PetscCall(VecDestroy(&sol->Nv_prev));
    if (sol->mesh->dim == 3)
        PetscCall(VecDestroy(&sol->Nw_prev));
    PetscCall(VecDestroy(&sol->u_tilde));
    PetscCall(VecDestroy(&sol->v_tilde));
    if (sol->mesh->dim == 3)
        PetscCall(VecDestroy(&sol->w_tilde));

    /* Decrease the reference counter. */
    FcObjectRestoreReference((FcObject *)&sol->mesh);

    /* Free memory. */
    PetscCall(PetscFree(sol));

    return 0;
}
