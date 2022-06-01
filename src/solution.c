#include <petscdm.h>
#include "../inc/private/meshimpl.h"
#include "../inc/private/solutionimpl.h"

PetscErrorCode FcSolutionCreate(FcMesh mesh, FcSolution *sol) {
    /* Allocate memory for the solution. */
    PetscCall(PetscNew(sol));

    /* Set the mesh. */
    (*sol)->mesh = mesh;

    /* Create vectors. */
    PetscCall(DMCreateLocalVector(mesh->da, &(*sol)->u));
    PetscCall(VecDuplicate((*sol)->u, &(*sol)->v));
    if (mesh->dim == 3)
        PetscCall(VecDuplicate((*sol)->u, &(*sol)->w));
    PetscCall(VecDuplicate((*sol)->u, &(*sol)->p_true));

    PetscCall(VecDuplicate((*sol)->u, &(*sol)->p));
    PetscCall(DMCreateLocalVector(mesh->stag, &(*sol)->UVW));
    PetscCall(VecDuplicate((*sol)->u, &(*sol)->u_star));
    PetscCall(VecDuplicate((*sol)->u, &(*sol)->v_star));
    if (mesh->dim == 3)
        PetscCall(VecDuplicate((*sol)->u, &(*sol)->w_star));
    PetscCall(VecDuplicate((*sol)->UVW, &(*sol)->UVW_star));
    PetscCall(VecDuplicate((*sol)->u, &(*sol)->p_prime));
    PetscCall(VecDuplicate((*sol)->u, &(*sol)->Nu));
    PetscCall(VecDuplicate((*sol)->u, &(*sol)->Nv));
    if (mesh->dim == 3)
        PetscCall(VecDuplicate((*sol)->u, &(*sol)->Nw));
    PetscCall(VecDuplicate((*sol)->u, &(*sol)->p_prev));
    PetscCall(VecDuplicate((*sol)->u, &(*sol)->Nu_prev));
    PetscCall(VecDuplicate((*sol)->u, &(*sol)->Nv_prev));
    if (mesh->dim == 3)
        PetscCall(VecDuplicate((*sol)->u, &(*sol)->Nw_prev));
    PetscCall(VecDuplicate((*sol)->u, &(*sol)->u_tilde));
    PetscCall(VecDuplicate((*sol)->u, &(*sol)->v_tilde));
    if (mesh->dim == 3)
        PetscCall(VecDuplicate((*sol)->u, &(*sol)->w_tilde));

    return 0;
}

PetscErrorCode FcSolutionDestroy(FcSolution *sol) {
    /* Destroy vectors. */
    PetscCall(VecDestroy(&(*sol)->u));
    PetscCall(VecDestroy(&(*sol)->v));
    if ((*sol)->mesh->dim == 3)
        PetscCall(VecDestroy(&(*sol)->w));
    PetscCall(VecDestroy(&(*sol)->p_true));

    PetscCall(VecDestroy(&(*sol)->p));
    PetscCall(VecDestroy(&(*sol)->UVW));
    PetscCall(VecDestroy(&(*sol)->u_star));
    PetscCall(VecDestroy(&(*sol)->v_star));
    if ((*sol)->mesh->dim == 3)
        PetscCall(VecDestroy(&(*sol)->w_star));
    PetscCall(VecDestroy(&(*sol)->UVW_star));
    PetscCall(VecDestroy(&(*sol)->p_prime));
    PetscCall(VecDestroy(&(*sol)->Nu));
    PetscCall(VecDestroy(&(*sol)->Nv));
    if ((*sol)->mesh->dim == 3)
        PetscCall(VecDestroy(&(*sol)->Nw));
    PetscCall(VecDestroy(&(*sol)->p_prev));
    PetscCall(VecDestroy(&(*sol)->Nu_prev));
    PetscCall(VecDestroy(&(*sol)->Nv_prev));
    if ((*sol)->mesh->dim == 3)
        PetscCall(VecDestroy(&(*sol)->Nw_prev));
    PetscCall(VecDestroy(&(*sol)->u_tilde));
    PetscCall(VecDestroy(&(*sol)->v_tilde));
    if ((*sol)->mesh->dim == 3)
        PetscCall(VecDestroy(&(*sol)->w_tilde));

    /* Free memory. */
    PetscCall(PetscFree(*sol));
    *sol = NULL;

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
