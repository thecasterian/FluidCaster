#include <petscdm.h>
#include "../inc/private/meshimpl.h"
#include "../inc/private/solutionimpl.h"
#include "../inc/solution.h"

PetscErrorCode FcSolutionCreate(FcMesh mesh, FcSolution *sol) {
    /* Allocate memory for the solution. */
    PetscCall(PetscNew(sol));
    FcObjectCreate(*sol, mesh->obj.comm, "FcSolution");

    /* Set the mesh. */
    (*sol)->mesh = mesh;

    /* Create vectors. */
    PetscCall(DMCreateLocalVector(mesh->da, &(*sol)->u));
    PetscCall(DMCreateLocalVector(mesh->da, &(*sol)->v));
    if (mesh->dim == 3)
        PetscCall(DMCreateLocalVector(mesh->da, &(*sol)->w));
    PetscCall(DMCreateLocalVector(mesh->da, &(*sol)->p));

    return 0;
}

PetscErrorCode FcSolutionDestroy(FcSolution *sol) {
    PetscCall(VecDestroy(&(*sol)->u));
    PetscCall(VecDestroy(&(*sol)->v));
    if ((*sol)->mesh->dim == 3)
        PetscCall(VecDestroy(&(*sol)->w));
    PetscCall(VecDestroy(&(*sol)->p));
    PetscCall(PetscFree(*sol));
    *sol = NULL;

    return 0;
}

PetscErrorCode FcSolutionGetVelocityPressureVec(FcSolution sol, Vec *u, Vec *v, Vec *w, Vec *p) {
    if (u)
        *u = sol->u;
    if (v)
        *v = sol->v;
    if (sol->mesh->dim == 3 && w)
        *w = sol->w;
    if (p)
        *p = sol->p;

    return 0;
}
