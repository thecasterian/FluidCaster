#include <petscsys.h>
#include "../inc/mesh.h"
#include "../inc/solution.h"

static const char *const help = "FluidCaster\n\n";

int main(int argc, char *argv[]) {
    FcMesh mesh = NULL;
    FcSolution sol = NULL;

    /* Initialize PETSc. */
    PetscInitialize(&argc, &argv, NULL, help);

    /* Create mesh. */
    FcMeshCreate2d(PETSC_COMM_WORLD, FC_MESH_BOUNDARY_NONE, FC_MESH_BOUNDARY_NONE, 3, 3, PETSC_DECIDE, PETSC_DECIDE,
                   NULL, NULL, &mesh);
    FcMeshSetFromOptions(mesh);
    FcMeshSetUp(mesh);

    /* Create solution. */
    FcSolutionCreate(mesh, &sol);

    /* Destroy mesh. */
    FcMeshDestory(&mesh);
    FcSolutionDestroy(&sol);

    /* Finalize PETSc. */
    PetscFinalize();

    return 0;
}
