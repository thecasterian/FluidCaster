#include <petscsys.h>
#include <petscdmda.h>
#include "../inc/mesh.h"
#include "../inc/solution.h"
#include "../inc/nsfsm.h"
#include "../inc/viewerpcgns.h"

static const char *const help = "FluidCaster\n\n";

int main(int argc, char *argv[]) {
    FcMesh mesh;
    FcSolution sol;
    FcMaterial mat;
    FcNS ns;
    FcViewer viewer;

    /* Initialize PETSc. */
    PetscCall(PetscInitialize(&argc, &argv, NULL, help));

    /* Create mesh. */
    PetscCall(FcMeshCreate2d(PETSC_COMM_WORLD, FC_MESH_BOUNDARY_NONE, FC_MESH_BOUNDARY_NONE, 3, 3, PETSC_DECIDE,
                             PETSC_DECIDE, NULL, NULL, &mesh));
    PetscCall(FcMeshSetFromOptions(mesh));
    PetscCall(FcMeshSetUp(mesh));

    /* Create solution. */
    PetscCall(FcSolutionCreate(mesh, &sol));

    /* Set material. */
    mat.rho = 1.0;
    mat.mu = 1.0e-3;

    /* Create solver. */
    PetscCall(FcNSFSMCreate(mesh, sol, mat, &ns));
    PetscCall(FcNSSetFromOptions(ns));

    /* Solve. */
    PetscCall(FcNSSolve(ns));

    /* View. */
    PetscCall(FcViewerPCGNSOpen(PETSC_COMM_WORLD, "cavity.cgns", &viewer));
    PetscCall(FcMeshView(mesh, viewer));
    PetscCall(FcSolutionView(sol, viewer));
    PetscCall(FcViewerClose(&viewer));

    /* Destroy. */
    PetscCall(FcMeshDestory(&mesh));
    PetscCall(FcSolutionDestroy(&sol));
    PetscCall(FcNSDestroy(&ns));

    /* Finalize PETSc. */
    PetscCall(PetscFinalize());

    return 0;
}
