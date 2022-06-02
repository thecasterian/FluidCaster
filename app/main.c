#include <petscsys.h>
#include <petscdmda.h>
#include "../inc/material.h"
#include "../inc/mesh.h"
#include "../inc/ns.h"
#include "../inc/solution.h"
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

    /* Create material. */
    PetscCall(FcMaterialCreate(1.0, 1.0e-3, &mat));

    /* Create NS solver. */
    PetscCall(FcNSCreate(mesh, sol, mat, &ns));
    PetscCall(FcNSSetFromOptions(ns));
    PetscCall(FcNSSetUp(ns));

    /* Solve. */
    PetscCall(FcNSSolve(ns));

    /* Save as CGNS. */
    PetscCall(FcViewerPCGNSCreate(PETSC_COMM_WORLD, "cavity.cgns", &viewer));
    PetscCall(FcMeshView(mesh, viewer));
    PetscCall(FcSolutionView(sol, viewer));
    PetscCall(FcViewerDestroy(&viewer));

    /* Destroy. */
    PetscCall(FcMeshDestory(&mesh));
    PetscCall(FcSolutionDestroy(&sol));
    PetscCall(FcMaterialDestroy(&mat));
    PetscCall(FcNSDestroy(&ns));

    /* Finalize PETSc. */
    PetscCall(PetscFinalize());

    return 0;
}
