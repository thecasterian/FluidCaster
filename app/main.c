#include <petscsys.h>
#include <petscdmda.h>
#include "../inc/mesh.h"
#include "../inc/solution.h"
#include "../inc/nsfsm.h"

static const char *const help = "FluidCaster\n\n";

int main(int argc, char *argv[]) {
    FcMesh mesh = NULL;
    FcSolution sol = NULL;
    FcMaterial mat;
    FcNS ns = NULL;

    /* Initialize PETSc. */
    PetscInitialize(&argc, &argv, NULL, help);

    /* Create mesh. */
    FcMeshCreate2d(PETSC_COMM_WORLD, FC_MESH_BOUNDARY_NONE, FC_MESH_BOUNDARY_NONE, 3, 3, PETSC_DECIDE, PETSC_DECIDE,
                   NULL, NULL, &mesh);
    FcMeshSetFromOptions(mesh);
    FcMeshSetUp(mesh);

    /* Create solution. */
    FcSolutionCreate(mesh, &sol);

    /* Set material. */
    mat.rho = 1.0;
    mat.mu = 1.0e-3;

    /* Create solver. */
    FcNSFSMCreate(mesh, sol, mat, &ns);
    FcNSSetFromOptions(ns);

    /* Solve. */
    FcNSSolve(ns);

    {
        DM da, stag;
        FcMeshInfo info;
        Vec u;
        FILE *fp;

        FcMeshGetDM(mesh, &da, &stag);
        FcMeshGetInfo(mesh, &info);
        FcSolutionGetVelocityPressureVec(sol, &u, NULL, NULL, NULL);

        fp = fopen("u.txt", "w");
        if (fp) {
            const PetscReal **arru;
            PetscInt i, j;

            PetscCall(DMDAVecGetArrayRead(da, u, &arru));
            for (i = 0; i < info.mx; i++) {
                for (j = 0; j < info.my; j++)
                    PetscFPrintf(PETSC_COMM_WORLD, fp, "%f ", (double)arru[j][i]);
                PetscFPrintf(PETSC_COMM_WORLD, fp, "\n");
            }
            PetscCall(DMDAVecRestoreArrayRead(da, u, &arru));

            fclose(fp);
        }
    }

    /* Destroy. */
    FcMeshDestory(&mesh);
    FcSolutionDestroy(&sol);
    FcNSDestroy(&ns);

    /* Finalize PETSc. */
    PetscFinalize();

    return 0;
}
