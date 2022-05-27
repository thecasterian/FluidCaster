#include <petscdmda.h>
#include <petscdmstag.h>
#include "../inc/mesh.h"
#include "../inc/private/meshimpl.h"

static PetscErrorCode ConvertBoundaryTypeFromFcMeshToDM(FcMeshBoundaryType fcb, DMBoundaryType *dmb);
static PetscErrorCode ConvertBoundaryTypeFromDMToFcMesh(DMBoundaryType dmb, FcMeshBoundaryType *fcb);

PetscErrorCode FcMeshCreate2d(MPI_Comm comm, FcMeshBoundaryType bx, FcMeshBoundaryType by, PetscInt mx, PetscInt my,
                              PetscInt px, PetscInt py, const PetscInt lx[], const PetscInt ly[], FcMesh *mesh) {
    DMBoundaryType dmbx, dmby;

    /* Allocate memory for the mesh. */
    PetscCall(PetscNew(mesh));
    (*mesh)->dim = 2;
    (*mesh)->da = NULL;
    (*mesh)->stag = NULL;

    /* Convert FcMeshBoundaryType to DMBoundaryType. */
    PetscCall(ConvertBoundaryTypeFromFcMeshToDM(bx, &dmbx));
    PetscCall(ConvertBoundaryTypeFromFcMeshToDM(by, &dmby));

    /* Create DMDA. */
    PetscCall(DMDACreate2d(comm, dmbx, dmby, DMDA_STENCIL_STAR, mx, my, px, py, 1, 1, lx, ly, &(*mesh)->da));

    return 0;
}

// TODO: Implement FcMeshCreate3d().

PetscErrorCode FcMeshSetFromOptions(FcMesh mesh) {
    PetscCall(DMSetFromOptions(mesh->da));

    return 0;
}

PetscErrorCode FcMeshSetUp(FcMesh mesh) {
    MPI_Comm comm;
    PetscInt dim, mx, my, mz, px, py, pz;
    const PetscInt *lx, *ly, *lz;
    DMBoundaryType bx, by, bz;

    /* Set up DMDA. */
    PetscCall(DMSetUp(mesh->da));

    /* Get DMDA informations. */
    PetscCall(PetscObjectGetComm((PetscObject)mesh->da, &comm));
    PetscCall(DMDAGetInfo(mesh->da, &dim, &mx, &my, &mz, &px, &py, &pz, NULL, NULL, &bx, &by, &bz, NULL));
    PetscCall(DMDAGetOwnershipRanges(mesh->da, &lx, &ly, &lz));

    /* Create DMStag. */
    if (dim == 2)
        DMStagCreate2d(comm, bx, by, mx, my, px, py, 0, 1, 0, DMSTAG_STENCIL_STAR, 1, lx, ly, &mesh->stag);
    // TODO: Handle 3d.

    return 0;
}

PetscErrorCode FcMeshDestory(FcMesh *mesh) {
    PetscCall(DMDestroy(&(*mesh)->da));
    PetscCall(DMDestroy(&(*mesh)->stag));
    PetscCall(PetscFree(*mesh));
    *mesh = NULL;

    return 0;
}

PetscErrorCode FcMeshGetInfo(FcMesh mesh, FcMeshInfo *info) {
    DMBoundaryType dmbx, dmby, dmbz;
    DMDALocalInfo localinfo;

    PetscCall(DMDAGetInfo(mesh->da, &info->dim, &info->mx, &info->my, &info->mz, &info->px, &info->py, &info->pz, NULL, NULL,
                          &dmbx, &dmby, &dmbz, NULL));
    PetscCall(DMDAGetLocalInfo(mesh->da, &localinfo));

    PetscCall(ConvertBoundaryTypeFromDMToFcMesh(dmbx, &info->bx));
    PetscCall(ConvertBoundaryTypeFromDMToFcMesh(dmby, &info->by));
    if (info->dim == 3)
        PetscCall(ConvertBoundaryTypeFromDMToFcMesh(dmbz, &info->bz));
    info->xs = localinfo.xs;
    info->ys = localinfo.ys;
    info->zs = localinfo.zs;

    return 0;
}

PetscErrorCode FcMeshGetDM(FcMesh mesh, DM *da, DM *stag) {
    *da = mesh->da;
    *stag = mesh->stag;

    return 0;
}

static PetscErrorCode ConvertBoundaryTypeFromFcMeshToDM(FcMeshBoundaryType fcb, DMBoundaryType *dmb) {
    switch (fcb) {
        case FC_MESH_BOUNDARY_NONE:
            *dmb = DM_BOUNDARY_GHOSTED;
            break;
        case FC_MESH_BOUNDARY_PERIODIC:
            *dmb = DM_BOUNDARY_PERIODIC;
            break;
        default:
            return PETSC_ERR_SUP;
    }

    return 0;
}

static PetscErrorCode ConvertBoundaryTypeFromDMToFcMesh(DMBoundaryType dmb, FcMeshBoundaryType *fcb) {
    switch (dmb) {
        case DM_BOUNDARY_GHOSTED:
            *fcb = FC_MESH_BOUNDARY_NONE;
            break;
        case DM_BOUNDARY_PERIODIC:
            *fcb = FC_MESH_BOUNDARY_PERIODIC;
            break;
        default:
            return PETSC_ERR_SUP;
    }

    return 0;
}
