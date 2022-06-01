#include <petscdmda.h>
#include <petscdmstag.h>
#include "../inc/private/meshimpl.h"

static const char *BoundaryTypes[] = {"none", "periodic", "FcMeshBoundaryType", "FC_MESH_BOUNDARY_", NULL};

static PetscErrorCode FcObjectDestroy_Mesh(FcObject *obj);
static PetscErrorCode ConvertBoundaryTypeFromFcMeshToDM(FcMeshBoundaryType fcb, DMBoundaryType *dmb);
static PetscErrorCode ConvertBoundaryTypeFromDMToFcMesh(DMBoundaryType dmb, FcMeshBoundaryType *fcb);

PetscErrorCode FcMeshCreate2d(MPI_Comm comm, FcMeshBoundaryType bx, FcMeshBoundaryType by, PetscInt mx, PetscInt my,
                              PetscInt px, PetscInt py, const PetscInt lx[], const PetscInt ly[], FcMesh *mesh) {
    FcMesh m;

    /* Create the new mesh. */
    PetscCall(PetscNew(&m));

    /* Initialize. */
    FcObjectInit((FcObject)m, comm, "FcMesh");
    m->obj.ops.destroy = FcObjectDestroy_Mesh;
    PetscCall(FcMeshSetDimension(m, 2));
    PetscCall(FcMeshSetBoundaryType(m, bx, by, FC_MESH_BOUNDARY_NONE));
    PetscCall(FcMeshSetSizes(m, mx, my, 1));
    PetscCall(FcMeshSetNumProcs(m, px, py, PETSC_DECIDE));
    PetscCall(FcMeshSetOwnershipRanges(m, lx, ly, NULL));
    m->da = NULL;
    m->stag = NULL;

    /* Get the reference. */
    PetscCall(FcObjectGetReference((FcObject)m, (FcObject *)mesh));

    return 0;
}

PetscErrorCode FcMeshCreate3d(MPI_Comm comm, FcMeshBoundaryType bx, FcMeshBoundaryType by, FcMeshBoundaryType bz,
                              PetscInt mx, PetscInt my, PetscInt mz, PetscInt px, PetscInt py, PetscInt pz,
                              const PetscInt lx[], const PetscInt ly[], const PetscInt lz[], FcMesh *mesh) {
    FcMesh m;

    /* Allocate memory for the mesh. */
    PetscCall(PetscNew(&m));

    /* Initialize. */
    FcObjectInit((FcObject)m, comm, "FcMesh");
    m->obj.ops.destroy = FcObjectDestroy_Mesh;
    PetscCall(FcMeshSetDimension(m, 3));
    PetscCall(FcMeshSetBoundaryType(m, bx, by, bz));
    PetscCall(FcMeshSetSizes(m, mx, my, mz));
    PetscCall(FcMeshSetNumProcs(m, px, py, pz));
    PetscCall(FcMeshSetOwnershipRanges(m, lx, ly, lz));
    m->da = NULL;
    m->stag = NULL;

    /* Get the reference. */
    PetscCall(FcObjectGetReference((FcObject)m, (FcObject *)mesh));

    return 0;
}

PetscErrorCode FcMeshSetDimension(FcMesh mesh, PetscInt dim) {
    mesh->dim = dim;

    return 0;
}

PetscErrorCode FcMeshSetBoundaryType(FcMesh mesh, FcMeshBoundaryType bx, FcMeshBoundaryType by, FcMeshBoundaryType bz) {
    mesh->bx = bx;
    mesh->by = by;
    mesh->bz = bz;

    return 0;
}

PetscErrorCode FcMeshSetSizes(FcMesh mesh, PetscInt mx, PetscInt my, PetscInt mz) {
    mesh->mx = mx;
    mesh->my = my;
    mesh->mz = mz;

    return 0;
}

PetscErrorCode FcMeshRefine(FcMesh mesh, PetscInt n) {
    PetscInt i;

    for (i = 0; i < n; i++) {
        mesh->mx = 2 * mesh->mx - 1;
        mesh->my = 2 * mesh->my - 1;
        mesh->mz = 2 * mesh->mz - 1;
    }

    return 0;
}

PetscErrorCode FcMeshSetNumProcs(FcMesh mesh, PetscInt px, PetscInt py, PetscInt pz) {
    mesh->px = px;
    mesh->py = py;
    mesh->pz = pz;

    return 0;
}

PetscErrorCode FcMeshSetOwnershipRanges(FcMesh mesh, const PetscInt lx[], const PetscInt ly[], const PetscInt lz[]) {
    mesh->lx = lx;
    mesh->ly = ly;
    mesh->lz = lz;

    return 0;
}

PetscErrorCode FcMeshSetFromOptions(FcMesh mesh) {
    PetscInt nrefs = 0;

    PetscOptionsBegin(mesh->obj.comm, "fc_mesh_", "Mesh (FcMesh) options", "FcMesh");

    PetscCall(PetscOptionsInt("-dim", "Dimension", "FcMeshSetDimension", mesh->dim, &mesh->dim, NULL));
    PetscCall(PetscOptionsEnum("-bndry_x", "Boundary type in x direction", "FcMeshSetBoundaryType", BoundaryTypes,
                               (PetscEnum)mesh->bx, (PetscEnum *)&mesh->bx, NULL));
    PetscCall(PetscOptionsEnum("-bndry_y", "Boundary type in y direction", "FcMeshSetBoundaryType", BoundaryTypes,
                               (PetscEnum)mesh->by, (PetscEnum *)&mesh->by, NULL));
    if (mesh->dim == 3)
        PetscCall(PetscOptionsEnum("-bndry_z", "Boundary type in z direction", "FcMeshSetBoundaryType", BoundaryTypes,
                                   (PetscEnum)mesh->bz, (PetscEnum *)&mesh->bz, NULL));
    PetscCall(PetscOptionsInt("-grid_x", "Number of grid points in x direction", "FcMeshSetSizes", mesh->mx, &mesh->mx,
                              NULL));
    PetscCall(PetscOptionsInt("-grid_y", "Number of grid points in y direction", "FcMeshSetSizes", mesh->my, &mesh->my,
                              NULL));
    if (mesh->dim == 3)
        PetscCall(PetscOptionsInt("-grid_z", "Number of grid points in z direction", "FcMeshSetSizes", mesh->mz,
                                  &mesh->mz, NULL));
    PetscCall(PetscOptionsInt("-refine", "Uniformly refine mesh one or more times", "FcMeshRefine", nrefs, &nrefs,
                              NULL));
    PetscCall(PetscOptionsInt("-processors_x", "Number of processors in x direction", "FcMeshSetNumProcs", mesh->px,
                              &mesh->px, NULL));
    PetscCall(PetscOptionsInt("-processors_y", "Number of processors in y direction", "FcMeshSetNumProcs", mesh->py,
                              &mesh->py, NULL));
    if (mesh->dim == 3)
        PetscCall(PetscOptionsInt("-processors_z", "Number of processors in z direction", "FcMeshSetNumProcs", mesh->pz,
                                  &mesh->pz, NULL));

    PetscOptionsEnd();

    FcMeshRefine(mesh, nrefs);

    return 0;
}

PetscErrorCode FcMeshSetUp(FcMesh mesh) {
    DMBoundaryType bx, by, bz;

    PetscCall(ConvertBoundaryTypeFromFcMeshToDM(mesh->bx, &bx));
    PetscCall(ConvertBoundaryTypeFromFcMeshToDM(mesh->by, &by));
    PetscCall(ConvertBoundaryTypeFromFcMeshToDM(mesh->bz, &bz));

    /* Create DMDA. */
    if (mesh->dim == 2)
        PetscCall(DMDACreate2d(mesh->obj.comm, bx, by, DMDA_STENCIL_STAR, mesh->mx, mesh->my, mesh->px, mesh->py, 1, 1,
                               mesh->lx, mesh->ly, &mesh->da));
    else
        PetscCall(DMDACreate3d(mesh->obj.comm, bx, by, bz, DMDA_STENCIL_STAR, mesh->mx, mesh->my, mesh->mz, mesh->px,
                               mesh->py, mesh->pz, 1, 1, mesh->lx, mesh->ly, mesh->lz, &mesh->da));
    PetscCall(DMSetUp(mesh->da));

    /* Clone DMDA. */
    PetscCall(DMClone(mesh->da, &mesh->dau));
    PetscCall(DMClone(mesh->da, &mesh->dav));
    if (mesh->dim == 3)
        PetscCall(DMClone(mesh->da, &mesh->daw));
    PetscCall(DMClone(mesh->da, &mesh->dap));

    /* Create DMStag. */
    if (mesh->dim == 2)
        PetscCall(DMStagCreate2d(mesh->obj.comm, bx, by, mesh->mx, mesh->my, mesh->px, mesh->py, 0, 1, 0,
                                 DMSTAG_STENCIL_STAR, 1, mesh->lx, mesh->ly, &mesh->stag));
    else
        PetscCall(DMStagCreate3d(mesh->obj.comm, bx, by, bz, mesh->mx, mesh->my, mesh->mz, mesh->px, mesh->py, mesh->pz,
                                 0, 0, 1, 0, DMSTAG_STENCIL_STAR, 1, mesh->lx, mesh->ly, mesh->lz, &mesh->stag));
    PetscCall(DMSetUp(mesh->stag));

    return 0;
}

PetscErrorCode FcMeshDestory(FcMesh *mesh) {
    PetscCall(FcObjectRestoreReference((FcObject *)mesh));

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
    PetscCall(ConvertBoundaryTypeFromDMToFcMesh(dmbz, &info->bz));
    info->xs = localinfo.xs;
    info->ys = localinfo.ys;
    info->zs = localinfo.zs;
    info->xm = localinfo.xm;
    info->ym = localinfo.ym;
    info->zm = localinfo.zm;

    return 0;
}

PetscErrorCode FcMeshGetDM(FcMesh mesh, DM *da, DM *stag) {
    *da = mesh->da;
    *stag = mesh->stag;

    return 0;
}

static PetscErrorCode FcObjectDestroy_Mesh(FcObject *obj) {
    FcMesh mesh = (FcMesh)(*obj);

    /* Destroy DMs. */
    PetscCall(DMDestroy(&mesh->da));
    PetscCall(DMDestroy(&mesh->stag));

    /* Free memory. */
    PetscCall(PetscFree(mesh));

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
