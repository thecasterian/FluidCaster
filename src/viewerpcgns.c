#include <petscdmda.h>
#include "../inc/private/meshimpl.h"
#include "../inc/private/solutionimpl.h"
#include "../inc/private/viewerpcgnsimpl.h"
#include "../inc/error.h"

#define CgnsCall(...) \
    do { \
        if (__VA_ARGS__) \
            SETERRQ(PETSC_COMM_SELF, FC_ERR_CGNS, "%s", cg_get_error()); \
    } while (0)

static PetscErrorCode FcViewerClose_PCGNS(FcViewer *viewer);
static PetscErrorCode FcViewerViewMesh_PCGNS(FcViewer viewer, FcMesh mesh);
static PetscErrorCode FcViewerViewSolution_PCGNS(FcViewer viewer, FcSolution sol);

PetscErrorCode FcViewerPCGNSOpen(MPI_Comm comm, const char *filename, FcViewer *viewer) {
    struct _FcViewerOps ops;
    VIEWER_PCGNS *pcgns;

    /* Set operations. */
    ops.close = FcViewerClose_PCGNS;
    ops.viewmesh = FcViewerViewMesh_PCGNS;
    ops.viewsol = FcViewerViewSolution_PCGNS;

    /* Create data. */
    PetscCall(PetscNew(&pcgns));

    /* Create viewer. */
    PetscCall(FcViewerCreate(comm, FC_VIEWER_PCGNS, &ops, pcgns, viewer));

    /* Open file. */
    CgnsCall(cgp_open(filename, CG_MODE_WRITE, &pcgns->fn));

    return 0;
}

static PetscErrorCode FcViewerClose_PCGNS(FcViewer *viewer) {
    VIEWER_PCGNS *pcgns = (*viewer)->data;

    /* Close file. */
    CgnsCall(cgp_close(pcgns->fn));

    /* Free data. */
    PetscCall(PetscFree(pcgns));

    /* Free viewer. */
    PetscCall(PetscFree(*viewer));
    *viewer = NULL;

    return 0;
}

static PetscErrorCode FcViewerViewMesh_PCGNS(FcViewer viewer, FcMesh mesh) {
    VIEWER_PCGNS *pcgns = viewer->data;
    FcMeshInfo info;
    cgsize_t zonesize[9] = {0}, rmin[3] = {0}, rmax[3] = {0};
    double *x, *y, *z, hx, hy, hz;
    cgsize_t i, j, k;
    PetscInt cnt;

    /* Get mesh info. */
    PetscCall(FcMeshGetInfo(mesh, &info));

    /* Zone size. */
    if (mesh->dim == 2) {
        zonesize[0] = info.mx + 1;
        zonesize[1] = info.my + 1;
        zonesize[2] = info.mx;
        zonesize[3] = info.my;
    } else {
        zonesize[0] = info.mx + 1;
        zonesize[1] = info.my + 1;
        zonesize[2] = info.mz + 1;
        zonesize[3] = info.mx;
        zonesize[4] = info.my;
        zonesize[5] = info.mz;
    }

    /* Range of vertices. Note that all indices in CGNS is 1-based. */
    if (mesh->dim == 2) {
        rmin[0] = info.xs + 1;
        rmin[1] = info.ys + 1;
        rmax[0] = info.xs + info.xm;
        rmax[1] = info.ys + info.ym;
    } else {
        rmin[0] = info.xs + 1;
        rmin[1] = info.ys + 1;
        rmin[2] = info.zs + 1;
        rmax[0] = info.xs + info.xm;
        rmax[1] = info.ys + info.ym;
        rmax[2] = info.zs + info.zm;
    }
    /* Since the number of vertices are greater than the number of elements, some process must write more. */
    if (rmax[0] == info.mx)
        rmax[0]++;
    if (rmax[1] == info.my)
        rmax[1]++;
    if (mesh->dim == 3 && rmax[2] == info.mz)
        rmax[2]++;

    /* Coordinate of vertices. */
    if (mesh->dim == 2) {
        PetscCall(PetscCalloc1((rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1), &x));
        PetscCall(PetscCalloc1((rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1), &y));
        hx = 1.0 / info.mx;
        hy = 1.0 / info.my;
        cnt = 0;
        for (j = rmin[1]; j <= rmax[1]; j++)
            for (i = rmin[0]; i <= rmax[0]; i++) {
                x[cnt] = (i - 1) * hx;
                y[cnt] = (j - 1) * hy;
                cnt++;
            }
    } else {
        PetscCall(PetscCalloc1((rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1) * (rmax[2] - rmin[2] + 1), &x));
        PetscCall(PetscCalloc1((rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1) * (rmax[2] - rmin[2] + 1), &y));
        PetscCall(PetscCalloc1((rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1) * (rmax[2] - rmin[2] + 1), &z));
        hx = 1.0 / info.mx;
        hy = 1.0 / info.my;
        hz = 1.0 / info.mz;
        cnt = 0;
        for (k = rmin[2]; k <= rmax[2]; k++)
            for (j = rmin[1]; j <= rmax[1]; j++)
                for (i = rmin[0]; i <= rmax[0]; i++) {
                    x[cnt] = (i - 1) * hx;
                    y[cnt] = (j - 1) * hy;
                    z[cnt] = (k - 1) * hz;
                    cnt++;
                }
    }

    /* Create base and zone. */
    CgnsCall(cg_base_write(pcgns->fn, "Base", mesh->dim, mesh->dim, &pcgns->B));
    CgnsCall(cg_zone_write(pcgns->fn, pcgns->B, "Zone", zonesize, Structured, &pcgns->Z));

    /* Write coordinate. */
    CgnsCall(cgp_coord_write(pcgns->fn, pcgns->B, pcgns->Z, RealDouble, "CoordinateX", &pcgns->Cx));
    CgnsCall(cgp_coord_write_data(pcgns->fn, pcgns->B, pcgns->Z, pcgns->Cx, rmin, rmax, x));
    CgnsCall(cgp_coord_write(pcgns->fn, pcgns->B, pcgns->Z, RealDouble, "CoordinateY", &pcgns->Cy));
    CgnsCall(cgp_coord_write_data(pcgns->fn, pcgns->B, pcgns->Z, pcgns->Cy, rmin, rmax, y));
    if (mesh->dim == 3) {
        CgnsCall(cgp_coord_write(pcgns->fn, pcgns->B, pcgns->Z, RealDouble, "CoordinateZ", &pcgns->Cz));
        CgnsCall(cgp_coord_write_data(pcgns->fn, pcgns->B, pcgns->Z, pcgns->Cz, rmin, rmax, z));
    }

    /* Free memory. */
    PetscCall(PetscFree(x));
    PetscCall(PetscFree(y));
    if (mesh->dim == 3)
        PetscCall(PetscFree(z));

    return 0;
}

static PetscErrorCode FcViewerViewSolution_PCGNS(FcViewer viewer, FcSolution sol) {
    VIEWER_PCGNS *pcgns = viewer->data;
    FcMesh mesh = sol->mesh;
    FcMeshInfo info;
    cgsize_t rmin[3] = {0}, rmax[3] = {0};
    double *arruraw, *arrvraw, *arrwraw, *arrpraw;
    PetscInt i, j, k;
    PetscInt cnt;

    /* For CGNS, solution data requires mesh data. */
    if (!viewer->mesh)
        SETERRQ(sol->obj.comm, PETSC_ERR_ARG_WRONGSTATE, "View mesh first for CGNS");

    /* Get mesh info. */
    PetscCall(FcMeshGetInfo(viewer->mesh, &info));

    /* Range of elements. Note that all indices in CGNS is 1-based. */
    if (mesh->dim == 2) {
        rmin[0] = info.xs + 1;
        rmin[1] = info.ys + 1;
        rmax[0] = info.xs + info.xm;
        rmax[1] = info.ys + info.ym;
    } else {
        rmin[0] = info.xs + 1;
        rmin[1] = info.ys + 1;
        rmin[2] = info.zs + 1;
        rmax[0] = info.xs + info.xm;
        rmax[1] = info.ys + info.ym;
        rmax[2] = info.zs + info.zm;
    }

    /* Velocity and pressure. */
    if (mesh->dim == 2) {
        const PetscReal **arru, **arrv, **arrp;
        PetscCall(DMDAVecGetArrayRead(mesh->da, sol->u, &arru));
        PetscCall(DMDAVecGetArrayRead(mesh->da, sol->v, &arrv));
        PetscCall(DMDAVecGetArrayRead(mesh->da, sol->p, &arrp));
        PetscCall(PetscCalloc1((rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1), &arruraw));
        PetscCall(PetscCalloc1((rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1), &arrvraw));
        PetscCall(PetscCalloc1((rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1), &arrpraw));
        cnt = 0;
        for (j = rmin[1]; j <= rmax[1]; j++)
            for (i = rmin[0]; i <= rmax[0]; i++) {
                arruraw[cnt] = arru[j][i];
                arrvraw[cnt] = arrv[j][i];
                arrpraw[cnt] = arrp[j][i];
                cnt++;
            }
        PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->u, &arru));
        PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->v, &arrv));
        PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->p, &arrp));
    } else {
        const PetscReal ***arru, ***arrv, ***arrw, ***arrp;
        PetscCall(DMDAVecGetArrayRead(mesh->da, sol->u, &arru));
        PetscCall(DMDAVecGetArrayRead(mesh->da, sol->v, &arrv));
        PetscCall(DMDAVecGetArrayRead(mesh->da, sol->w, &arrw));
        PetscCall(DMDAVecGetArrayRead(mesh->da, sol->p, &arrp));
        PetscCall(PetscCalloc1((rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1) * (rmax[2] - rmin[2] + 1), &arruraw));
        PetscCall(PetscCalloc1((rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1) * (rmax[2] - rmin[2] + 1), &arrvraw));
        PetscCall(PetscCalloc1((rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1) * (rmax[2] - rmin[2] + 1), &arrwraw));
        PetscCall(PetscCalloc1((rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1) * (rmax[2] - rmin[2] + 1), &arrpraw));
        cnt = 0;
        for (k = rmin[2]; k <= rmax[2]; k++)
            for (j = rmin[1]; j <= rmax[1]; j++)
                for (i = rmin[0]; i <= rmax[0]; i++) {
                    arruraw[cnt] = arru[k][j][i];
                    arrvraw[cnt] = arrv[k][j][i];
                    arrwraw[cnt] = arrw[k][j][i];
                    arrpraw[cnt] = arrp[k][j][i];
                    cnt++;
                }
        PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->u, &arru));
        PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->v, &arrv));
        PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->w, &arrw));
        PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->p, &arrp));
    }

    /* Create section. */
    CgnsCall(cg_sol_write(pcgns->fn, pcgns->B, pcgns->Z, "Solution", CellCenter, &pcgns->S));

    /* Write solution. */
    CgnsCall(cgp_field_write(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, RealDouble, "VelocityX", &pcgns->F));
    CgnsCall(cgp_field_write_data(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, pcgns->F, rmin, rmax, arruraw));
    CgnsCall(cgp_field_write(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, RealDouble, "VelocityY", &pcgns->F));
    CgnsCall(cgp_field_write_data(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, pcgns->F, rmin, rmax, arrvraw));
    if (mesh->dim == 3) {
        CgnsCall(cgp_field_write(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, RealDouble, "VelocityZ", &pcgns->F));
        CgnsCall(cgp_field_write_data(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, pcgns->F, rmin, rmax, arrwraw));
    }
    CgnsCall(cgp_field_write(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, RealDouble, "Pressure", &pcgns->F));
    CgnsCall(cgp_field_write_data(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, pcgns->F, rmin, rmax, arrpraw));

    return 0;
}
