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

static PetscErrorCode FcObjectDestroy_ViewerPCGNS(FcObject *obj);
static PetscErrorCode FcViewerViewMesh_PCGNS(FcViewer viewer, FcMesh mesh);
static PetscErrorCode FcViewerViewSolution_PCGNS(FcViewer viewer, FcSolution sol);

PetscErrorCode FcViewerPCGNSCreate(MPI_Comm comm, const char *filename, FcViewer *viewer) {
    VIEWER_PCGNS *pcgns;

    /* Create data. */
    PetscCall(PetscNew(&pcgns));

    /* Create viewer. */
    PetscCall(FcViewerCreate(comm, FC_VIEWER_PCGNS, pcgns, viewer));

    /* Set operations. */
    (*viewer)->obj.ops.destroy = FcObjectDestroy_ViewerPCGNS;
    (*viewer)->ops.viewmesh = FcViewerViewMesh_PCGNS;
    (*viewer)->ops.viewsol = FcViewerViewSolution_PCGNS;

    /* Open file. */
    CgnsCall(cgp_open(filename, CG_MODE_WRITE, &pcgns->fn));

    return 0;
}

static PetscErrorCode FcObjectDestroy_ViewerPCGNS(FcObject *obj) {
    FcViewer viewer = (FcViewer)(*obj);
    VIEWER_PCGNS *pcgns = viewer->data;

    /* Close file. */
    CgnsCall(cgp_close(pcgns->fn));

    /* Free data. */
    PetscCall(PetscFree(pcgns));

    /* Restore the references. */
    PetscCall(FcObjectRestoreReference((FcObject *)&viewer->mesh));
    PetscCall(FcObjectRestoreReference((FcObject *)&viewer->sol));

    /* Free viewer. */
    PetscCall(PetscFree(viewer));

    return 0;
}

static PetscErrorCode FcViewerViewMesh_PCGNS(FcViewer viewer, FcMesh mesh) {
    VIEWER_PCGNS *pcgns = viewer->data;
    FcMeshInfo info;
    cgsize_t zonesize[9] = {0}, rmin[3] = {0}, rmax[3] = {0};
    double *x, *y, *z;
    int Cx, Cy, Cz;
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
        cnt = 0;
        for (j = rmin[1]; j <= rmax[1]; j++)
            for (i = rmin[0]; i <= rmax[0]; i++) {
                x[cnt] = mesh->xf[i - 1];
                y[cnt] = mesh->yf[j - 1];
                cnt++;
            }
    } else {
        PetscCall(PetscCalloc1((rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1) * (rmax[2] - rmin[2] + 1), &x));
        PetscCall(PetscCalloc1((rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1) * (rmax[2] - rmin[2] + 1), &y));
        PetscCall(PetscCalloc1((rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1) * (rmax[2] - rmin[2] + 1), &z));
        cnt = 0;
        for (k = rmin[2]; k <= rmax[2]; k++)
            for (j = rmin[1]; j <= rmax[1]; j++)
                for (i = rmin[0]; i <= rmax[0]; i++) {
                    x[cnt] = mesh->xf[i - 1];
                    y[cnt] = mesh->yf[j - 1];
                    z[cnt] = mesh->zf[k - 1];
                    cnt++;
                }
    }

    /* Create base and zone. */
    CgnsCall(cg_base_write(pcgns->fn, "Base", mesh->dim, mesh->dim, &pcgns->B));
    CgnsCall(cg_zone_write(pcgns->fn, pcgns->B, "Zone", zonesize, Structured, &pcgns->Z));

    /* Write coordinate. */
    CgnsCall(cgp_coord_write(pcgns->fn, pcgns->B, pcgns->Z, RealDouble, "CoordinateX", &Cx));
    CgnsCall(cgp_coord_write_data(pcgns->fn, pcgns->B, pcgns->Z, Cx, rmin, rmax, x));
    CgnsCall(cgp_coord_write(pcgns->fn, pcgns->B, pcgns->Z, RealDouble, "CoordinateY", &Cy));
    CgnsCall(cgp_coord_write_data(pcgns->fn, pcgns->B, pcgns->Z, Cy, rmin, rmax, y));
    if (mesh->dim == 3) {
        CgnsCall(cgp_coord_write(pcgns->fn, pcgns->B, pcgns->Z, RealDouble, "CoordinateZ", &Cz));
        CgnsCall(cgp_coord_write_data(pcgns->fn, pcgns->B, pcgns->Z, Cz, rmin, rmax, z));
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
    int Fu, Fv, Fw, Fp;
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
                arruraw[cnt] = arru[j-1][i-1];
                arrvraw[cnt] = arrv[j-1][i-1];
                arrpraw[cnt] = arrp[j-1][i-1];
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
                    arruraw[cnt] = arru[k-1][j-1][i-1];
                    arrvraw[cnt] = arrv[k-1][j-1][i-1];
                    arrwraw[cnt] = arrw[k-1][j-1][i-1];
                    arrpraw[cnt] = arrp[k-1][j-1][i-1];
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
    CgnsCall(cgp_field_write(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, RealDouble, "VelocityX", &Fu));
    CgnsCall(cgp_field_write_data(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, Fu, rmin, rmax, arruraw));
    CgnsCall(cgp_field_write(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, RealDouble, "VelocityY", &Fv));
    CgnsCall(cgp_field_write_data(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, Fv, rmin, rmax, arrvraw));
    if (mesh->dim == 3) {
        CgnsCall(cgp_field_write(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, RealDouble, "VelocityZ", &Fw));
        CgnsCall(cgp_field_write_data(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, Fw, rmin, rmax, arrwraw));
    }
    CgnsCall(cgp_field_write(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, RealDouble, "Pressure", &Fp));
    CgnsCall(cgp_field_write_data(pcgns->fn, pcgns->B, pcgns->Z, pcgns->S, Fp, rmin, rmax, arrpraw));

    return 0;
}
