#include <petscdmda.h>
#include <petscdmstag.h>
#include <petscoptions.h>
#include <petscsys.h>
#include "../inc/private/materialimpl.h"
#include "../inc/private/meshimpl.h"
#include "../inc/private/nsimpl.h"
#include "../inc/private/solutionimpl.h"

PetscErrorCode CalculateConvection(FcNS ns) {
    FcMesh mesh = ns->mesh;
    FcSolution sol = ns->sol;

    DMDALocalInfo info;
    PetscReal hx, hy;
    PetscReal **arrNu, **arrNv;
    const PetscReal **arru, **arrv, ***arrUV;
    PetscInt iU, iV;
    PetscReal up, uw, ue, us, un, vp, vw, ve, vs, vn;
    PetscInt i, j;

    PetscCall(DMDAGetLocalInfo(mesh->da, &info));
    hx = 1.0 / info.mx;
    hy = 1.0 / info.my;

    PetscCall(DMDAVecGetArray(mesh->da, sol->Nu, &arrNu));
    PetscCall(DMDAVecGetArray(mesh->da, sol->Nv, &arrNv));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->u, &arru));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->v, &arrv));
    PetscCall(DMStagVecGetArrayRead(mesh->stag, sol->UVW, &arrUV));

    PetscCall(DMStagGetLocationSlot(mesh->stag, DMSTAG_LEFT, 0, &iU));
    PetscCall(DMStagGetLocationSlot(mesh->stag, DMSTAG_DOWN, 0, &iV));

    for (j = info.ys; j < info.ys + info.ym; j++)
        for (i = info.xs; i < info.xs + info.xm; i++) {
            up = arru[j][i];
            vp = arrv[j][i];
            /* Left wall. */
            if (i == 0) {
                uw = -arru[j][i];
                vw = -arrv[j][i];
            } else {
                uw = arru[j][i-1];
                vw = arrv[j][i-1];
            }
            /* Right wall. */
            if (i == info.mx - 1) {
                ue = -arru[j][i];
                ve = -arrv[j][i];
            } else {
                ue = arru[j][i+1];
                ve = arrv[j][i+1];
            }
            /* Bottom wall. */
            if (j == 0) {
                us = -arru[j][i];
                vs = -arrv[j][i];
            } else {
                us = arru[j-1][i];
                vs = arrv[j-1][i];
            }
            /* Top wall. */
            if (j == info.my - 1) {
                un = 2.0 - arru[j][i];
                vn = -arrv[j][i];
            } else {
                un = arru[j+1][i];
                vn = arrv[j+1][i];
            }

            arrNu[j][i] = (arrUV[j][i+1][iU] * (up + ue) / 2.0 - arrUV[j][i][iU] * (uw + up) / 2.0) / hx
                          + (arrUV[j+1][i][iV] * (up + un) / 2.0 - arrUV[j][i][iV] * (us + up) / 2.0) / hy;
            arrNv[j][i] = (arrUV[j][i+1][iU] * (vp + ve) / 2.0 - arrUV[j][i][iU] * (vw + vp) / 2.0) / hx
                          + (arrUV[j+1][i][iV] * (vp + vn) / 2.0 - arrUV[j][i][iV] * (vs + vp) / 2.0) / hy;
        }

    PetscCall(DMDAVecRestoreArray(mesh->da, sol->Nu, &arrNu));
    PetscCall(DMDAVecRestoreArray(mesh->da, sol->Nv, &arrNv));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->u, &arru));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->v, &arrv));
    PetscCall(DMStagVecRestoreArrayRead(mesh->stag, sol->UVW, &arrUV));

    return 0;
}

PetscErrorCode CalculateIntermediateVelocity(FcNS ns) {
    FcMesh mesh = ns->mesh;
    FcSolution sol = ns->sol;
    FcMaterial mat = ns->mat;

    Vec x;
    DMDALocalInfo info;
    PetscReal hx, hy;
    PetscReal **arru_tilde, **arrv_tilde, ***arrUV_star;
    const PetscReal **arru_star, **arrv_star, **arrp, **arru_tilde_rd, **arrv_tilde_rd;
    PetscInt iU, iV;
    PetscReal pw, pe, ps, pn;
    PetscInt i, j;

    /* Calculate cell-centered intermediate velocity. */
    PetscCall(KSPSolve(ns->kspu, NULL, NULL));
    PetscCall(KSPGetSolution(ns->kspu, &x));
    PetscCall(DMGlobalToLocal(mesh->dau, x, INSERT_VALUES, sol->u_star));

    PetscCall(KSPSolve(ns->kspv, NULL, NULL));
    PetscCall(KSPGetSolution(ns->kspv, &x));
    PetscCall(DMGlobalToLocal(mesh->dav, x, INSERT_VALUES, sol->v_star));

    /* Calculate face-centered intermediate velocity. */
    PetscCall(DMDAGetLocalInfo(mesh->da, &info));
    hx = 1.0 / info.mx;
    hy = 1.0 / info.my;

    PetscCall(DMDAVecGetArray(mesh->da, sol->u_tilde, &arru_tilde));
    PetscCall(DMDAVecGetArray(mesh->da, sol->v_tilde, &arrv_tilde));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->u_star, &arru_star));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->v_star, &arrv_star));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->p, &arrp));

    PetscCall(DMStagGetLocationSlot(mesh->stag, DMSTAG_LEFT, 0, &iU));
    PetscCall(DMStagGetLocationSlot(mesh->stag, DMSTAG_DOWN, 0, &iV));

    for (j = info.ys; j < info.ys + info.ym; j++)
        for (i = info.xs; i < info.xs + info.xm; i++) {
            /* Left wall. */
            if (i == 0)
                pw = arrp[j][i];
            else
                pw = arrp[j][i-1];
            /* Right wall. */
            if (i == info.mx - 1)
                pe = arrp[j][i];
            else
                pe = arrp[j][i+1];
            /* Bottom wall. */
            if (j == 0)
                ps = arrp[j][i];
            else
                ps = arrp[j-1][i];
            /* Top wall. */
            if (j == info.my - 1)
                pn = arrp[j][i];
            else
                pn = arrp[j+1][i];

            arru_tilde[j][i] = arru_star[j][i] + ns->dt/mat->rho * (pe - pw) / (2.0 * hx);
            arrv_tilde[j][i] = arrv_star[j][i] + ns->dt/mat->rho * (pn - ps) / (2.0 * hy);
        }

    PetscCall(DMDAVecRestoreArray(mesh->da, sol->u_tilde, &arru_tilde));
    PetscCall(DMDAVecRestoreArray(mesh->da, sol->v_tilde, &arrv_tilde));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->u_star, &arru_star));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->v_star, &arrv_star));

    PetscCall(DMLocalToLocalBegin(mesh->da, sol->u_tilde, INSERT_VALUES, sol->u_tilde));
    PetscCall(DMLocalToLocalEnd(mesh->da, sol->u_tilde, INSERT_VALUES, sol->u_tilde));
    PetscCall(DMLocalToLocalBegin(mesh->da, sol->v_tilde, INSERT_VALUES, sol->v_tilde));
    PetscCall(DMLocalToLocalEnd(mesh->da, sol->v_tilde, INSERT_VALUES, sol->v_tilde));

    PetscCall(DMStagVecGetArray(mesh->stag, sol->UVW_star, &arrUV_star));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->u_tilde, &arru_tilde_rd));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->v_tilde, &arrv_tilde_rd));

    for (j = info.ys; j < info.ys + info.ym + 1; j++)
        for (i = info.xs; i < info.xs + info.xm + 1; i++) {
            if (j < info.my) {
                /* Left wall. */
                if (i == 0)
                    arrUV_star[j][i][iU] = 0.0;
                /* Right wall. */
                else if (i == info.mx)
                    arrUV_star[j][i][iU] = 0.0;
                else
                    arrUV_star[j][i][iU] = (arru_tilde_rd[j][i-1] + arru_tilde_rd[j][i]) / 2.0
                                           - ns->dt/mat->rho * (arrp[j][i] - arrp[j][i-1]) / hx;
            }
            if (i < info.mx) {
                /* Bottom wall. */
                if (j == 0)
                    arrUV_star[j][i][iV] = 0.0;
                /* Top wall. */
                else if (j == info.my)
                    arrUV_star[j][i][iV] = 0.0;
                else
                    arrUV_star[j][i][iV] = (arrv_tilde_rd[j-1][i] + arrv_tilde_rd[j][i]) / 2.0
                                           - ns->dt/mat->rho * (arrp[j][i] - arrp[j-1][i]) / hy;
            }
        }

    PetscCall(DMStagVecRestoreArray(mesh->stag, sol->UVW_star, &arrUV_star));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->u_tilde, &arru_tilde_rd));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->v_tilde, &arrv_tilde_rd));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->p, &arrp));

    return 0;
}

PetscErrorCode CalculatePressureCorrection(FcNS ns) {
    Vec x;

    PetscCall(KSPSolve(ns->kspp, NULL, NULL));
    PetscCall(KSPGetSolution(ns->kspp, &x));
    PetscCall(DMGlobalToLocal(ns->mesh->dap, x, INSERT_VALUES, ns->sol->p_prime));

    return 0;
}

PetscErrorCode Update(FcNS ns) {
    FcMesh mesh = ns->mesh;
    FcSolution sol = ns->sol;
    FcMaterial mat = ns->mat;

    DMDALocalInfo info;
    PetscReal hx, hy;
    PetscReal **arru, **arrv, **arrp, ***arrUV;
    const PetscReal **arru_star, **arrv_star, **arrp_prime, ***arrUV_star;
    PetscInt iU, iV;
    PetscReal ppp, ppw, ppe, pps, ppn;
    PetscInt i, j;

    PetscCall(DMDAGetLocalInfo(mesh->da, &info));
    hx = 1.0 / info.mx;
    hy = 1.0 / info.my;

    PetscCall(DMDAVecGetArray(mesh->da, sol->u, &arru));
    PetscCall(DMDAVecGetArray(mesh->da, sol->v, &arrv));
    PetscCall(DMDAVecGetArray(mesh->da, sol->p, &arrp));
    PetscCall(DMStagVecGetArray(mesh->stag, sol->UVW, &arrUV));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->u_star, &arru_star));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->v_star, &arrv_star));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->p_prime, &arrp_prime));
    PetscCall(DMStagVecGetArrayRead(mesh->stag, sol->UVW_star, &arrUV_star));

    PetscCall(DMStagGetLocationSlot(mesh->stag, DMSTAG_LEFT, 0, &iU));
    PetscCall(DMStagGetLocationSlot(mesh->stag, DMSTAG_DOWN, 0, &iV));

    for (j = info.ys; j < info.ys + info.ym; j++)
        for (i = info.xs; i < info.xs + info.xm; i++) {
            ppp = arrp_prime[j][i];
            /* Left wall. */
            if (i == 0)
                ppw = arrp_prime[j][i];
            else
                ppw = arrp_prime[j][i-1];
            /* Right wall. */
            if (i == info.mx - 1)
                ppe = arrp_prime[j][i];
            else
                ppe = arrp_prime[j][i+1];
            /* Bottom wall. */
            if (j == 0)
                pps = arrp_prime[j][i];
            else
                pps = arrp_prime[j-1][i];
            /* Top wall. */
            if (j == info.my - 1)
                ppn = arrp_prime[j][i];
            else
                ppn = arrp_prime[j+1][i];

            arru[j][i] = arru_star[j][i] - ns->dt/mat->rho * (ppe - ppw) / (2.0*hx);
            arrv[j][i] = arrv_star[j][i] - ns->dt/mat->rho * (ppn - pps) / (2.0*hy);
            arrp[j][i] += arrp_prime[j][i]
                          - mat->mu*ns->dt/(2.0*mat->rho) * ((ppe - 2.0*ppp + ppw) / (hx*hx)
                                                             + ((ppn - 2.0*ppp + pps) / (hy*hy)));
        }

    for (j = info.ys; j < info.ys + info.ym + 1; j++)
        for (i = info.xs; i < info.xs + info.xm + 1; i++) {
            if (j < info.my) {
                /* Left wall. */
                if (i == 0)
                    arrUV[j][i][iU] = 0.0;
                /* Right wall. */
                else if (i == info.mx)
                    arrUV[j][i][iU] = 0.0;
                else
                    arrUV[j][i][iU] = arrUV_star[j][i][iU]
                                      - ns->dt/mat->rho * (arrp_prime[j][i] - arrp_prime[j][i-1]) / hx;
            }
            if (i < info.mx) {
                /* Bottom wall. */
                if (j == 0)
                    arrUV[j][i][iV] = 0.0;
                /* Top wall. */
                else if (j == info.my)
                    arrUV[j][i][iV] = 0.0;
                else
                    arrUV[j][i][iV] = arrUV_star[j][i][iV]
                                      - ns->dt/mat->rho * (arrp_prime[j][i] - arrp_prime[j-1][i]) / hy;
            }
        }

    PetscCall(DMDAVecRestoreArray(mesh->da, sol->u, &arru));
    PetscCall(DMDAVecRestoreArray(mesh->da, sol->v, &arrv));
    PetscCall(DMDAVecRestoreArray(mesh->da, sol->p, &arrp));
    PetscCall(DMStagVecRestoreArray(mesh->stag, sol->UVW, &arrUV));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->u_star, &arru_star));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->v_star, &arrv_star));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->p_prime, &arrp_prime));
    PetscCall(DMStagVecRestoreArrayRead(mesh->stag, sol->UVW_star, &arrUV_star));

    PetscCall(DMLocalToLocalBegin(mesh->da, sol->u, INSERT_VALUES, sol->u));
    PetscCall(DMLocalToLocalEnd(mesh->da, sol->u, INSERT_VALUES, sol->u));
    PetscCall(DMLocalToLocalBegin(mesh->da, sol->v, INSERT_VALUES, sol->v));
    PetscCall(DMLocalToLocalEnd(mesh->da, sol->v, INSERT_VALUES, sol->v));
    PetscCall(DMLocalToLocalBegin(mesh->da, sol->p, INSERT_VALUES, sol->p));
    PetscCall(DMLocalToLocalEnd(mesh->da, sol->p, INSERT_VALUES, sol->p));

    PetscCall(VecCopy(sol->Nu, sol->Nu_prev));
    PetscCall(VecCopy(sol->Nv, sol->Nv_prev));

    return 0;
}

PetscErrorCode ComputeRHSUstar2d(KSP ksp, Vec b, void *ctx) {
    FcNS ns = ctx;
    FcMesh mesh = ns->mesh;
    FcSolution sol = ns->sol;
    FcMaterial mat = ns->mat;

    DMDALocalInfo info;
    PetscReal hx, hy;
    PetscReal **arrb;
    const PetscReal **arrNu, **arrNu_prev, **arru, **arrv, **arrp;
    PetscReal up, uw, ue, us, un, pw, pe;
    PetscInt i, j;

    PetscCall(DMDAGetLocalInfo(mesh->da, &info));
    hx = 1.0 / info.mx;
    hy = 1.0 / info.my;

    PetscCall(DMDAVecGetArray(mesh->dau, b, &arrb));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->Nu, &arrNu));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->Nu_prev, &arrNu_prev));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->u, &arru));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->v, &arrv));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->p, &arrp));

    for (j = info.ys; j < info.ys + info.ym; j++)
        for (i = info.xs; i < info.xs + info.xm; i++) {
            up = arru[j][i];
            /* Left wall. */
            if (i == 0) {
                uw = -arru[j][i];
                pw = arrp[j][i];
            } else {
                uw = arru[j][i-1];
                pw = arrp[j][i-1];
            }
            /* Right wall. */
            if (i == info.mx - 1) {
                ue = -arru[j][i];
                pe = arrp[j][i];
            } else {
                ue = arru[j][i+1];
                pe = arrp[j][i+1];
            }
            /* Bottom wall. */
            if (j == 0)
                us = -arru[j][i];
            else
                us = arru[j-1][i];
            /* Top wall. */
            if (j == info.my - 1)
                un = 2.0 - arru[j][i];
            else
                un = arru[j+1][i];

            arrb[j][i] = up - 1.5*ns->dt * arrNu[j][i] + 0.5*ns->dt * arrNu_prev[j][i]
                         - ns->dt/mat->rho * (pe-pw)/(2.0*hx)
                         + mat->mu*ns->dt/(2.0*mat->rho) * ((uw-2.0*up+ue)/(hx*hx) + (us-2.0*up+un)/(hy*hy));
            if (j == info.my - 1)
                arrb[j][i] += mat->mu*ns->dt/(mat->rho*hy*hy);
        }

    PetscCall(DMDAVecRestoreArray(mesh->dau, b, &arrb));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->Nu, &arrNu));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->Nu_prev, &arrNu_prev));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->u, &arru));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->v, &arrv));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->p, &arrp));

    return 0;
}

PetscErrorCode ComputeRHSVstar2d(KSP ksp, Vec b, void *ctx) {
    FcNS ns = ctx;
    FcMesh mesh = ns->mesh;
    FcSolution sol = ns->sol;
    FcMaterial mat = ns->mat;

    DMDALocalInfo info;
    PetscReal hx, hy;
    PetscReal **arrb;
    const PetscReal **arrNv, **arrNv_prev, **arru, **arrv, **arrp;
    PetscReal vp, vw, ve, vs, vn, ps, pn;
    PetscInt i, j;

    PetscCall(DMDAGetLocalInfo(mesh->da, &info));
    hx = 1.0 / info.mx;
    hy = 1.0 / info.my;

    PetscCall(DMDAVecGetArray(mesh->dav, b, &arrb));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->Nv, &arrNv));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->Nv_prev, &arrNv_prev));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->u, &arru));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->v, &arrv));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->p, &arrp));

    for (j = info.ys; j < info.ys + info.ym; j++)
        for (i = info.xs; i < info.xs + info.xm; i++) {
            vp = arrv[j][i];
            /* Left wall. */
            if (i == 0)
                vw = -arrv[j][i];
            else
                vw = arrv[j][i-1];
            /* Right wall. */
            if (i == info.mx - 1)
                ve = -arrv[j][i];
            else
                ve = arrv[j][i+1];
            /* Bottom wall. */
            if (j == 0) {
                vs = -arrv[j][i];
                ps = arrp[j][i];
            } else {
                vs = arrv[j-1][i];
                ps = arrp[j-1][i];
            }
            /* Top wall. */
            if (j == info.my - 1) {
                vn = -arrv[j][i];
                pn = arrp[j][i];
            } else {
                vn = arrv[j+1][i];
                pn = arrp[j+1][i];
            }

            arrb[j][i] = vp - 1.5*ns->dt * arrNv[j][i] + 0.5*ns->dt * arrNv_prev[j][i]
                         - ns->dt/mat->rho * (pn-ps)/(2.0*hy)
                         + mat->mu*ns->dt/(2.0*mat->rho) * ((vw-2.0*vp+ve)/(hx*hx) + (vs-2.0*vp+vn)/(hy*hy));
        }

    PetscCall(DMDAVecRestoreArray(mesh->dav, b, &arrb));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->Nv, &arrNv));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->Nv_prev, &arrNv_prev));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->u, &arru));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->v, &arrv));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->p, &arrp));

    return 0;
}

PetscErrorCode ComputeRHSPprime2d(KSP ksp, Vec b, void *ctx) {
    FcNS ns = ctx;
    FcMesh mesh = ns->mesh;
    FcSolution sol = ns->sol;
    FcMaterial mat = ns->mat;

    DMDALocalInfo info;
    PetscReal hx, hy;
    PetscReal **arrb;
    const PetscReal ***arrUV_star;
    PetscInt iU, iV;
    MatNullSpace nullspace;
    PetscInt i, j;

    DMDAGetLocalInfo(mesh->da, &info);
    hx = 1.0 / info.mx;
    hy = 1.0 / info.my;

    PetscCall(DMDAVecGetArray(mesh->dap, b, &arrb));
    PetscCall(DMStagVecGetArrayRead(mesh->stag, sol->UVW_star, &arrUV_star));

    PetscCall(DMStagGetLocationSlot(mesh->stag, DMSTAG_LEFT, 0, &iU));
    PetscCall(DMStagGetLocationSlot(mesh->stag, DMSTAG_DOWN, 0, &iV));

    for (j = info.ys; j < info.ys + info.ym; j++)
        for (i = info.xs; i < info.xs + info.xm; i++)
            arrb[j][i] = -mat->rho/ns->dt * (hy * (arrUV_star[j][i+1][iU] - arrUV_star[j][i][iU])
                                             + hx * (arrUV_star[j+1][i][iV] - arrUV_star[j][i][iV]));

    PetscCall(DMDAVecRestoreArray(mesh->dap, b, &arrb));
    PetscCall(DMStagVecRestoreArrayRead(mesh->stag, sol->UVW_star, &arrUV_star));

    /* Remove null space. */
    PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullspace));
    PetscCall(MatNullSpaceRemove(nullspace, b));
    PetscCall(MatNullSpaceDestroy(&nullspace));

    return 0;
}

PetscErrorCode ComputeOperatorsUVstar2d(KSP ksp, Mat J, Mat Jpre, void *ctx) {
    FcNS ns = ctx;
    FcMaterial mat = ns->mat;

    DM da;
    DMDALocalInfo info;
    PetscReal hx, hy;
    MatStencil row, col[5];
    PetscReal v[5];
    PetscInt ncols;
    PetscInt i, j;

    PetscCall(KSPGetDM(ksp, &da));
    PetscCall(DMDAGetLocalInfo(da, &info));
    hx = 1.0 / info.mx;
    hy = 1.0 / info.my;

    for (j = info.ys; j < info.ys + info.ym; j++)
        for (i = info.xs; i < info.xs + info.xm; i++) {
            row.i = col[0].i = i;
            row.j = col[0].j = j;
            v[0] = 1.0 + 2.0*mat->mu*ns->dt/mat->rho * (1.0/(hx*hx) + 1.0/(hy*hy));
            v[1] = v[2] = v[3] = v[4] = 0.0;
            ncols = 1;

            /* Not left wall. */
            if (i != 0) {
                col[ncols].i = i - 1;
                col[ncols].j = j;
                v[0] -= mat->mu*ns->dt/(2.0*mat->rho*hx*hx);
                v[ncols] = -mat->mu*ns->dt/(2.0*mat->rho*hx*hx);
                ncols++;
            }
            /* Not right wall. */
            if (i != info.mx - 1) {
                col[ncols].i = i + 1;
                col[ncols].j = j;
                v[0] -= mat->mu*ns->dt/(2.0*mat->rho*hx*hx);
                v[ncols] = -mat->mu*ns->dt/(2.0*mat->rho*hx*hx);
                ncols++;
            }
            /* Not bottom wall. */
            if (j != 0) {
                col[ncols].i = i;
                col[ncols].j = j - 1;
                v[0] -= mat->mu*ns->dt/(2.0*mat->rho*hy*hy);
                v[ncols] = -mat->mu*ns->dt/(2.0*mat->rho*hy*hy);
                ncols++;
            }
            /* Not top wall. */
            if (j != info.my - 1) {
                col[ncols].i = i;
                col[ncols].j = j + 1;
                v[0] -= mat->mu*ns->dt/(2.0*mat->rho*hy*hy);
                v[ncols] = -mat->mu*ns->dt/(2.0*mat->rho*hy*hy);
                ncols++;
            }

            PetscCall(MatSetValuesStencil(Jpre, 1, &row, ncols, col, v, INSERT_VALUES));
        }

    PetscCall(MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY));

    return 0;
}

PetscErrorCode ComputeOperatorsPprime2d(KSP ksp, Mat J, Mat Jpre, void *ctx) {
    DM da;
    DMDALocalInfo info;
    PetscReal hx, hy;
    MatStencil row, col[5];
    PetscReal v[5];
    PetscInt ncols;
    MatNullSpace nullspace;
    PetscInt i, j;

    PetscCall(KSPGetDM(ksp, &da));
    PetscCall(DMDAGetLocalInfo(da, &info));
    hx = 1.0 / info.mx;
    hy = 1.0 / info.my;

    for (j = info.ys; j < info.ys + info.ym; j++)
        for (i = info.xs; i < info.xs + info.xm; i++) {
            row.i = col[0].i = i;
            row.j = col[0].j = j;
            v[0] = v[1] = v[2] = v[3] = v[4] = 0.0;
            ncols = 1;

            /* Not left wall. */
            if (i != 0) {
                col[ncols].i = i - 1;
                col[ncols].j = j;
                v[0] += hy / hx;
                v[ncols] -= hy / hx;
                ncols++;
            }
            /* Not right wall. */
            if (i != info.mx - 1) {
                col[ncols].i = i + 1;
                col[ncols].j = j;
                v[0] += hy / hx;
                v[ncols] -= hy / hx;
                ncols++;
            }
            /* Not bottom wall. */
            if (j != 0) {
                col[ncols].i = i;
                col[ncols].j = j - 1;
                v[0] += hx / hy;
                v[ncols] -= hx / hy;
                ncols++;
            }
            /* Not top wall. */
            if (j != info.my - 1) {
                col[ncols].i = i;
                col[ncols].j = j + 1;
                v[0] += hx / hy;
                v[ncols] -= hx / hy;
                ncols++;
            }

            PetscCall(MatSetValuesStencil(Jpre, 1, &row, ncols, col, v, INSERT_VALUES));
        }

    PetscCall(MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY));

    /* Remove null space. */
    PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullspace));
    PetscCall(MatSetNullSpace(J, nullspace));
    PetscCall(MatNullSpaceDestroy(&nullspace));

    return 0;
}

PetscErrorCode ComputeRHSUstar3d(KSP ksp, Vec b, void *ctx) {
    return 0;
}

PetscErrorCode ComputeRHSVstar3d(KSP ksp, Vec b, void *ctx) {
    return 0;
}

PetscErrorCode ComputeRHSWstar3d(KSP ksp, Vec b, void *ctx) {
    return 0;
}

PetscErrorCode ComputeRHSPprime3d(KSP ksp, Vec b, void *ctx) {
    return 0;
}

PetscErrorCode ComputeOperatorsUVWstar3d(KSP ksp, Mat J, Mat Jpre, void *ctx) {
    return 0;
}

PetscErrorCode ComputeOperatorsPprime3d(KSP ksp, Mat J, Mat Jpre, void *ctx) {
    return 0;
}

