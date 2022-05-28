#include <petscdmda.h>
#include <petscdmstag.h>
#include <petscoptions.h>
#include <petscsys.h>
#include "../inc/nsfsm.h"
#include "../inc/private/meshimpl.h"
#include "../inc/private/nsfsmimpl.h"
#include "../inc/private/solutionimpl.h"

static PetscErrorCode FcNSDestroy_FSM(FcNS *ns);
static PetscErrorCode FcNSSetFromOptions_FSM(FcNS ns, PetscOptionItems *PetscOptionsObject);
static PetscErrorCode FcNSSetUp_FSM(FcNS ns);
static PetscErrorCode FcNSSolve_FSM(FcNS ns);

static PetscErrorCode CreateDM(NS_FSM *fsm, FcMesh mesh);
static PetscErrorCode CreateKSP(NS_FSM *fsm, FcMesh mesh, FcNS ns);
static PetscErrorCode CreateVec(NS_FSM *fsm, FcMesh mesh);
static PetscErrorCode Update(FcNS ns);

static PetscErrorCode CalculateConvection(FcNS ns);
static PetscErrorCode CalculateIntermediateVelocity(FcNS ns);
static PetscErrorCode CalculatePressureCorrection(FcNS ns);

static PetscErrorCode ComputeRHSUstar2d(KSP ksp, Vec b, void *ctx);
static PetscErrorCode ComputeRHSVstar2d(KSP ksp, Vec b, void *ctx);
static PetscErrorCode ComputeRHSPprime2d(KSP ksp, Vec b, void *ctx);
static PetscErrorCode ComputeOperatorsUVstar2d(KSP ksp, Mat J, Mat Jpre, void *ctx);
static PetscErrorCode ComputeOperatorsPprime2d(KSP ksp, Mat J, Mat Jpre, void *ctx);

static PetscErrorCode ComputeRHSUstar3d(KSP ksp, Vec b, void *ctx);
static PetscErrorCode ComputeRHSVstar3d(KSP ksp, Vec b, void *ctx);
static PetscErrorCode ComputeRHSWstar3d(KSP ksp, Vec b, void *ctx);
static PetscErrorCode ComputeRHSPprime3d(KSP ksp, Vec b, void *ctx);
static PetscErrorCode ComputeOperatorsUVWstar3d(KSP ksp, Mat J, Mat Jpre, void *ctx);
static PetscErrorCode ComputeOperatorsPprime3d(KSP ksp, Mat J, Mat Jpre, void *ctx);

PetscErrorCode FcNSFSMCreate(FcMesh mesh, FcSolution sol, FcMaterial mat, FcNS *ns) {
    struct _FcNSOps ops;
    NS_FSM *fsm;

    /* Set operations. */
    ops.destroy = FcNSDestroy_FSM;
    ops.setfromoptions = FcNSSetFromOptions_FSM;
    ops.setup = FcNSSetUp_FSM;
    ops.solve = FcNSSolve_FSM;

    /* Create data. */
    PetscCall(PetscNew(&fsm));

    /* Create Navier-Stokes solver. */
    PetscCall(FcNSCreate(mesh, sol, mat, FC_NS_FSM, &ops, fsm, ns));

    /* Set data. */
    PetscCall(CreateDM(fsm, mesh));
    PetscCall(CreateKSP(fsm, mesh, *ns));
    PetscCall(CreateVec(fsm, mesh));

    return 0;
}

static PetscErrorCode FcNSDestroy_FSM(FcNS *ns) {
    NS_FSM *fsm;
    FcMesh mesh;

    fsm = (NS_FSM *)((*ns)->data);
    mesh = (*ns)->mesh;

    PetscCall(DMDestroy(&fsm->dau));
    PetscCall(DMDestroy(&fsm->dav));
    if (mesh->dim == 3)
        PetscCall(DMDestroy(&fsm->daw));
    PetscCall(DMDestroy(&fsm->dap));

    PetscCall(KSPDestroy(&fsm->kspu));
    PetscCall(KSPDestroy(&fsm->kspv));
    if (mesh->dim == 3)
        PetscCall(KSPDestroy(&fsm->kspw));
    PetscCall(KSPDestroy(&fsm->kspp));

    PetscCall(VecDestroy(&fsm->p));
    PetscCall(VecDestroy(&fsm->UVW));
    PetscCall(VecDestroy(&fsm->u_star));
    PetscCall(VecDestroy(&fsm->v_star));
    if (mesh->dim == 3)
        PetscCall(VecDestroy(&fsm->w_star));
    PetscCall(VecDestroy(&fsm->UVW_star));
    PetscCall(VecDestroy(&fsm->p_prime));
    PetscCall(VecDestroy(&fsm->Nu));
    PetscCall(VecDestroy(&fsm->Nv));
    if (mesh->dim == 3)
        PetscCall(VecDestroy(&fsm->Nw));
    PetscCall(VecDestroy(&fsm->p_prev));
    PetscCall(VecDestroy(&fsm->Nu_prev));
    PetscCall(VecDestroy(&fsm->Nv_prev));
    if (mesh->dim == 3)
        PetscCall(VecDestroy(&fsm->Nw_prev));
    PetscCall(VecDestroy(&fsm->u_tilde));
    PetscCall(VecDestroy(&fsm->v_tilde));
    if (mesh->dim == 3)
        PetscCall(VecDestroy(&fsm->w_tilde));

    PetscCall(PetscFree(fsm));

    return 0;
}

static PetscErrorCode FcNSSetFromOptions_FSM(FcNS ns, PetscOptionItems *PetscOptionsObject) {
    PetscOptionsHead(PetscOptionsObject, "Fractional step method options");

    PetscOptionsTail();

    return 0;
}

static PetscErrorCode FcNSSetUp_FSM(FcNS ns) {
    if (ns->setup)
        return 0;

    return 0;
}

static PetscErrorCode FcNSSolve_FSM(FcNS ns) {
    PetscInt i;

    PetscCall(FcNSSetUp_FSM(ns));

    for (i = 0; i < ns->maxsteps; i++) {
        PetscCall(CalculateConvection(ns));
        PetscCall(CalculateIntermediateVelocity(ns));
        PetscCall(CalculatePressureCorrection(ns));
        PetscCall(Update(ns));
    }

    return 0;
}

static PetscErrorCode CreateDM(NS_FSM *fsm, FcMesh mesh) {
    PetscCall(DMClone(mesh->da, &fsm->dau));
    PetscCall(DMClone(mesh->da, &fsm->dav));
    if (mesh->dim == 3)
        PetscCall(DMClone(mesh->da, &fsm->daw));
    PetscCall(DMClone(mesh->da, &fsm->dap));

    return 0;
}

static PetscErrorCode CreateKSP(NS_FSM *fsm, FcMesh mesh, FcNS ns) {
    PetscCall(KSPCreate(mesh->obj.comm, &fsm->kspu));
    PetscCall(KSPSetDM(fsm->kspu, fsm->dau));
    if (mesh->dim == 2) {
        PetscCall(KSPSetComputeRHS(fsm->kspu, ComputeRHSUstar2d, ns));
        PetscCall(KSPSetComputeOperators(fsm->kspu, ComputeOperatorsUVstar2d, ns));
    } else {
        PetscCall(KSPSetComputeRHS(fsm->kspu, ComputeRHSUstar3d, ns));
        PetscCall(KSPSetComputeOperators(fsm->kspu, ComputeOperatorsUVWstar3d, ns));
    }
    PetscCall(KSPSetFromOptions(fsm->kspu));

    PetscCall(KSPCreate(mesh->obj.comm, &fsm->kspv));
    PetscCall(KSPSetDM(fsm->kspv, fsm->dav));
    if (mesh->dim == 2) {
        PetscCall(KSPSetComputeRHS(fsm->kspv, ComputeRHSVstar2d, ns));
        PetscCall(KSPSetComputeOperators(fsm->kspv, ComputeOperatorsUVstar2d, ns));
    } else {
        PetscCall(KSPSetComputeRHS(fsm->kspv, ComputeRHSVstar3d, ns));
        PetscCall(KSPSetComputeOperators(fsm->kspv, ComputeOperatorsUVWstar3d, ns));
    }
    PetscCall(KSPSetFromOptions(fsm->kspv));

    if (mesh->dim == 3) {
        PetscCall(KSPCreate(mesh->obj.comm, &fsm->kspw));
        PetscCall(KSPSetDM(fsm->kspw, fsm->daw));
        PetscCall(KSPSetComputeRHS(fsm->kspw, ComputeRHSWstar3d, ns));
        PetscCall(KSPSetComputeOperators(fsm->kspw, ComputeOperatorsUVWstar3d, ns));
        PetscCall(KSPSetFromOptions(fsm->kspw));
    }

    PetscCall(KSPCreate(mesh->obj.comm, &fsm->kspp));
    PetscCall(KSPSetDM(fsm->kspp, fsm->dap));
    if (mesh->dim == 2) {
        PetscCall(KSPSetComputeRHS(fsm->kspp, ComputeRHSPprime2d, ns));
        PetscCall(KSPSetComputeOperators(fsm->kspp, ComputeOperatorsPprime2d, ns));
    } else {
        PetscCall(KSPSetComputeRHS(fsm->kspp, ComputeRHSPprime3d, ns));
        PetscCall(KSPSetComputeOperators(fsm->kspp, ComputeOperatorsPprime3d, ns));
    }
    PetscCall(KSPSetFromOptions(fsm->kspp));

    return 0;
}

static PetscErrorCode CreateVec(NS_FSM *fsm, FcMesh mesh) {
    PetscCall(DMCreateLocalVector(mesh->da, &fsm->p));
    PetscCall(DMCreateLocalVector(mesh->stag, &fsm->UVW));
    PetscCall(DMCreateLocalVector(mesh->da, &fsm->u_star));
    PetscCall(DMCreateLocalVector(mesh->da, &fsm->v_star));
    if (mesh->dim == 3)
        PetscCall(DMCreateLocalVector(mesh->da, &fsm->w_star));
    PetscCall(DMCreateLocalVector(mesh->stag, &fsm->UVW_star));
    PetscCall(DMCreateLocalVector(mesh->da, &fsm->p_prime));
    PetscCall(DMCreateLocalVector(mesh->da, &fsm->Nu));
    PetscCall(DMCreateLocalVector(mesh->da, &fsm->Nv));
    if (mesh->dim == 3)
        PetscCall(DMCreateLocalVector(mesh->da, &fsm->Nw));
    PetscCall(DMCreateLocalVector(mesh->da, &fsm->p_prev));
    PetscCall(DMCreateLocalVector(mesh->da, &fsm->Nu_prev));
    PetscCall(DMCreateLocalVector(mesh->da, &fsm->Nv_prev));
    if (mesh->dim == 3)
        PetscCall(DMCreateLocalVector(mesh->da, &fsm->Nw_prev));
    PetscCall(DMCreateLocalVector(mesh->da, &fsm->u_tilde));
    PetscCall(DMCreateLocalVector(mesh->da, &fsm->v_tilde));
    if (mesh->dim == 3)
        PetscCall(DMCreateLocalVector(mesh->da, &fsm->w_tilde));

    return 0;
}

static PetscErrorCode CalculateConvection(FcNS ns) {
    NS_FSM *fsm = ns->data;
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

    PetscCall(DMDAVecGetArray(mesh->da, fsm->Nu, &arrNu));
    PetscCall(DMDAVecGetArray(mesh->da, fsm->Nv, &arrNv));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->u, &arru));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->v, &arrv));
    PetscCall(DMStagVecGetArrayRead(mesh->stag, fsm->UVW, &arrUV));

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

    PetscCall(DMDAVecRestoreArray(mesh->da, fsm->Nu, &arrNu));
    PetscCall(DMDAVecRestoreArray(mesh->da, fsm->Nv, &arrNv));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->u, &arru));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->v, &arrv));
    PetscCall(DMStagVecRestoreArrayRead(mesh->stag, fsm->UVW, &arrUV));

    return 0;
}

static PetscErrorCode CalculateIntermediateVelocity(FcNS ns) {
    NS_FSM *fsm = ns->data;
    FcMesh mesh = ns->mesh;
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
    PetscCall(KSPSolve(fsm->kspu, NULL, NULL));
    PetscCall(KSPGetSolution(fsm->kspu, &x));
    PetscCall(DMGlobalToLocal(fsm->dau, x, INSERT_VALUES, fsm->u_star));

    PetscCall(KSPSolve(fsm->kspv, NULL, NULL));
    PetscCall(KSPGetSolution(fsm->kspv, &x));
    PetscCall(DMGlobalToLocal(fsm->dav, x, INSERT_VALUES, fsm->v_star));

    /* Calculate face-centered intermediate velocity. */
    PetscCall(DMDAGetLocalInfo(mesh->da, &info));
    hx = 1.0 / info.mx;
    hy = 1.0 / info.my;

    PetscCall(DMDAVecGetArray(mesh->da, fsm->u_tilde, &arru_tilde));
    PetscCall(DMDAVecGetArray(mesh->da, fsm->v_tilde, &arrv_tilde));
    PetscCall(DMDAVecGetArrayRead(mesh->da, fsm->u_star, &arru_star));
    PetscCall(DMDAVecGetArrayRead(mesh->da, fsm->v_star, &arrv_star));
    PetscCall(DMDAVecGetArrayRead(mesh->da, fsm->p, &arrp));

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

            arru_tilde[j][i] = arru_star[j][i] + ns->dt/mat.rho * (pe - pw) / (2.0 * hx);
            arrv_tilde[j][i] = arrv_star[j][i] + ns->dt/mat.rho * (pn - ps) / (2.0 * hy);
        }

    PetscCall(DMDAVecRestoreArray(mesh->da, fsm->u_tilde, &arru_tilde));
    PetscCall(DMDAVecRestoreArray(mesh->da, fsm->v_tilde, &arrv_tilde));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, fsm->u_star, &arru_star));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, fsm->v_star, &arrv_star));

    PetscCall(DMLocalToLocalBegin(mesh->da, fsm->u_tilde, INSERT_VALUES, fsm->u_tilde));
    PetscCall(DMLocalToLocalEnd(mesh->da, fsm->u_tilde, INSERT_VALUES, fsm->u_tilde));
    PetscCall(DMLocalToLocalBegin(mesh->da, fsm->v_tilde, INSERT_VALUES, fsm->v_tilde));
    PetscCall(DMLocalToLocalEnd(mesh->da, fsm->v_tilde, INSERT_VALUES, fsm->v_tilde));

    PetscCall(DMStagVecGetArray(mesh->stag, fsm->UVW_star, &arrUV_star));
    PetscCall(DMDAVecGetArrayRead(mesh->da, fsm->u_tilde, &arru_tilde_rd));
    PetscCall(DMDAVecGetArrayRead(mesh->da, fsm->v_tilde, &arrv_tilde_rd));

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
                                           - ns->dt/mat.rho * (arrp[j][i] - arrp[j][i-1]) / hx;
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
                                           - ns->dt/mat.rho * (arrp[j][i] - arrp[j-1][i]) / hy;
            }
        }

    PetscCall(DMStagVecRestoreArray(mesh->stag, fsm->UVW_star, &arrUV_star));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, fsm->u_tilde, &arru_tilde_rd));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, fsm->v_tilde, &arrv_tilde_rd));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, fsm->p, &arrp));

    return 0;
}

static PetscErrorCode CalculatePressureCorrection(FcNS ns) {
    NS_FSM *fsm = ns->data;

    Vec x;

    PetscCall(KSPSolve(fsm->kspp, NULL, NULL));
    PetscCall(KSPGetSolution(fsm->kspp, &x));
    PetscCall(DMGlobalToLocal(fsm->dap, x, INSERT_VALUES, fsm->p_prime));

    return 0;
}

static PetscErrorCode Update(FcNS ns) {
    NS_FSM *fsm = ns->data;
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
    PetscCall(DMDAVecGetArray(mesh->da, fsm->p, &arrp));
    PetscCall(DMStagVecGetArray(mesh->stag, fsm->UVW, &arrUV));
    PetscCall(DMDAVecGetArrayRead(mesh->da, fsm->u_star, &arru_star));
    PetscCall(DMDAVecGetArrayRead(mesh->da, fsm->v_star, &arrv_star));
    PetscCall(DMDAVecGetArrayRead(mesh->da, fsm->p_prime, &arrp_prime));
    PetscCall(DMStagVecGetArrayRead(mesh->stag, fsm->UVW_star, &arrUV_star));

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

            arru[j][i] = arru_star[j][i] - ns->dt/mat.rho * (ppe - ppw) / (2.0*hx);
            arrv[j][i] = arrv_star[j][i] - ns->dt/mat.rho * (ppn - pps) / (2.0*hy);
            arrp[j][i] += arrp_prime[j][i]
                          - mat.mu*ns->dt/(2.0*mat.rho) * ((ppe - 2.0*ppp + ppw) / (hx*hx)
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
                                      - ns->dt/mat.rho * (arrp_prime[j][i] - arrp_prime[j][i-1]) / hx;
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
                                      - ns->dt/mat.rho * (arrp_prime[j][i] - arrp_prime[j-1][i]) / hy;
            }
        }

    PetscCall(DMDAVecRestoreArray(mesh->da, sol->u, &arru));
    PetscCall(DMDAVecRestoreArray(mesh->da, sol->v, &arrv));
    PetscCall(DMDAVecRestoreArray(mesh->da, fsm->p, &arrp));
    PetscCall(DMStagVecRestoreArray(mesh->stag, fsm->UVW, &arrUV));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, fsm->u_star, &arru_star));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, fsm->v_star, &arrv_star));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, fsm->p_prime, &arrp_prime));
    PetscCall(DMStagVecRestoreArrayRead(mesh->stag, fsm->UVW_star, &arrUV_star));

    PetscCall(DMLocalToLocalBegin(mesh->da, sol->u, INSERT_VALUES, sol->u));
    PetscCall(DMLocalToLocalEnd(mesh->da, sol->u, INSERT_VALUES, sol->u));
    PetscCall(DMLocalToLocalBegin(mesh->da, sol->v, INSERT_VALUES, sol->v));
    PetscCall(DMLocalToLocalEnd(mesh->da, sol->v, INSERT_VALUES, sol->v));
    PetscCall(DMLocalToLocalBegin(mesh->da, fsm->p, INSERT_VALUES, fsm->p));
    PetscCall(DMLocalToLocalEnd(mesh->da, fsm->p, INSERT_VALUES, fsm->p));

    PetscCall(VecCopy(fsm->Nu, fsm->Nu_prev));
    PetscCall(VecCopy(fsm->Nv, fsm->Nv_prev));

    return 0;
}

static PetscErrorCode ComputeRHSUstar2d(KSP ksp, Vec b, void *ctx) {
    FcNS ns = ctx;
    FcMesh mesh = ns->mesh;
    FcSolution sol = ns->sol;
    NS_FSM *fsm = ns->data;
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

    PetscCall(DMDAVecGetArray(fsm->dau, b, &arrb));
    PetscCall(DMDAVecGetArrayRead(mesh->da, fsm->Nu, &arrNu));
    PetscCall(DMDAVecGetArrayRead(mesh->da, fsm->Nu_prev, &arrNu_prev));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->u, &arru));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->v, &arrv));
    PetscCall(DMDAVecGetArrayRead(mesh->da, fsm->p, &arrp));

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
                         - ns->dt/mat.rho * (pe-pw)/(2.0*hx)
                         + mat.mu*ns->dt/(2.0*mat.rho) * ((uw-2.0*up+ue)/(hx*hx) + (us-2.0*up+un)/(hy*hy));
            if (j == info.my - 1)
                arrb[j][i] += mat.mu*ns->dt/(mat.rho*hy*hy);
        }

    PetscCall(DMDAVecRestoreArray(fsm->dau, b, &arrb));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, fsm->Nu, &arrNu));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, fsm->Nu_prev, &arrNu_prev));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->u, &arru));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->v, &arrv));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, fsm->p, &arrp));

    return 0;
}

static PetscErrorCode ComputeRHSVstar2d(KSP ksp, Vec b, void *ctx) {
    FcNS ns = ctx;
    FcMesh mesh = ns->mesh;
    FcSolution sol = ns->sol;
    NS_FSM *fsm = ns->data;
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

    PetscCall(DMDAVecGetArray(fsm->dav, b, &arrb));
    PetscCall(DMDAVecGetArrayRead(mesh->da, fsm->Nv, &arrNv));
    PetscCall(DMDAVecGetArrayRead(mesh->da, fsm->Nv_prev, &arrNv_prev));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->u, &arru));
    PetscCall(DMDAVecGetArrayRead(mesh->da, sol->v, &arrv));
    PetscCall(DMDAVecGetArrayRead(mesh->da, fsm->p, &arrp));

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
                         - ns->dt/mat.rho * (pn-ps)/(2.0*hy)
                         + mat.mu*ns->dt/(2.0*mat.rho) * ((vw-2.0*vp+ve)/(hx*hx) + (vs-2.0*vp+vn)/(hy*hy));
        }

    PetscCall(DMDAVecRestoreArray(fsm->dav, b, &arrb));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, fsm->Nv, &arrNv));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, fsm->Nv_prev, &arrNv_prev));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->u, &arru));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, sol->v, &arrv));
    PetscCall(DMDAVecRestoreArrayRead(mesh->da, fsm->p, &arrp));

    return 0;
}

static PetscErrorCode ComputeRHSPprime2d(KSP ksp, Vec b, void *ctx) {
    FcNS ns = ctx;
    FcMesh mesh = ns->mesh;
    NS_FSM *fsm = ns->data;
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

    PetscCall(DMDAVecGetArray(fsm->dap, b, &arrb));
    PetscCall(DMStagVecGetArrayRead(mesh->stag, fsm->UVW_star, &arrUV_star));

    PetscCall(DMStagGetLocationSlot(mesh->stag, DMSTAG_LEFT, 0, &iU));
    PetscCall(DMStagGetLocationSlot(mesh->stag, DMSTAG_DOWN, 0, &iV));

    for (j = info.ys; j < info.ys + info.ym; j++)
        for (i = info.xs; i < info.xs + info.xm; i++)
            arrb[j][i] = -mat.rho/ns->dt * (hy * (arrUV_star[j][i+1][iU] - arrUV_star[j][i][iU])
                                            + hx * (arrUV_star[j+1][i][iV] - arrUV_star[j][i][iV]));

    PetscCall(DMDAVecRestoreArray(fsm->dap, b, &arrb));
    PetscCall(DMStagVecRestoreArrayRead(mesh->stag, fsm->UVW_star, &arrUV_star));

    /* Remove null space. */
    PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullspace));
    PetscCall(MatNullSpaceRemove(nullspace, b));
    PetscCall(MatNullSpaceDestroy(&nullspace));

    return 0;
}

static PetscErrorCode ComputeOperatorsUVstar2d(KSP ksp, Mat J, Mat Jpre, void *ctx) {
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
            v[0] = 1.0 + 2.0*mat.mu*ns->dt/mat.rho * (1.0/(hx*hx) + 1.0/(hy*hy));
            v[1] = v[2] = v[3] = v[4] = 0.0;
            ncols = 1;

            /* Not left wall. */
            if (i != 0) {
                col[ncols].i = i - 1;
                col[ncols].j = j;
                v[0] -= mat.mu*ns->dt/(2.0*mat.rho*hx*hx);
                v[ncols] = -mat.mu*ns->dt/(2.0*mat.rho*hx*hx);
                ncols++;
            }
            /* Not right wall. */
            if (i != info.mx - 1) {
                col[ncols].i = i + 1;
                col[ncols].j = j;
                v[0] -= mat.mu*ns->dt/(2.0*mat.rho*hx*hx);
                v[ncols] = -mat.mu*ns->dt/(2.0*mat.rho*hx*hx);
                ncols++;
            }
            /* Not bottom wall. */
            if (j != 0) {
                col[ncols].i = i;
                col[ncols].j = j - 1;
                v[0] -= mat.mu*ns->dt/(2.0*mat.rho*hy*hy);
                v[ncols] = -mat.mu*ns->dt/(2.0*mat.rho*hy*hy);
                ncols++;
            }
            /* Not top wall. */
            if (j != info.my - 1) {
                col[ncols].i = i;
                col[ncols].j = j + 1;
                v[0] -= mat.mu*ns->dt/(2.0*mat.rho*hy*hy);
                v[ncols] = -mat.mu*ns->dt/(2.0*mat.rho*hy*hy);
                ncols++;
            }

            PetscCall(MatSetValuesStencil(Jpre, 1, &row, ncols, col, v, INSERT_VALUES));
        }

    PetscCall(MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY));

    return 0;
}

static PetscErrorCode ComputeOperatorsPprime2d(KSP ksp, Mat J, Mat Jpre, void *ctx) {
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

static PetscErrorCode ComputeRHSUstar3d(KSP ksp, Vec b, void *ctx) {
    return 0;
}

static PetscErrorCode ComputeRHSVstar3d(KSP ksp, Vec b, void *ctx) {
    return 0;
}

static PetscErrorCode ComputeRHSWstar3d(KSP ksp, Vec b, void *ctx) {
    return 0;
}

static PetscErrorCode ComputeRHSPprime3d(KSP ksp, Vec b, void *ctx) {
    return 0;
}

static PetscErrorCode ComputeOperatorsUVWstar3d(KSP ksp, Mat J, Mat Jpre, void *ctx) {
    return 0;
}

static PetscErrorCode ComputeOperatorsPprime3d(KSP ksp, Mat J, Mat Jpre, void *ctx) {
    return 0;
}

