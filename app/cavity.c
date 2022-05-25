#include <petsc.h>

const char *const help = "Solves a 2D lid-driven cavity flow.\n\n";

DM da, dau, dav, dap, stag;
DMDALocalInfo info;
Vec u, v, p, UV, u_star, v_star, u_tilde, v_tilde, p_prime, UV_star;
Vec Nu, Nv, Nu_prev, Nv_prev;

KSP kspu, kspv, kspp;

PetscReal rho = 1.0, mu = 1.0, dt = 0.01;
PetscInt nt = 1000;

PetscErrorCode CreateDM(void);
PetscErrorCode CreateVec(void);
PetscErrorCode CreateSolver(void);
PetscErrorCode CalculateConvection(void);
PetscErrorCode CalculateIntermediateVelocity(void);
PetscErrorCode CalculatePressureCorrection(void);
PetscErrorCode Update(void);

PetscErrorCode ComputeRHSUstar(KSP ksp, Vec b, void *ctx);
PetscErrorCode ComputeRHSVstar(KSP ksp, Vec b, void *ctx);
PetscErrorCode ComputeRHSPprime(KSP ksp, Vec b, void *ctx);
PetscErrorCode ComputeOperatorsUVstar(KSP ksp, Mat J, Mat Jpre, void *ctx);
PetscErrorCode ComputeOperatorsPprime(KSP ksp, Mat J, Mat Jpre, void *ctx);

int main(int argc, char *argv[]) {
    /* Initialize PETSc. */
    PetscCall(PetscInitialize(&argc, &argv, NULL, help));

    /* Options. */
    PetscOptionsBegin(PETSC_COMM_WORLD, "cavity_", "Cavity flow options", "cavity");
    PetscOptionsReal("-rho", "Density", NULL, rho, &rho, NULL);
    PetscOptionsReal("-mu", "Viscosity", NULL, mu, &mu, NULL);
    PetscOptionsReal("-dt", "Time step size", NULL, dt, &dt, NULL);
    PetscOptionsInt("-nt", "Number of time steps", NULL, nt, &nt, NULL);
    PetscOptionsEnd();

    /* Create DMs. */
    PetscCall(CreateDM());
    /* Create vectors. */
    PetscCall(CreateVec());
    /* Create solvers. */
    PetscCall(CreateSolver());

    for (PetscInt t = 0; t < nt; t++) {
        /* Calculate convection term. */
        PetscCall(CalculateConvection());
        /* Calculate intermediate velocity. */
        PetscCall(CalculateIntermediateVelocity());
        /* Calculate pressure correction. */
        PetscCall(CalculatePressureCorrection());
        /* Update to next time step. */
        PetscCall(Update());
    }

    FILE *fp;

    fp = fopen("u.txt", "w");
    if (fp) {
        const PetscReal **arru;
        PetscInt i, j;

        PetscCall(DMDAVecGetArrayRead(da, u, &arru));
        for (i = 0; i < info.mx; i++) {
            for (j = 0; j < info.my; j++) {
                PetscFPrintf(PETSC_COMM_WORLD, fp, "%f ", arru[j][i]);
            }
            PetscFPrintf(PETSC_COMM_WORLD, fp, "\n");
        }
        PetscCall(DMDAVecRestoreArrayRead(da, u, &arru));

        fclose(fp);
    }

    /* Destroy. */
    PetscCall(DMDestroy(&da));
    PetscCall(DMDestroy(&dau));
    PetscCall(DMDestroy(&dav));
    PetscCall(DMDestroy(&dap));
    PetscCall(DMDestroy(&stag));

    PetscCall(VecDestroy(&u));
    PetscCall(VecDestroy(&v));
    PetscCall(VecDestroy(&p));
    PetscCall(VecDestroy(&UV));
    PetscCall(VecDestroy(&u_star));
    PetscCall(VecDestroy(&v_star));
    PetscCall(VecDestroy(&u_tilde));
    PetscCall(VecDestroy(&v_tilde));
    PetscCall(VecDestroy(&p_prime));
    PetscCall(VecDestroy(&Nu));
    PetscCall(VecDestroy(&Nv));
    PetscCall(VecDestroy(&Nu_prev));
    PetscCall(VecDestroy(&Nv_prev));

    PetscCall(KSPDestroy(&kspu));
    PetscCall(KSPDestroy(&kspv));
    PetscCall(KSPDestroy(&kspp));

    PetscCall(PetscFinalize());

    return 0;
}

PetscErrorCode CreateDM(void) {
    PetscInt M, N, m, n;
    const PetscInt *lx, *ly;

    /* Create DMDA. */
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, 3, 3, PETSC_DECIDE,
                           PETSC_DECIDE, 1, 1, NULL, NULL, &da));
    PetscCall(DMSetFromOptions(da));
    PetscCall(DMSetUp(da));

    PetscCall(DMDAGetInfo(da, NULL, &M, &N, NULL, &m, &n, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
    PetscCall(DMDAGetOwnershipRanges(da, &lx, &ly, NULL));
    PetscCall(DMDAGetLocalInfo(da, &info));

    PetscCall(DMClone(da, &dau));
    PetscCall(DMClone(da, &dav));
    PetscCall(DMClone(da, &dap));

    /* Create DMStag. */
    PetscCall(DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, M, N, m, n, 0, 1, 0,
                             DMSTAG_STENCIL_STAR, 1, lx, ly, &stag));
    PetscCall(DMSetUp(stag));

    return 0;
}

PetscErrorCode CreateVec(void) {
    PetscCall(DMCreateLocalVector(da, &u));
    PetscCall(DMCreateLocalVector(da, &v));
    PetscCall(DMCreateLocalVector(da, &p));
    PetscCall(DMCreateLocalVector(stag, &UV));
    PetscCall(DMCreateLocalVector(da, &u_star));
    PetscCall(DMCreateLocalVector(da, &v_star));
    PetscCall(DMCreateLocalVector(da, &u_tilde));
    PetscCall(DMCreateLocalVector(da, &v_tilde));
    PetscCall(DMCreateLocalVector(da, &p_prime));
    PetscCall(DMCreateLocalVector(stag, &UV_star));

    PetscCall(DMCreateLocalVector(da, &Nu));
    PetscCall(DMCreateLocalVector(da, &Nv));
    PetscCall(DMCreateLocalVector(da, &Nu_prev));
    PetscCall(DMCreateLocalVector(da, &Nv_prev));

    return 0;
}

PetscErrorCode CreateSolver(void) {
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &kspu));
    PetscCall(KSPSetDM(kspu, dau));
    PetscCall(KSPSetComputeRHS(kspu, ComputeRHSUstar, NULL));
    PetscCall(KSPSetComputeOperators(kspu, ComputeOperatorsUVstar, NULL));
    PetscCall(KSPSetFromOptions(kspu));

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &kspv));
    PetscCall(KSPSetDM(kspv, dav));
    PetscCall(KSPSetComputeRHS(kspv, ComputeRHSVstar, NULL));
    PetscCall(KSPSetComputeOperators(kspv, ComputeOperatorsUVstar, NULL));
    PetscCall(KSPSetFromOptions(kspv));

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &kspp));
    PetscCall(KSPSetDM(kspp, dap));
    PetscCall(KSPSetComputeRHS(kspp, ComputeRHSPprime, NULL));
    PetscCall(KSPSetComputeOperators(kspp, ComputeOperatorsPprime, NULL));
    PetscCall(KSPSetFromOptions(kspp));

    return 0;
}

PetscErrorCode CalculateConvection(void) {
    PetscReal hx, hy;
    PetscReal **arrNu, **arrNv;
    const PetscReal **arru, **arrv, ***arrUV;
    PetscInt iU, iV;
    PetscReal up, uw, ue, us, un, vp, vw, ve, vs, vn;
    PetscInt i, j;

    hx = 1.0 / info.mx;
    hy = 1.0 / info.my;

    PetscCall(DMDAVecGetArray(da, Nu, &arrNu));
    PetscCall(DMDAVecGetArray(da, Nv, &arrNv));
    PetscCall(DMDAVecGetArrayRead(da, u, &arru));
    PetscCall(DMDAVecGetArrayRead(da, v, &arrv));
    PetscCall(DMStagVecGetArrayRead(stag, UV, &arrUV));

    PetscCall(DMStagGetLocationSlot(stag, DMSTAG_LEFT, 0, &iU));
    PetscCall(DMStagGetLocationSlot(stag, DMSTAG_DOWN, 0, &iV));

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

    PetscCall(DMDAVecRestoreArray(da, Nu, &arrNu));
    PetscCall(DMDAVecRestoreArray(da, Nv, &arrNv));
    PetscCall(DMDAVecRestoreArrayRead(da, u, &arru));
    PetscCall(DMDAVecRestoreArrayRead(da, v, &arrv));
    PetscCall(DMStagVecRestoreArrayRead(stag, UV, &arrUV));

    return 0;
}

PetscErrorCode CalculateIntermediateVelocity(void) {
    Vec x;
    PetscReal hx, hy;
    PetscReal **arru_tilde, **arrv_tilde, ***arrUV_star;
    const PetscReal **arru_star, **arrv_star, **arrp, **arru_tilde_rd, **arrv_tilde_rd;
    PetscInt iU, iV;
    PetscReal pw, pe, ps, pn;
    PetscInt i, j;

    /* Calculate cell-centered intermediate velocity. */
    PetscCall(KSPSolve(kspu, NULL, NULL));
    PetscCall(KSPGetSolution(kspu, &x));
    PetscCall(DMGlobalToLocal(dau, x, INSERT_VALUES, u_star));

    PetscCall(KSPSolve(kspv, NULL, NULL));
    PetscCall(KSPGetSolution(kspv, &x));
    PetscCall(DMGlobalToLocal(dav, x, INSERT_VALUES, v_star));

    /* Calculate face-centered intermediate velocity. */
    hx = 1.0 / info.mx;
    hy = 1.0 / info.my;

    PetscCall(DMDAVecGetArray(da, u_tilde, &arru_tilde));
    PetscCall(DMDAVecGetArray(da, v_tilde, &arrv_tilde));
    PetscCall(DMDAVecGetArrayRead(da, u_star, &arru_star));
    PetscCall(DMDAVecGetArrayRead(da, v_star, &arrv_star));
    PetscCall(DMDAVecGetArrayRead(da, p, &arrp));

    PetscCall(DMStagGetLocationSlot(stag, DMSTAG_LEFT, 0, &iU));
    PetscCall(DMStagGetLocationSlot(stag, DMSTAG_DOWN, 0, &iV));

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

            arru_tilde[j][i] = arru_star[j][i] + dt/rho * (pe - pw) / (2.0 * hx);
            arrv_tilde[j][i] = arrv_star[j][i] + dt/rho * (pn - ps) / (2.0 * hy);
        }

    PetscCall(DMDAVecRestoreArray(da, u_tilde, &arru_tilde));
    PetscCall(DMDAVecRestoreArray(da, v_tilde, &arrv_tilde));
    PetscCall(DMDAVecRestoreArrayRead(da, u_star, &arru_star));
    PetscCall(DMDAVecRestoreArrayRead(da, v_star, &arrv_star));

    DMLocalToLocalBegin(da, u_tilde, INSERT_VALUES, u_tilde);
    DMLocalToLocalEnd(da, u_tilde, INSERT_VALUES, u_tilde);
    DMLocalToLocalBegin(da, v_tilde, INSERT_VALUES, v_tilde);
    DMLocalToLocalEnd(da, v_tilde, INSERT_VALUES, v_tilde);

    PetscCall(DMStagVecGetArray(stag, UV_star, &arrUV_star));
    PetscCall(DMDAVecGetArrayRead(da, u_tilde, &arru_tilde_rd));
    PetscCall(DMDAVecGetArrayRead(da, v_tilde, &arrv_tilde_rd));

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
                                           - dt/rho * (arrp[j][i] - arrp[j][i-1]) / hx;
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
                                           - dt/rho * (arrp[j][i] - arrp[j-1][i]) / hy;
            }
        }

    PetscCall(DMStagVecRestoreArray(stag, UV_star, &arrUV_star));
    PetscCall(DMDAVecRestoreArrayRead(da, u_tilde, &arru_tilde_rd));
    PetscCall(DMDAVecRestoreArrayRead(da, v_tilde, &arrv_tilde_rd));
    PetscCall(DMDAVecRestoreArrayRead(da, p, &arrp));

    return 0;
}

PetscErrorCode CalculatePressureCorrection(void) {
    Vec x;

    PetscCall(KSPSolve(kspp, NULL, NULL));
    PetscCall(KSPGetSolution(kspp, &x));
    PetscCall(DMGlobalToLocal(dap, x, INSERT_VALUES, p_prime));

    return 0;
}

PetscErrorCode Update(void) {
    PetscReal hx, hy;
    PetscReal **arru, **arrv, **arrp, ***arrUV;
    const PetscReal **arru_star, **arrv_star, **arrp_prime, ***arrUV_star;
    PetscInt iU, iV;
    PetscReal ppp, ppw, ppe, pps, ppn, pw, pe, ps, pn;
    PetscInt i, j;

    hx = 1.0 / info.mx;
    hy = 1.0 / info.my;

    PetscCall(DMDAVecGetArray(da, u, &arru));
    PetscCall(DMDAVecGetArray(da, v, &arrv));
    PetscCall(DMDAVecGetArray(da, p, &arrp));
    PetscCall(DMStagVecGetArray(stag, UV, &arrUV));
    PetscCall(DMDAVecGetArrayRead(da, u_star, &arru_star));
    PetscCall(DMDAVecGetArrayRead(da, v_star, &arrv_star));
    PetscCall(DMDAVecGetArrayRead(da, p_prime, &arrp_prime));
    PetscCall(DMStagVecGetArrayRead(stag, UV_star, &arrUV_star));

    PetscCall(DMStagGetLocationSlot(stag, DMSTAG_LEFT, 0, &iU));
    PetscCall(DMStagGetLocationSlot(stag, DMSTAG_DOWN, 0, &iV));

    for (j = info.ys; j < info.ys + info.ym; j++)
        for (i = info.xs; i < info.xs + info.xm; i++) {
            ppp = arrp_prime[j][i];
            /* Left wall. */
            if (i == 0) {
                ppw = arrp_prime[j][i];
                pw = arrp[j][i];
            } else {
                ppw = arrp_prime[j][i-1];
                pw = arrp[j][i-1];
            }
            /* Right wall. */
            if (i == info.mx - 1) {
                ppe = arrp_prime[j][i];
                pe = arrp[j][i];
            } else {
                ppe = arrp_prime[j][i+1];
                pe = arrp[j][i+1];
            }
            /* Bottom wall. */
            if (j == 0) {
                pps = arrp_prime[j][i];
                ps = arrp[j][i];
            } else {
                pps = arrp_prime[j-1][i];
                ps = arrp[j-1][i];
            }
            /* Top wall. */
            if (j == info.my - 1) {
                ppn = arrp_prime[j][i];
                pn = arrp[j][i];
            } else {
                ppn = arrp_prime[j+1][i];
                pn = arrp[j+1][i];
            }

            arru[j][i] = arru_star[j][i] - dt/rho * (ppe - ppw) / (2.0*hx);
            arrv[j][i] = arrv_star[j][i] - dt/rho * (ppn - pps) / (2.0*hy);
            arrp[j][i] += arrp_prime[j][i]
                          - mu*dt/(2.0*rho) * ((ppe - 2.0*ppp + ppw) / (hx*hx) + ((ppn - 2.0*ppp + pps) / (hy*hy)));
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
                    arrUV[j][i][iU] = arrUV_star[j][i][iU] - dt/rho * (arrp_prime[j][i] - arrp_prime[j][i-1]) / hx;
            }
            if (i < info.mx) {
                /* Bottom wall. */
                if (j == 0)
                    arrUV[j][i][iV] = 0.0;
                /* Top wall. */
                else if (j == info.my)
                    arrUV[j][i][iV] = 0.0;
                else
                    arrUV[j][i][iV] = arrUV_star[j][i][iV] - dt/rho * (arrp_prime[j][i] - arrp_prime[j-1][i]) / hy;
            }
        }

    PetscCall(DMDAVecRestoreArray(da, u, &arru));
    PetscCall(DMDAVecRestoreArray(da, v, &arrv));
    PetscCall(DMDAVecRestoreArray(da, p, &arrp));
    PetscCall(DMStagVecRestoreArray(stag, UV, &arrUV));
    PetscCall(DMDAVecRestoreArrayRead(da, u_star, &arru_star));
    PetscCall(DMDAVecRestoreArrayRead(da, v_star, &arrv_star));
    PetscCall(DMDAVecRestoreArrayRead(da, p_prime, &arrp_prime));
    PetscCall(DMStagVecRestoreArrayRead(stag, UV_star, &arrUV_star));

    DMLocalToLocalBegin(da, u, INSERT_VALUES, u);
    DMLocalToLocalEnd(da, u, INSERT_VALUES, u);
    DMLocalToLocalBegin(da, v, INSERT_VALUES, v);
    DMLocalToLocalEnd(da, v, INSERT_VALUES, v);
    DMLocalToLocalBegin(da, p, INSERT_VALUES, p);
    DMLocalToLocalEnd(da, p, INSERT_VALUES, p);

    PetscCall(VecCopy(Nu, Nu_prev));
    PetscCall(VecCopy(Nv, Nv_prev));

    return 0;
}

PetscErrorCode ComputeRHSUstar(KSP ksp, Vec b, void *ctx) {
    PetscReal hx, hy;
    PetscReal **arrb;
    const PetscReal **arrNu, **arrNu_prev, **arru, **arrv, **arrp;
    PetscReal up, uw, ue, us, un, pw, pe;
    PetscInt i, j;

    hx = 1.0 / info.mx;
    hy = 1.0 / info.my;

    PetscCall(DMDAVecGetArray(dau, b, &arrb));
    PetscCall(DMDAVecGetArrayRead(da, Nu, &arrNu));
    PetscCall(DMDAVecGetArrayRead(da, Nu_prev, &arrNu_prev));
    PetscCall(DMDAVecGetArrayRead(da, u, &arru));
    PetscCall(DMDAVecGetArrayRead(da, v, &arrv));
    PetscCall(DMDAVecGetArrayRead(da, p, &arrp));

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

            arrb[j][i] = up - 1.5*dt * arrNu[j][i] + 0.5*dt * arrNu_prev[j][i] - dt/rho * (pe-pw)/(2.0*hx)
                         + mu*dt/(2.0*rho) * ((uw-2.0*up+ue)/(hx*hx) + (us-2.0*up+un)/(hy*hy));
            if (j == info.my - 1)
                arrb[j][i] += mu*dt/(rho*hy*hy);
        }

    PetscCall(DMDAVecRestoreArray(dau, b, &arrb));
    PetscCall(DMDAVecRestoreArrayRead(da, Nu, &arrNu));
    PetscCall(DMDAVecRestoreArrayRead(da, Nu_prev, &arrNu_prev));
    PetscCall(DMDAVecRestoreArrayRead(da, u, &arru));
    PetscCall(DMDAVecRestoreArrayRead(da, v, &arrv));
    PetscCall(DMDAVecRestoreArrayRead(da, p, &arrp));

    return 0;
}

PetscErrorCode ComputeRHSVstar(KSP ksp, Vec b, void *ctx) {
    PetscReal hx, hy;
    PetscReal **arrb;
    const PetscReal **arrNv, **arrNv_prev, **arru, **arrv, **arrp;
    PetscReal vp, vw, ve, vs, vn, ps, pn;
    PetscInt i, j;

    hx = 1.0 / info.mx;
    hy = 1.0 / info.my;

    PetscCall(DMDAVecGetArray(dav, b, &arrb));
    PetscCall(DMDAVecGetArrayRead(da, Nv, &arrNv));
    PetscCall(DMDAVecGetArrayRead(da, Nv_prev, &arrNv_prev));
    PetscCall(DMDAVecGetArrayRead(da, u, &arru));
    PetscCall(DMDAVecGetArrayRead(da, v, &arrv));
    PetscCall(DMDAVecGetArrayRead(da, p, &arrp));

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

            arrb[j][i] = vp - 1.5*dt * arrNv[j][i] + 0.5*dt * arrNv_prev[j][i] - dt/rho * (pn-ps)/(2.0*hy)
                         + mu*dt/(2.0*rho) * ((vw-2.0*vp+ve)/(hx*hx) + (vs-2.0*vp+vn)/(hy*hy));
        }

    PetscCall(DMDAVecRestoreArray(dav, b, &arrb));
    PetscCall(DMDAVecRestoreArrayRead(da, Nv, &arrNv));
    PetscCall(DMDAVecRestoreArrayRead(da, Nv_prev, &arrNv_prev));
    PetscCall(DMDAVecRestoreArrayRead(da, u, &arru));
    PetscCall(DMDAVecRestoreArrayRead(da, v, &arrv));
    PetscCall(DMDAVecRestoreArrayRead(da, p, &arrp));

    return 0;
}

PetscErrorCode ComputeRHSPprime(KSP ksp, Vec b, void *ctx) {
    PetscReal hx, hy;
    PetscReal **arrb;
    const PetscReal ***arrUV_star;
    PetscInt iU, iV;
    MatNullSpace ns;
    PetscInt i, j;

    hx = 1.0 / info.mx;
    hy = 1.0 / info.my;

    PetscCall(DMDAVecGetArray(dap, b, &arrb));
    PetscCall(DMStagVecGetArrayRead(stag, UV_star, &arrUV_star));

    PetscCall(DMStagGetLocationSlot(stag, DMSTAG_LEFT, 0, &iU));
    PetscCall(DMStagGetLocationSlot(stag, DMSTAG_DOWN, 0, &iV));

    for (j = info.ys; j < info.ys + info.ym; j++)
        for (i = info.xs; i < info.xs + info.xm; i++)
            arrb[j][i] = -rho/dt * (hy * (arrUV_star[j][i+1][iU] - arrUV_star[j][i][iU])
                                    + hx * (arrUV_star[j+1][i][iV] - arrUV_star[j][i][iV]));

    PetscCall(DMDAVecRestoreArray(dap, b, &arrb));
    PetscCall(DMStagVecRestoreArrayRead(stag, UV_star, &arrUV_star));

    /* Remove null space. */
    PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &ns));
    PetscCall(MatNullSpaceRemove(ns, b));
    PetscCall(MatNullSpaceDestroy(&ns));

    return 0;
}

PetscErrorCode ComputeOperatorsUVstar(KSP ksp, Mat J, Mat Jpre, void *ctx) {
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
            v[0] = 1.0 + 2.0*mu*dt/rho * (1.0/(hx*hx) + 1.0/(hy*hy));
            v[1] = v[2] = v[3] = v[4] = 0.0;
            ncols = 1;

            /* Not left wall. */
            if (i != 0) {
                col[ncols].i = i - 1;
                col[ncols].j = j;
                v[0] -= mu*dt/(2.0*rho*hx*hx);
                v[ncols] = -mu*dt/(2.0*rho*hx*hx);
                ncols++;
            }
            /* Not right wall. */
            if (i != info.mx - 1) {
                col[ncols].i = i + 1;
                col[ncols].j = j;
                v[0] -= mu*dt/(2.0*rho*hx*hx);
                v[ncols] = -mu*dt/(2.0*rho*hx*hx);
                ncols++;
            }
            /* Not bottom wall. */
            if (j != 0) {
                col[ncols].i = i;
                col[ncols].j = j - 1;
                v[0] -= mu*dt/(2.0*rho*hy*hy);
                v[ncols] = -mu*dt/(2.0*rho*hy*hy);
                ncols++;
            }
            /* Not top wall. */
            if (j != info.my - 1) {
                col[ncols].i = i;
                col[ncols].j = j + 1;
                v[0] -= mu*dt/(2.0*rho*hy*hy);
                v[ncols] = -mu*dt/(2.0*rho*hy*hy);
                ncols++;
            }

            PetscCall(MatSetValuesStencil(Jpre, 1, &row, ncols, col, v, INSERT_VALUES));
        }

    PetscCall(MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY));

    return 0;
}

PetscErrorCode ComputeOperatorsPprime(KSP ksp, Mat J, Mat Jpre, void *ctx) {
    DM da;
    DMDALocalInfo info;
    PetscReal hx, hy;
    MatStencil row, col[5];
    PetscReal v[5];
    PetscInt ncols;
    MatNullSpace ns;
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
    PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &ns));
    PetscCall(MatSetNullSpace(J, ns));
    PetscCall(MatNullSpaceDestroy(&ns));

    return 0;
}
