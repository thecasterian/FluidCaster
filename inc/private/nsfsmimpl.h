#ifndef NS_FSM_IMPL_H
#define NS_FSM_IMPL_H

#include <petscdm.h>
#include <petscksp.h>
#include <petscvec.h>
#include "nsimpl.h"

typedef struct {
    /** DMDA for solving the intermediate x-velocity. */
    DM dau;
    /** DMDA for solving the intermediate y-velocity. */
    DM dav;
    /** DMDA for solving the intermediate z-velocity. */
    DM daw;
    /** DMDA for solving the pressure correction. */
    DM dap;

    /** KSP for solving the intermediate x-velocity. */
    KSP kspu;
    /** KSP for solving the intermediate y-velocity. */
    KSP kspv;
    /** KSP for solving the intermediate z-velocity. */
    KSP kspw;
    /** KSP for solving the pressure correction. */
    KSP kspp;

    /* Pressure at the half time step. */
    Vec p;
    /** Face-centered velocity. */
    Vec UVW;
    /** Intermediate x-velocity. */
    Vec u_star;
    /** Intermediate y-velocity. */
    Vec v_star;
    /** Intermediate z-velocity. */
    Vec w_star;
    /** Face-centered intermediate velocity. */
    Vec UVW_star;
    /** Pressure correction. */
    Vec p_prime;
    /** X-convection term. */
    Vec Nu;
    /** Y-convection term. */
    Vec Nv;
    /** Z-convection term. */
    Vec Nw;
    /** Pressure at the previous half time step. */
    Vec p_prev;
    /** X-convection term at the previous time step. */
    Vec Nu_prev;
    /** Y-convection term at the previous time step. */
    Vec Nv_prev;
    /** Z-convection term at the previous time step. */
    Vec Nw_prev;

    /** Temporaty values used in velocity interpolation. @{ */
    Vec u_tilde, v_tilde, w_tilde;
    /** @} */
} NS_FSM;

#endif
