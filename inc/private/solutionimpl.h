#ifndef SOLUTION_IMPL_H
#define SOLUTION_IMPL_H

#include <petscvec.h>
#include "../mesh.h"
#include "../solution.h"

struct _p_FcSolution {
    /* Mesh. */
    FcMesh mesh;
    /* X-velocity. */
    Vec u;
    /* Y-velocity. */
    Vec v;
    /* Z-velocity. */
    Vec w;
    /* True pressure. */
    Vec p_true;

    /* Pressure at the half time step. */
    Vec p;
    /* Face-centered velocity. */
    Vec UVW;
    /* Intermediate x-velocity. */
    Vec u_star;
    /* Intermediate y-velocity. */
    Vec v_star;
    /* Intermediate z-velocity. */
    Vec w_star;
    /* Face-centered intermediate velocity. */
    Vec UVW_star;
    /* Pressure correction. */
    Vec p_prime;
    /* X-convection term. */
    Vec Nu;
    /* Y-convection term. */
    Vec Nv;
    /* Z-convection term. */
    Vec Nw;
    /* Pressure at the previous half time step. */
    Vec p_prev;
    /* X-convection term at the previous time step. */
    Vec Nu_prev;
    /* Y-convection term at the previous time step. */
    Vec Nv_prev;
    /* Z-convection term at the previous time step. */
    Vec Nw_prev;

    /* Temporaty values used in velocity interpolation. */
    Vec u_tilde, v_tilde, w_tilde;
};

#endif
