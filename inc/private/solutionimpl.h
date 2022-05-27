#ifndef SOLUTION_IMPL_H
#define SOLUTION_IMPL_H

#include <petscvec.h>

typedef struct _p_FcMesh *FcMesh;

struct _p_FcSolution {
    /** Mesh. */
    FcMesh mesh;
    /** X-velocity. */
    Vec u;
    /** Y-velocity. */
    Vec v;
    /** Z-velocity. */
    Vec w;
    /** Pressure. */
    Vec p;
};

#endif
