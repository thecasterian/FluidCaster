#ifndef SOLUTION_IMPL_H
#define SOLUTION_IMPL_H

#include <petscvec.h>
#include "../mesh.h"
#include "objectimpl.h"

struct _p_FcSolution {
    /** Object. */
    struct _p_FcObject obj;

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
