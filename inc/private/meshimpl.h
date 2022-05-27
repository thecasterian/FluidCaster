#ifndef MESH_IMPL_H
#define MESH_IMPL_H

#include <petscdm.h>
#include "objectimpl.h"

struct _p_FcMesh {
    /** Object. */
    struct _p_FcObject obj;

    /** Dimension. */
    PetscInt dim;
    /** DMDA. */
    DM da;
    /** DMStag. */
    DM stag;
};

#endif
