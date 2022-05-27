#ifndef MESH_IMPL_H
#define MESH_IMPL_H

#include <petscdm.h>

struct _p_FcMesh {
    /** Dimension. */
    PetscInt dim;
    /** DMDA. */
    DM da;
    /** DMStag. */
    DM stag;
};

#endif
