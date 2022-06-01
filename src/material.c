#include <petscsys.h>
#include "../inc/private/materialimpl.h"
#include "../inc/private/objectimpl.h"

static PetscErrorCode FcObjecDestroy_Material(FcObject *obj);

PetscErrorCode FcMaterialCreate(PetscReal rho, PetscReal mu, FcMaterial *mat) {
    FcMaterial m;

    /* Create the new material. */
    PetscCall(PetscNew(&m));

    /* Initialize. */
    FcObjectInit((FcObject)m, PETSC_COMM_SELF, "FcMaterial");
    m->obj.ops.destroy = FcObjecDestroy_Material;
    m->rho = rho;
    m->mu = mu;

    /* Get the reference. */
    PetscCall(FcObjectGetReference((FcObject)m, (FcObject *)mat));

    return 0;
}

PetscErrorCode FcMaterialDestroy(FcMaterial *mat) {
    PetscCall(FcObjectRestoreReference((FcObject *)mat));

    return 0;
}

static PetscErrorCode FcObjecDestroy_Material(FcObject *obj) {
    FcMaterial mat = (FcMaterial)(*obj);

    /* Free memory. */
    PetscCall(PetscFree(mat));

    return 0;
}
