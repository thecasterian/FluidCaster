#ifndef FSM_IMPL_H
#define FSM_IMPL_H

#include <petscksp.h>
#include <petscsystypes.h>
#include "../ns.h"

PetscErrorCode FSMCalculateConvection2d(FcNS ns);
PetscErrorCode FSMCalculateIntermediateVelocity2d(FcNS ns);
PetscErrorCode FSMCalculatePressureCorrection2d(FcNS ns);
PetscErrorCode FSMUpdate2d(FcNS ns);

PetscErrorCode FSMComputeRHSUstar2d(KSP ksp, Vec b, void *ctx);
PetscErrorCode FSMComputeRHSVstar2d(KSP ksp, Vec b, void *ctx);
PetscErrorCode FSMComputeRHSPprime2d(KSP ksp, Vec b, void *ctx);
PetscErrorCode FSMComputeOperatorsUVstar2d(KSP ksp, Mat J, Mat Jpre, void *ctx);
PetscErrorCode FSMComputeOperatorsPprime2d(KSP ksp, Mat J, Mat Jpre, void *ctx);

PetscErrorCode FSMComputeRHSUstar3d(KSP ksp, Vec b, void *ctx);
PetscErrorCode FSMComputeRHSVstar3d(KSP ksp, Vec b, void *ctx);
PetscErrorCode FSMComputeRHSWstar3d(KSP ksp, Vec b, void *ctx);
PetscErrorCode FSMComputeRHSPprime3d(KSP ksp, Vec b, void *ctx);
PetscErrorCode FSMComputeOperatorsUVWstar3d(KSP ksp, Mat J, Mat Jpre, void *ctx);
PetscErrorCode FSMComputeOperatorsPprime3d(KSP ksp, Mat J, Mat Jpre, void *ctx);

#endif
