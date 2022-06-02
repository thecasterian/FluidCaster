#ifndef FSM_IMPL_H
#define FSM_IMPL_H

#include <petscksp.h>
#include <petscsystypes.h>
#include "../ns.h"

PetscErrorCode CalculateConvection(FcNS ns);
PetscErrorCode CalculateIntermediateVelocity(FcNS ns);
PetscErrorCode CalculatePressureCorrection(FcNS ns);
PetscErrorCode Update(FcNS ns);

PetscErrorCode ComputeRHSUstar2d(KSP ksp, Vec b, void *ctx);
PetscErrorCode ComputeRHSVstar2d(KSP ksp, Vec b, void *ctx);
PetscErrorCode ComputeRHSPprime2d(KSP ksp, Vec b, void *ctx);
PetscErrorCode ComputeOperatorsUVstar2d(KSP ksp, Mat J, Mat Jpre, void *ctx);
PetscErrorCode ComputeOperatorsPprime2d(KSP ksp, Mat J, Mat Jpre, void *ctx);

PetscErrorCode ComputeRHSUstar3d(KSP ksp, Vec b, void *ctx);
PetscErrorCode ComputeRHSVstar3d(KSP ksp, Vec b, void *ctx);
PetscErrorCode ComputeRHSWstar3d(KSP ksp, Vec b, void *ctx);
PetscErrorCode ComputeRHSPprime3d(KSP ksp, Vec b, void *ctx);
PetscErrorCode ComputeOperatorsUVWstar3d(KSP ksp, Mat J, Mat Jpre, void *ctx);
PetscErrorCode ComputeOperatorsPprime3d(KSP ksp, Mat J, Mat Jpre, void *ctx);

#endif
