/**
 * @file mesh.h
 * @brief Structured mesh.
 */

#ifndef MESH_H
#define MESH_H

#include <mpi.h>
#include <petscdmtypes.h>
#include <petscsystypes.h>

/**
 * @brief Mesh.
 */
typedef struct _p_FcMesh *FcMesh;

/**
 * @brief Type of a physical domain boundary in a mesh.
 */
typedef enum {
    /** Not periodic. */
    FC_MESH_BOUNDARY_NONE,
    /** Periodic. */
    FC_MESH_BOUNDARY_PERIODIC,
} FcMeshBoundaryType;

/**
 * @brief Mesh informations.
 */
typedef struct {
    /** Dimension. */
    PetscInt dim;
    /** Boundary type in the x-direction. */
    FcMeshBoundaryType bx;
    /** Boundary type in the y-direction. */
    FcMeshBoundaryType by;
    /** Boundary type in the z-direction. */
    FcMeshBoundaryType bz;
    /** Global number of elements in the x-direction. */
    PetscInt mx;
    /** Global number of elements in the y-direction. */
    PetscInt my;
    /** Global number of elements in the z-direction. */
    PetscInt mz;
    /** Local number of elements in the x-direction. */
    PetscInt xm;
    /** Local number of elements in the y-direction. */
    PetscInt ym;
    /** Local number of elements in the z-direction. */
    PetscInt zm;
    /** X-index of the first element in this process. */
    PetscInt xs;
    /** Y-index of the first element in this process. */
    PetscInt ys;
    /** Z-index of the first element in this process. */
    PetscInt zs;
    /** Number of processes in the x-direction. */
    PetscInt px;
    /** Number of processes in the y-direction. */
    PetscInt py;
    /** Number of processes in the z-direction. */
    PetscInt pz;
    /** Mesh. */
    FcMesh mesh;
} FcMeshInfo;

/**
 * @brief Creates a mesh.
 *
 * @param comm MPI communicator.
 * @param bx Type of boundary in the x-direction.
 * @param by Type of boundary in the y-direction.
 * @param mx Global number of elements in the x-direction.
 * @param my Global number of elements in the y-direction.
 * @param px Number of processes in the x-direction, or PETSC_DECIDE to have calculated.
 * @param py Number of processes in the y-direction, or PETSC_DECIDE to have calculated.
 * @param[in] lx Array containing the number of elements in each process along the x-direction, or NULL. If not NULL,
 * the array must be of length as @p px, which cannot be PETSC_DECIDE, and the sum of the entries must be @p mx.
 * @param[in] ly Array containing the number of elements in each process along the y-direction, or NULL. If not NULL,
 * the array must be of length as @p py, which cannot be PETSC_DECIDE, and the sum of the entries must be @p my.
 * @param[out] mesh Resulting mesh.
 */
PetscErrorCode FcMeshCreate2d(MPI_Comm comm, FcMeshBoundaryType bx, FcMeshBoundaryType by, PetscInt mx, PetscInt my,
                              PetscInt px, PetscInt py, const PetscInt lx[], const PetscInt ly[], FcMesh *mesh);

/**
 * @brief Creates a mesh.
 *
 * @param comm MPI communicator.
 * @param bx Type of boundary in the x-direction.
 * @param by Type of boundary in the y-direction.
 * @param bz Type of boundary in the z-direction.
 * @param mx Global number of elements in the x-direction.
 * @param my Global number of elements in the y-direction.
 * @param mz Global number of elements in the z-direction.
 * @param px Number of processes in the x-direction, or PETSC_DECIDE to have calculated.
 * @param py Number of processes in the y-direction, or PETSC_DECIDE to have calculated.
 * @param pz Number of processes in the z-direction, or PETSC_DECIDE to have calculated.
 * @param[in] lx Array containing the number of elements in each process along the x-direction, or NULL. If not NULL,
 * the array must be of length as @p px, which cannot be PETSC_DECIDE, and the sum of the entries must be @p mx.
 * @param[in] ly Array containing the number of elements in each process along the y-direction, or NULL. If not NULL,
 * the array must be of length as @p py, which cannot be PETSC_DECIDE, and the sum of the entries must be @p my.
 * @param[in] lz Array containing the number of elements in each process along the z-direction, or NULL. If not NULL,
 * the array must be of length as @p pz, which cannot be PETSC_DECIDE, and the sum of the entries must be @p mz.
 * @param[out] mesh Resulting mesh.
 */
PetscErrorCode FcMeshCreate3d(MPI_Comm comm, FcMeshBoundaryType bx, FcMeshBoundaryType by, FcMeshBoundaryType bz,
                              PetscInt mx, PetscInt my, PetscInt mz, PetscInt px, PetscInt py, PetscInt pz,
                              const PetscInt lx[], const PetscInt ly[], const PetscInt lz[], FcMesh *mesh);

/**
 * @brief Sets parameters in a mesh from the options database.
 *
 * @param mesh Mesh to set options for.
 */
PetscErrorCode FcMeshSetFromOptions(FcMesh mesh);

/**
 * @brief Sets up the data structures inside a mesh.
 *
 * @param mesh Mesh to set up.
 */
PetscErrorCode FcMeshSetUp(FcMesh mesh);

/**
 * @brief Destroys a mesh.
 *
 * @param[in] mesh Mesh to destroy.
 */
PetscErrorCode FcMeshDestory(FcMesh *mesh);

/**
 * @brief Gets informations about a mesh.
 *
 * @param mesh Mesh.
 * @param[out] info Struct containing the informations.
 */
PetscErrorCode FcMeshGetInfo(FcMesh mesh, FcMeshInfo *info);

/**
 * @brief Gets the number of elements in the x-, y-, and z-directions that are owned by each process.
 *
 * @param mesh Mesh.
 * @param[out] lx Ownership along the x-direction, or NULL.
 * @param[out] ly Ownership along the y-direction, or NULL.
 * @param[out] lz Ownership along the z-direction, or NULL.
 *
 * @note These correspond to the arguments passed to FcMeshCreate2d() and FcMeshCreate3d().
 * @warning User must not free these arrays nor chage their entries.
 */
PetscErrorCode FcMeshGetOwnershipRanges(FcMesh mesh, const PetscInt *lx[], const PetscInt *ly[], const PetscInt *lz[]);

/**
 * @brief Gets the DMDA and DMStag used by a mesh.
 *
 * @param mesh Mesh.
 * @param da DMDA.
 * @param stag DMStag.
 *
 * @warning User must not modify nor destroy these DM objects.
 */
PetscErrorCode FcMeshGetDM(FcMesh mesh, DM *da, DM *stag);

#endif
