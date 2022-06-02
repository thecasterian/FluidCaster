/**
 * @file viewer.h
 * @brief Mesh and solution viewer.
 */

#ifndef VIEWER_H
#define VIEWER_H

#include <petscsystypes.h>

/**
 * @brief Viewer interface that helps view mesh and solution.
 */
typedef struct _p_FcViewer *FcViewer;

/**
 * @brief Type of a viewer.
 */
typedef const char *FcViewerType;

/** Parallel CGNS viewer. */
#define FC_VIEWER_PCGNS "PCGNS"

/**
 * @brief Destroys a viewer.
 *
 * @param viewer Viewer to close.
 */
PetscErrorCode FcViewerDestroy(FcViewer *viewer);

#endif
