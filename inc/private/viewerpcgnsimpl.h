#ifndef VIEWER_PCGNS_IMPL_H
#define VIEWER_PCGNS_IMPL_H

#include <pcgnslib.h>
#include "../viewerpcgns.h"
#include "viewerimpl.h"

typedef struct {
    /** File number. */
    int fn;
    /** Base index. */
    int B;
    /** Zone index. */
    int Z;
    /** Solution index. */
    int S;
} VIEWER_PCGNS;

#endif
