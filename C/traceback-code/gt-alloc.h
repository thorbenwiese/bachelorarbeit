#ifndef GT_ALLOC_H
#define GT_ALLOC_H
#include <assert.h>

#define gt_assert(X)    assert(X)
#define gt_malloc(X)    malloc(X)
#define gt_realloc(X,Y) realloc(X,Y)
#define gt_calloc(X,Y)  calloc(X,Y)
#define gt_free(X)      free(X)

#endif
