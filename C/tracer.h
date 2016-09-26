#ifndef TRACER_H
#define TRACER_H
#include <stdint.h>
#include "gt-defs.h"
#include "eoplist.h"

typedef struct GtFrontTrace GtFrontTrace;

GtFrontTrace *front_trace_new(void);

void front_trace_delete(GtFrontTrace *front_trace);

void front_trace_reset(GtFrontTrace *front_trace);

void front_trace_add_trace(GtFrontTrace *front_trace,
                           uint8_t backreference,
                           uint32_t localmatch_count);

#define FT_EOP_MISMATCH  1
#define FT_EOP_INSERTION (1 << 1)
#define FT_EOP_DELETION  (1 << 2)

void front_trace2eoplist_directed(GtEoplist *eoplist,
                                  const GtFrontTrace *front_trace,
                                  GtUword distance,
                                  const GtUchar *useq,
                                  GtUword ulen,
                                  const GtUchar *vseq,
                                  GtUword vlen);

#endif
