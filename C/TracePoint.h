#ifndef TRACEPOINT_H
#define TRACEPOINT_H

#include "gt-defs.h"

typedef struct TracePointData TracePointData;

/* The constructor method */
TracePointData *tracepoint_data_new(void);

/* the destructor */
void gt_tracepoint_data_delete(TracePointData *tp_data);

/* reset the data to empty it  */
void gt_tracepoint_data_reset(TracePointData *tp_data);

/* set the data */
void gt_tracepoint_data_set(TracePointData *tp_data, 
                            const GtUchar *useq,
                            const GtUchar *vseq,
                            GtUword ulen,
                            GtUword vlen,
                            GtUword start1,
                            GtUword end1,
                            GtUword start2,
                            GtUword end2,
                            GtUword delta);

/* encode function to generate TracePoint Array */
void encode(const TracePointData *tp_data);
#endif
