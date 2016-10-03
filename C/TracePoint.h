#ifndef TRACEPOINT_H
#define TRACEPOINT_H

#include "gt-defs.h"

typedef struct TracePointList TracePointList;

/* The constructor method */
TracePointList *gt_tracepoint_list_new(void);

/* the destructor */
void gt_tracepoint_list_delete(TracePointList *tp_list);

/* reset the data to empty it  */
void gt_tracepoint_list_reset(TracePointList *tp_list);

/* set the data */
void gt_tracepoint_list_set(TracePointList *tp_list, 
                            const GtUchar *useq,
                            const GtUchar *vseq,
                            GtUword *TP,
                            GtUword TP_len,
                            GtUword start1,
                            GtUword end1,
                            GtUword start2,
                            GtUword end2,
                            GtUword delta);

/* print Trace Point List */
void gt_print_tracepoint_list(const TracePointList *tp_list);

/* encode function to generate TracePoint Array */
void gt_tracepoint_encode(TracePointList *tp_list, GtEoplist *eoplist);

#endif
