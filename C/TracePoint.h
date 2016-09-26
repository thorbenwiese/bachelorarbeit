#ifndef TRACEPOINT_H
#define TRACEPOINT_H

typedef struct TracePointData TracePointData;

/* The constructor method */
TracePointData *tracepoint_data_new(void);

/* the destructor */
void gt_tracepoint_data_delete(TracePointData *tp_data);

/* reset the data to empty it  */
void gt_tracepoint_data_reset(TracePointData *tp_data);

/* encode function to generate TracePoint Array */
void encode(const TracePointData *tp_data);
#endif
