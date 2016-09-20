#ifndef TRACEPOINT_H
#define TRACEPOINT_H
typedef struct TracePointAlignment TracePointAlignment;

typedef struct TracePointList TracePointList;

void encode(const TracePointAlignment *tp_aln);
#endif
