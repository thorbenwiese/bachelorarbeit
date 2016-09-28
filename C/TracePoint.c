#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "gt-alloc.h"
#include "front-with-trace.h"
#include "TracePoint.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

struct TracePointData 
{
  const GtUchar *useq, *vseq;
  GtUword ulen, vlen, start1, end1, start2, end2, delta;
}; 


void encode(const TracePointData *tp_data)
{
  FrontEdistTrace *fet;
  GtEoplist *eoplist;
  GtEoplistReader *eoplist_reader;
  GtCigarOp co;

  GtUword q = 0, edist, count = 0, num_chars_in_v = 0, num_chars_in_u = 0, 
          v_len = 0;

  // p is a factor to dynamically adjust the interval borders
  // if the sequence starts at pos 0 then p should be 1
  GtUword p = MAX(1,ceil(tp_data->start1/tp_data->delta));

  // number of intervals
  GtUword tau = ceil(tp_data->end1 / tp_data->delta) - floor(
            tp_data->start1 / tp_data->delta);
  
  // Trace Points in u sequence
  // tau - 1 because the last interval has no Trace Point
  gt_assert(tau > 1);
  GtUword *u_tp = gt_malloc((tau - 1) * sizeof *u_tp);

  gt_assert(u_tp != NULL);
  for(q = 0; q <= tau - 1; q++)
  {
    u_tp[q] = (p + q) * tp_data->delta - 1;
  }

  // Trace Points in v sequence
  // tau - 1 because the last interval has no Trace Point
  GtUword *v_tp = gt_malloc((tau - 1) * sizeof *v_tp);

  fet = front_edist_trace_new();
  eoplist = gt_eoplist_new();
  edist = front_edist_trace_eoplist(eoplist,
                                    fet,
                                    tp_data->useq,
                                    tp_data->ulen,
                                    tp_data->vseq,
                                    tp_data->vlen,
                                    false);
  gt_assert(edist == gt_eoplist_unit_cost(eoplist));
  eoplist_reader = gt_eoplist_reader_new(eoplist);

  while (gt_eoplist_reader_next_cigar(&co, eoplist_reader))
  {
    //printf("%lu%c",co.iteration, gt_eoplist_pretty_print(co.eoptype, false));
    GtUword i;
    for(i = 0; i < co.iteration; i++)
    {
      if(gt_eoplist_pretty_print(co.eoptype, false) == 'I')
      {
        num_chars_in_v++;
      }
      else if(gt_eoplist_pretty_print(co.eoptype, false) == 'D')
      {
        num_chars_in_u++;
      }
      else
      {
        num_chars_in_v++;
        num_chars_in_u++;
      }
      //printf("uchars: %lu, vchars: %lu\n",num_chars_in_u, num_chars_in_v);
      //gt_assert(count < (tau - 1));
      if(num_chars_in_u == u_tp[count])
      {
        v_tp[v_len++] = num_chars_in_v;

        // TODO do not increment count in the last interval
        if(count == tau - 1)
        {
          //front_edist_trace_delete(fet);
          //gt_eoplist_reader_delete(eoplist_reader);
          //gt_eoplist_delete(eoplist);
          GtUword i = 0;
          for(i = 0; i <= tau - 1; i++)
          {
            printf("Trace Point %lu: %lu\n", i + 1, v_tp[i]);
          }
          //return v_tp;
        }
        else
        {
          count++;
        }
      }
    }
  }
  //fprintf(stderr, "Encode Failed.");
  //exit(EXIT_FAILURE);
}

void gt_tracepoint_data_set(TracePointData *tp_data, 
                            const GtUchar *useq,
                            const GtUchar *vseq,
                            GtUword ulen,
                            GtUword vlen,
                            GtUword start1,
                            GtUword end1,
                            GtUword start2,
                            GtUword end2,
                            GtUword delta)
{
  if(tp_data != NULL)
  {
    tp_data->useq = useq;
    tp_data->vseq = vseq;
    tp_data->ulen = ulen;
    tp_data->vlen = vlen;
    tp_data->start1 = start1;
    tp_data->end1 = end1;
    tp_data->start2 = start2;
    tp_data->end2 = end2;
    tp_data->delta = delta;
    printf("useq: %s\n",useq);
    printf("vseq: %s\n",vseq);
    printf("ulen: %lu\n",tp_data->ulen);
    printf("vlen: %lu\n",tp_data->vlen);
    printf("start1: %lu\n",tp_data->start1);
    printf("end1: %lu\n",tp_data->end1);
    printf("start2: %lu\n",tp_data->start2);
    printf("end2: %lu\n",tp_data->end2);
    printf("delta: %lu\n",tp_data->delta);
  }
}

void gt_tracepoint_data_reset(TracePointData *tp_data)
{
  if(tp_data != NULL)
  {
    tp_data->useq = NULL;
    tp_data->vseq = NULL;
    tp_data->ulen = 0;
    tp_data->vlen = 0;
    tp_data->start1 = 0;
    tp_data->end1 = 0;
    tp_data->start2 = 0;
    tp_data->end2 = 0;
    tp_data->delta = 0;
  }
}

TracePointData *tracepoint_data_new(void)
{
  TracePointData *tp_data = gt_malloc(sizeof *tp_data);

  gt_assert(tp_data != NULL);
  gt_tracepoint_data_reset(tp_data);

  return tp_data;
}

void gt_tracepoint_data_delete(TracePointData *tp_data)
{
  if(tp_data != NULL)
  {
    gt_free(tp_data);
  }
}
