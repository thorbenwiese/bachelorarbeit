#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "gt-alloc.h"
#include "front-with-trace.h"
#include "TracePoint.h"
#include "substring.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

struct TracePointData 
{
  const GtUchar *useq, *vseq;
  GtUword ulen, vlen, start1, end1, start2, end2, delta;
}; 

/* function to create a TracePoint Array from TracePointData*/
GtUword * encode(const TracePointData *tp_data)
{
  FrontEdistTrace *fet;
  GtEoplist *eoplist;
  GtEoplistReader *eoplist_reader;
  GtCigarOp co;

  GtUword q = 0, edist, count = 0, num_chars_in_v = 0, num_chars_in_u = 0, 
          v_len = 0;

  /* p is a factor to dynamically adjust the interval borders */
  /* if the sequence starts at pos 0 then p should be 1 */
  GtUword p = MAX(1,ceil(tp_data->start1/tp_data->delta));

  /* number of intervals */
  GtUword tau = ceil(tp_data->end1 / tp_data->delta) - floor(
            tp_data->start1 / tp_data->delta);
  
  /* Trace Points in useq */
  /* tau - 1 because the last interval has no Trace Point */
  gt_assert(tau > 1);
  GtUword *u_tp = gt_malloc((tau - 1) * sizeof *u_tp);

  gt_assert(u_tp != NULL);
  for(q = 0; q <= tau - 1; q++)
  {
    u_tp[q] = (p + q) * tp_data->delta - 1;
  }

  /* Trace Points in vseq */
  /* tau - 1 because the last interval has no Trace Point */
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
      gt_assert(count <= tau - 1);
      if(num_chars_in_u == u_tp[count])
      {
        v_tp[v_len++] = num_chars_in_v;

        // do not increment count in the last interval
        if(count == tau - 1)
        {
          gt_free(u_tp);
          gt_eoplist_reader_delete(eoplist_reader);
          gt_eoplist_delete(eoplist);
          return v_tp;
        }
        else
        {
          count++;
        }
      }
    }
  }
  gt_free(u_tp);
  gt_eoplist_reader_delete(eoplist_reader);
  gt_eoplist_delete(eoplist);
  fprintf(stderr, "Encode Failed.\n");
  exit(EXIT_FAILURE);
}


/* function to create a GtEoplist from TracePointData and TracePoint Array */
/*
GtEoplist * decode(const GtUword *TP, GtUword TPlen, TracePointData *tp_data)
{
  FrontEdistTrace *fet;
  GtEoplist *eoplist;

  fet = front_edist_trace_new();
  eoplist = gt_eoplist_new();

  GtUchar *usub, *vsub;
  GtUword i, usublen, vsublen;
  for(i = 1; i <= TPlen; i++)
  {
    if(i == 0)
    {
      usub = substring(tp_data->useq, 0, tp_data->delta);
      vsub = substring(tp_data->vseq, 0, TP[0] + 1);
    }
    else if(i == TPlen - 1)
    {
      usub = substring(tp_data->useq, i * tp_data->delta, tp_data->ulen);
      vsub = substring(tp_data->vseq, TP[i - 1] + 1, tp_data->vlen);
    }
    else
    {
      usub = substring(tp_data->useq, i * tp_data->delta, (i + 1) * tp_data->delta);
      vsub = substring(tp_data->vseq, TP[i - 1] + 1, TP[i] + 1);
    }
    usublen = strlen(usub);
    vsublen = strlen(vsub);
    edist = front_edist_trace_eoplist(eoplist,
                                      fet,
                                      tp_data->usub,
                                      tp_data->usublen,
                                      tp_data->vsub,
                                      tp_data->vsublen,
                                      false);
    gt_assert(edist == gt_eoplist_unit_cost(eoplist));
    eoplist_reader = gt_eoplist_reader_new(eoplist);
    // ... read cigar and concatenate?
    while (gt_eoplist_reader_next_cigar(&co, eoplist_reader))
    {
      printf("%lu%c",co.iteration, gt_eoplist_pretty_print(co.eoptype, false));
    }
   
  return eoplist;
}
*/

/* function to set TracePointData */
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
  }
}

/* function to reset TracePointData */
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

/* function to create new TracePointData */
TracePointData *tracepoint_data_new(void)
{
  TracePointData *tp_data = gt_malloc(sizeof *tp_data);

  gt_assert(tp_data != NULL);
  gt_tracepoint_data_reset(tp_data);

  return tp_data;
}

/* function to delete and free TracePointData */
void gt_tracepoint_data_delete(TracePointData *tp_data)
{
  if(tp_data != NULL)
  {
    gt_free(tp_data);
  }
}
