#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "gt-alloc.h"
#include "front-with-trace.h"
#include "TracePoint.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

struct TracePointList 
{
  const GtUchar *useq, *vseq;
  GtUword *TP, TP_len, start1, end1, start2, end2, delta;
}; 

/* function to create a TracePoint Array from eoplist*/
void gt_tracepoint_encode(TracePointList *tp_list, 
                          GtEoplist *eoplist)
{
  GtEoplistReader *eoplist_reader = gt_eoplist_reader_new(eoplist);
  GtCigarOp co;
  GtUword q, count = 0, num_chars_in_v = 0, num_chars_in_u = 0, v_len = 0;

  /* number of intervals */
  GtUword tau = ceil((float)(tp_list->end1 - tp_list->start1) / tp_list->delta);
  
  /* Trace Points in useq */
  gt_assert(tau > 1);
  GtUword *u_tp = gt_malloc((tau) * sizeof *u_tp);
  gt_assert(u_tp != NULL);

  for(q = 0; q <= tau-1; q++)
  {
    u_tp[q] = (q + 1) * tp_list->delta - 1;
  }

  /* Trace Point Array */
  tp_list->TP = gt_malloc((tau) * sizeof tp_list->TP);
  gt_assert(tp_list->TP != NULL);

  while (gt_eoplist_reader_next_cigar(&co, eoplist_reader))
  {
    GtUword i;
    for(i = 0; i < co.iteration; i++)
    {
      if(co.eoptype == GtInsertionOp)
      {
        num_chars_in_v++;
      }
      else if(co.eoptype == GtDeletionOp)
      {
        num_chars_in_u++;
      }
      else
      {
        num_chars_in_v++;
        num_chars_in_u++;
      }
      gt_assert(count < tau);
      if(num_chars_in_u == u_tp[count] + 1 && count < tau - 1)
      {
        tp_list->TP[v_len++] = num_chars_in_v - 1;
        count++;
      }
      if(count == tau - 1)
      {
        tp_list->TP_len = v_len;
        break;
      }
    }
  }

  /* print Trace Points -> for debugging only */
  gt_eoplist_reader_delete(eoplist_reader);
  gt_free(u_tp);
}

/* function to generate GtEoplist from TracePointList */
GtEoplist *gt_tracepoint_decode(const TracePointList *tp_list)
{
  GtEoplist *eoplist = gt_eoplist_new();
  FrontEdistTrace *fet = front_edist_trace_new();
  GtUword i, final_edist = 0;

  for(i = 0; i <= tp_list->TP_len; i++)
  {
    const GtUchar *usub = NULL, *vsub = NULL;
    GtUword ulen = tp_list->delta, vlen = 0;

    if(i == 0)
    {
      usub = tp_list->useq;
      vsub = tp_list->vseq;
      vlen = tp_list->TP[0] + 1;
    }
    else if(i == tp_list->TP_len)
    {
      usub = tp_list->useq + (i * tp_list->delta);
      vsub = tp_list->vseq + tp_list->TP[i - 1] + 1;

      gt_assert(tp_list->end1 + 1 >= tp_list->start1 + i * tp_list->delta);
      ulen = tp_list->end1 - tp_list->start1 - i * tp_list->delta + 1;
      gt_assert(tp_list->end2 >= tp_list->start2 + tp_list->TP[i - 1]);
      vlen = tp_list->end2 - tp_list->start2 - tp_list->TP[i - 1];
    }
    else
    {
      usub = tp_list->useq + (i * tp_list->delta);
      vsub = tp_list->vseq + tp_list->TP[i - 1] + 1;
      vlen = tp_list->TP[i] - tp_list->TP[i - 1];
    }

    final_edist += front_edist_trace_eoplist(eoplist,
                                             fet,
                                             usub,
                                             ulen,
                                             vsub,
                                             vlen,
                                             false);
  }
  gt_eoplist_reader_verify(eoplist,
                           tp_list->useq,
                           tp_list->end1 - tp_list->start1 + 1,
                           tp_list->vseq,
                           tp_list->end2 - tp_list->start2 + 1,
                           final_edist,
                           false);
  gt_assert(final_edist == gt_eoplist_unit_cost(eoplist));
  front_edist_trace_delete(fet);
  return eoplist;
}

/* function to set TracePointData */
void gt_tracepoint_list_set(TracePointList *tp_list, 
                            const GtUchar *useq,
                            const GtUchar *vseq,
                            GtUword start1,
                            GtUword end1,
                            GtUword start2,
                            GtUword end2,
                            GtUword delta)
{
  if(tp_list != NULL)
  {
    tp_list->useq = useq;
    tp_list->vseq = vseq;
    tp_list->start1 = start1;
    tp_list->end1 = end1;
    tp_list->start2 = start2;
    tp_list->end2 = end2;
    tp_list->delta = delta;

    gt_assert(useq != NULL);
    gt_assert(vseq != NULL);
    gt_assert(start1 >= 0 && start1 < end1);
    gt_assert(start2 >= 0 && start2 < end2);
    gt_assert(delta > 0);
  }
}

/* function to reset TracePointData */
void gt_tracepoint_list_reset(TracePointList *tp_list)
{
  if(tp_list != NULL)
  {
    tp_list->useq = NULL;
    tp_list->vseq = NULL;
    tp_list->TP = NULL;
    tp_list->TP_len = 0;
    tp_list->start1 = 0;
    tp_list->end1 = 0;
    tp_list->start2 = 0;
    tp_list->end2 = 0;
    tp_list->delta = 0;
  }
}

/* function to create new TracePointList */
TracePointList *gt_tracepoint_list_new(void)
{
  TracePointList *tp_list = gt_malloc(sizeof *tp_list);

  gt_assert(tp_list != NULL);
  gt_tracepoint_list_reset(tp_list);

  return tp_list;
}

/* function to delete and free TracePointList */
void gt_tracepoint_list_delete(TracePointList *tp_list)
{
  if(tp_list != NULL)
  {
    gt_free(tp_list->TP);
    gt_free(tp_list);
  }
}

/* function to print Trace Points of TracePointList */
void gt_print_tracepoint_list(const TracePointList *tp_list)
{
  if(tp_list != NULL)
  {
    GtUword i;
    printf("Trace Points: ");
    for(i = 0; i < tp_list->TP_len; i++)
    {
      printf("%lu ", tp_list->TP[i]);
    }
    printf("\n\n");
  }
}
