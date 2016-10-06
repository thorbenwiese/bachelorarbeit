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
  GtEoplistReader *eoplist_reader;
  GtCigarOp co;
  GtUword p, q, count = 0, num_chars_in_v = 0, num_chars_in_u = 0, 
          v_len = 0;

  eoplist_reader = gt_eoplist_reader_new(eoplist);

  /* p is a factor to dynamically adjust the interval borders */
  /* if the sequence starts at pos 0 then p should be 1 */
  p = MAX(1, ceil(tp_list->start1 / tp_list->delta));

  /* number of intervals */
  GtUword tau = ceil(tp_list->end1 / tp_list->delta) - 
                floor(tp_list->start1 / tp_list->delta);
  
  /* Trace Points in useq */
  gt_assert(tau > 1);
  GtUword *u_tp = gt_malloc((tau) * sizeof *u_tp);
  gt_assert(u_tp != NULL);

  for(q = 0; q <= tau - 1; q++)
  {
    u_tp[q] = (p + q) * tp_list->delta - 1;
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
      if(num_chars_in_u == u_tp[count])
      {
        tp_list->TP[v_len++] = num_chars_in_v;

        // do not increment count in the last interval
        if(count == tau - 1)
        {
          tp_list->TP_len = v_len;

          /* print Trace Points -> for debugging only */
          gt_print_tracepoint_list(tp_list);

          break;
        }
        else
        {
          count++;
        }
      }
    }
  }
  gt_eoplist_reader_delete(eoplist_reader);
  gt_free(u_tp);
}


/* function to generate GtEoplist from TracePointList */
GtEoplist *gt_tracepoint_decode(TracePointList *tp_list)
{
  GtEoplist *eoplist = NULL;
  GtEoplist *final_eoplist = NULL;
  GtUchar *usub = NULL, *vsub = NULL;
  GtUword edist = 0, final_edist = 0, ulen, vlen, i;
  FrontEdistTrace *fet;

  fet = front_edist_trace_new();
  eoplist = gt_eoplist_new();
  final_eoplist = gt_eoplist_new();

  for(i = 0; i < tp_list->TP_len; i++)
  {
    gt_eoplist_reset(eoplist);
    if(i == 0)
    {
      usub = (unsigned char *) tp_list->useq;
      vsub = (unsigned char *) tp_list->vseq;
      ulen = tp_list->delta;
      vlen = tp_list->TP[0] + 1;
    }
    else if(i == tp_list->TP_len - 1)
    {
      usub = (unsigned char *) tp_list->useq + (i * tp_list->delta);
      vsub = (unsigned char *) tp_list->vseq + (tp_list->TP[i - 1] + 1);
      ulen = tp_list->end1 - i * tp_list->delta;
      vlen = tp_list->end2 - tp_list->TP[i - 1] + 1;
    }
    else
    {
      usub = (unsigned char *) (tp_list->useq + (i * tp_list->delta));
      vsub = (unsigned char *) (tp_list->vseq + (tp_list->TP[i - 1] + 1));
      ulen = (i+1) * tp_list->delta;
      vlen = tp_list->TP[i] + 1;
    }

    edist = front_edist_trace_eoplist(eoplist,
                                      fet,
                                      usub,
                                      ulen,
                                      vsub,
                                      vlen,
                                      false);
    gt_assert(edist == gt_eoplist_unit_cost(eoplist));

    final_edist += edist;
    gt_eoplist_append(final_eoplist, eoplist);

  }
  gt_assert(final_edist == gt_eoplist_unit_cost(final_eoplist));

  front_edist_trace_delete(fet);

  return final_eoplist;
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
    printf("\n");
  }
}
