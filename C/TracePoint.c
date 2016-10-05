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
  GtUword *TP;
  GtUword TP_len, start1, end1, start2, end2, delta;
}; 

/* function to create a TracePoint Array from eoplist*/
void gt_tracepoint_encode(TracePointList *tp_list, GtEoplist *eoplist)
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

  /* Trace Points in vseq */
  GtUword *v_tp = gt_malloc((tau) * sizeof *v_tp);
  gt_assert(v_tp != NULL);

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
        v_tp[v_len++] = num_chars_in_v;

        // do not increment count in the last interval
        if(count == tau - 1)
        {
          tp_list->TP_len = v_len;
          tp_list->TP = v_tp;

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
  gt_free(v_tp);
}


/* function to generate GtEoplist from TracePointList */
GtEoplist *gt_tracepoint_decode(TracePointList *tp_list)
{
  GtEoplist *eoplist = NULL;
  GtEoplistReader *eoplist_reader = NULL;
  GtUchar *usub = NULL, *vsub = NULL;
  GtUword edist, ulen, vlen, total_edist = 0, i;
  GtCigarOp co;
  FrontEdistTrace *fet;


  gt_print_tracepoint_list(tp_list);
  printf("BEFORE FET\n");
  //printf("%s\n%s\n%lu %lu %lu %lu %lu %lu %lu\n", tp_list->useq, tp_list->vseq, tp_list->start1, tp_list->end1, tp_list->start2,
  //                                          tp_list->end2, tp_list->delta, tp_list->TP_len, tp_list->TP[0]);
  //TP falsch
  fet = front_edist_trace_new();
  eoplist = gt_eoplist_new();

  printf("TRACE POINT ENCODE: ");
  //for(i = 0; i < tp_list->TP_len; i++)
  //{
  //  printf("%lu ",tp_list->TP[i]);
  //}
  printf("\n");
  for(i = 0; i < tp_list->TP_len; i++)
  {
    printf("I: %lu\n", i);
    //printf("TP %lu\n", tp_list->TP[i]);
    if(i == 0)
    {
      usub = (unsigned char *) tp_list->useq;
      vsub = (unsigned char *) tp_list->vseq;
      ulen = tp_list->delta;
      vlen = tp_list->TP[0] + 1;
      printf("i=0\n%s\n%s\n", usub, vsub);
    }
    else if(i == tp_list->TP_len - 1)
    {
      usub = (unsigned char *) tp_list->useq + (i * tp_list->delta);
      vsub = (unsigned char *) tp_list->vseq + (tp_list->TP[i - 1] + 1);
      ulen = tp_list->end1 - i * tp_list->delta;
      vlen = tp_list->end2 - tp_list->TP[i - 1] + 1;
      printf("i=TPlen-1 %s\n%s\n", usub, vsub);
    }
    else
    {
      //printf("%s\n",tp_list->useq + (i * tp_list->delta));
      //printf("%lu\n", tp_list->TP[0]);
      //printf("%s\n",tp_list->vseq + (tp_list->TP[i - 1] + 1));
      usub = (unsigned char *) (tp_list->useq + (i * tp_list->delta));
      vsub = (unsigned char *) (tp_list->vseq + (tp_list->TP[i - 1] + 1));
      ulen = tp_list->end1 - (i+1) * tp_list->delta;
      vlen = tp_list->end2 - tp_list->TP[i] + 1;
      printf("i=... %s\n%s\n", usub, vsub);
    }
    printf("AFTER FET\n");
    //printf("%s\n%s\n%lu %lu %lu %lu %lu %lu %lu\n", tp_list->useq, tp_list->vseq, tp_list->start1, tp_list->end1, 
    //                                                tp_list->start2,tp_list->end2, tp_list->delta, tp_list->TP_len, 
    //                                                tp_list->TP[0]);
    gt_print_tracepoint_list(tp_list);

    edist = front_edist_trace_eoplist(eoplist,
                                      fet,
                                      usub,
                                      ulen,
                                      vsub,
                                      vlen,
                                      false);

    total_edist += edist;
  } 
  gt_assert(total_edist == gt_eoplist_unit_cost(eoplist));
  eoplist_reader = gt_eoplist_reader_new(eoplist);
  // ... read cigar and concatenate?
  while (gt_eoplist_reader_next_cigar(&co, eoplist_reader))
  {
    printf("%lu%c",co.iteration, gt_eoplist_pretty_print(co.eoptype, false));
  }
  printf("\n");

  gt_eoplist_reader_delete(eoplist_reader);
  //gt_free(usub);
  //gt_free(vsub);

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
    gt_free(tp_list);
  }
}

void gt_print_tracepoint_list(const TracePointList *tp_list)
{
  if(tp_list != NULL)
  {
    GtUword i;
    printf("%s\n%s\n", tp_list->useq, tp_list->vseq);
    printf("TP_len: %lu\n",tp_list->TP_len);
    printf("Trace Points: ");
    for(i = 0; i < tp_list->TP_len; i++)
    {
      printf("%lu ", tp_list->TP[i]);
    }
    printf("\n");
  }
}
