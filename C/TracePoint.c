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

  /* Trace Points in v sequence */
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
      gt_assert(count <= tau - 1);
      if(num_chars_in_u == u_tp[count])
      {
        v_tp[v_len++] = num_chars_in_v;

        // TODO do not increment count in the last interval
        if(count == tau - 1)
        {
          return v_tp;
        }
        else
        {
          count++;
        }
      }
    }
  }
  fprintf(stderr, "Encode Failed.\n");
  exit(EXIT_FAILURE);
}


/* function to create a GtEoplist from TracePointData and TracePoint Array */
/*
GtEoplist * decode(const *GtUword TP, TracePointData *tp_data)
{
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
   
  return eoplist;
}
*/

/*
# create new intervals from TracePoints and calculate new alignment
  def decode(self, tp):

    assert self.seq1, "First sequence for decode function is empty."
    assert self.seq2, "Second sequence for decode function is empty."
    assert self.delta > 0, "Delta for decode function is <= 0."
    assert tp, "TracePoint Array for decode function is empty."
    assert self.start_seq1 >= 0, \
      "Starting position for first sequence in decode function is < 0."
    assert self.start_seq2 >= 0, \
      "Starting position for second sequence in decode function is < 0."
    assert self.end_seq1 > 0, \
      "End position for first sequence in decode function is <= 0."
    assert self.end_seq2 > 0, \
      "End position for second sequence in decode function is <= 0."

    # calculate CIGAR of intervals
    cigar = ""
    
    aln = Alignment.Alignment(self.seq1, self.seq2, self.start_seq1,
                              self.end_seq1,self.start_seq2,self.end_seq2)

    for i in range(0,len(tp)):

      if i == 0:

        cigar = aln.calc_cigar(self.seq1[0:self.delta], self.seq2[0:tp[i] + 1])
      
      elif i == len(tp) - 1:
 
        cigar += aln.calc_cigar(self.seq1[i*self.delta:len(self.seq1)],
                                self.seq2[tp[i - 1] + 1:len(self.seq2)])

      else:
        
        cigar += aln.calc_cigar(self.seq1[i * self.delta:(i + 1) * self.delta],
                                self.seq2[tp[i - 1] + 1:tp[i] + 1])

    cigar = Cigar_Pattern.combine_cigar(cigar)

    return cigar
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
