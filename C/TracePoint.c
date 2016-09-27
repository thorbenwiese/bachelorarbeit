#include "TracePoint.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

struct TracePointData {
  const GtUchar *useq, *vseq;
  GtUword ulen, vlen, start1, end1, start2, end2;
}; 

TracePointData *tracepoint_data_new()...
void tracepoint_data_delete(TracePointData *tp_data)...

GtUword * encode(const TracePointData *tp_data)
{
  // p is a factor to dynamically adjust the interval borders
  // if the sequence starts at pos 0 then p should be 1
  GtUword p = MAX(1,ceil(tp_data->start1/tp_data->delta));

  // number of intervals
  GtUword tau = ceil(tp_data->end1 / tp_data->delta) - floor(
            tp_data->start1 / tp_data->delta);
  
  // Trace Points in u sequence
  // tau - 1 because the last interval has no Trace Point
  assert(tau > 1);
  GtUword *u_tp = malloc((tau - 1) * sizeof *u_tp);

  assert(u_tp != NULL);
  for(int q = 0; q <= tau - 1; q++){
    u_tp[q] = (p + q) * tp_data->delta - 1;
  }

  GtUword cig_count = 0, count = 0, num_chars_in_v = 0, num_chars_in_u = 0, v_len = 0;

  // Trace Points in v sequence
  // tau - 1 because the last interval has no Trace Point
  int *v_tp = (int*)malloc((tau - 1) * sizeof(int));

  // first and last chars are "s
  GtUword i = 0;
  for(i = 0; i < tp_data->ciglen; i++){
    GtUchar c = (char)tp_data->cigar[i];
    if(isdigit(c)){
      cig_count += c - '0';
    }
    else{
      GtUword i = 0;
      for(i = 0; i < cig_count; i++){
        if(c == 'I'){
          num_chars_in_v++;
        }
        else if(c == 'D'){
          num_chars_in_u++;
        }
        else{
          num_chars_in_u++;
          num_chars_in_v++;
        }
        assert(count < (tau - 1));
        if(num_chars_in_u == u_tp[count]){
          v_tp[v_len++] = num_chars_in_v;

          // do not increment count in the last interval
          if(count == tau - 1){
            return v_tp;
            GtUWord i = 0;
            for(i = 0; i <= tau - 1; i++){
              printf("Trace Point %d: %d\n", i + 1, v_tp[i]);
            }
            break;
          }
          else{
            count++;
          }
        }
      }
      cig_count = 0;
    }
  } 
}

/*
void decode(struct TracePointAlignment tp_aln){

}
*/
void gt_tracepoint_data_reset(TracePointData *tp_data)
{
  if(tp_data != NULL)
  {
    tp_data->useq = NULL;
    tp_data->vseq = NULL;
    tp_data->start1 = 0;
    printf("START.11 %s\n",tp_data->start1);
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
