#include "TracePoint.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

struct TracePointAlignment {
  const unsigned char *useq;
  const unsigned char *vseq;
  unsigned long ulen;
  unsigned long vlen;
  unsigned long start1;
  unsigned long end1;
  unsigned long start2;
  unsigned long end2;
  unsigned long delta; /* nur parameter */
  const unsigned char *cigar;
  unsigned long ciglen;
}; 

TracePointAlignment *tracepoint_alignment_new()...
void tracepoint_alignment_delete(TracePointAlignment *tp_aln)...

void encode(const TracePointAlignment *tp_aln)
{
  // p is a factor to dynamically adjust the interval borders
  // if the sequence starts at pos 0 then p should be 1
  int p = MAX(1,ceil(tp_aln->start1/tp_aln->delta));

  // number of intervals
  int tau = ceil(tp_aln->end1 / tp_aln->delta) - floor(
            tp_aln->start1 / tp_aln->delta);
  
  // Trace Points in u sequence
  // tau - 1 because the last interval has no Trace Point
  assert(tau > 1);
  unsigned long *u_tp = malloc((tau - 1) * sizeof *u_tp);

  assert(u_tp != NULL);
  for(int q = 0; q <= tau - 1; q++){
    u_tp[q] = (p + q) * tp_aln->delta - 1;
  }

  int cig_count = 0, count = 0;
  int num_chars_in_v = 0, num_chars_in_u = 0;
  int v_len = 0;

  // Trace Points in v sequence
  // tau - 1 because the last interval has no Trace Point
  int *v_tp = (int*)malloc((tau - 1) * sizeof(int));

  // first and last chars are "s
  for(int i = 0; i < tp_aln->ciglen; i++){
    char c = (char)tp_aln->cigar[i];
    if(isdigit(c)){
      cig_count += c - '0';
    }
    else{
      for(int i = 0; i < cig_count; i++){
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
            for(int i = 0; i <= tau - 1; i++){
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

