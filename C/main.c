#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "TracePoint.h"
#include "gt-alloc.h"

int main(int argc, char *argv[])
{

  unsigned char *useq, *vseq;
  GtUword ulen, vlen, start1, end1, start2, end2, delta;

  if (argc != 8)
  {
    fprintf(stderr,"Usage: %s <sequence1> <sequence2> <start1> <end1> "
                             "<start2> <end2> <delta>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  /* read and check parameters */
  useq = (unsigned char *) strdup(argv[1]);
  gt_assert(useq != NULL);

  ulen = strlen(argv[1]);
  gt_assert(ulen > 0);

  vseq = (unsigned char *) strdup(argv[2]);
  gt_assert(vseq != NULL);

  vlen = strlen(argv[2]);
  gt_assert(vlen > 0);

  start1 = atoi(argv[3]);

  end1 = atoi(argv[4]);
  gt_assert(end1 > 0);
  gt_assert(start1 < end1);

  start2 = atoi(argv[5]);

  end2 = atoi(argv[6]);
  gt_assert(end2 > 0);
  gt_assert(start2 < end2);

  delta = atoi(argv[7]);
  gt_assert(delta > 0);

  /* create TracePointData */
  TracePointData *tp_data;
  tp_data = tracepoint_data_new();
  gt_tracepoint_data_set(tp_data, useq, vseq, ulen, vlen, 
                         start1, end1, start2, end2, delta);

  /* encode data to Trace Point Array */
  GtUword *TP;
  TP = encode(tp_data);

  gt_tracepoint_data_delete(tp_data);
  gt_free(useq);
  gt_free(vseq);

  return 0;
}