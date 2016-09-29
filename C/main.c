#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//TODO nur fuer die Berechnung von tau fuer die Ausgabe der TPs
#include <math.h>

#include "TracePoint.h"
#include "gt-alloc.h"
#include "substring.h"

int main(int argc, char *argv[])
{

  GtUchar *long_useq, *long_vseq, *useq, *vseq;
  GtUword ulen, vlen, start1, end1, start2, end2, delta;
  GtUword *TP;
  TracePointData *tp_data;

  if (argc != 8)
  {
    fprintf(stderr,"Usage: %s <sequence1> <sequence2> <start1> <end1> "
                             "<start2> <end2> <delta>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  /* read and check parameters */
  long_useq = (unsigned char *) strdup(argv[1]);
  gt_assert(long_useq != NULL);

  long_vseq = (unsigned char *) strdup(argv[2]);
  gt_assert(long_vseq != NULL);

  ulen = strlen(argv[1]);
  assert(ulen > 0);

  vlen = strlen(argv[2]);
  assert(vlen > 0);

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

  /* get substring with regard to start and end of sequence */
  useq = substring(long_useq, start1 + 1, end1 - start1 + 1);
  vseq = substring(long_vseq, start2 + 1, end2 -start2 + 1);

  /* create TracePointData */
  tp_data = tracepoint_data_new();
  gt_tracepoint_data_set(tp_data, useq, vseq, ulen, vlen, 
                         start1, end1, start2, end2, delta);

  /* number of intervals */
  GtUword tau = ceil(end1 / delta) - floor(start1 / delta);

  /* encode data to Trace Point Array */
  GtUword i;
  TP = encode(tp_data);
  printf("Trace Points: ");
  for(i = 0; i < tau; i++)
  {
    printf("%lu ",TP[i]);
  }
  printf("\n");

  /* decode TracePoint Array and TracePointData to GtEoplist */
  /*
  GtEoplist *eoplist;
  eoplist = decode(TP, tau - 1, tp_data);
  */

  gt_tracepoint_data_delete(tp_data);
  gt_free(useq);
  gt_free(vseq);

  return 0;
}
