#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "TracePoint.h"

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
  
  TracePointData *tp_data;
  tp_data = tracepoint_data_new();

  useq = (unsigned char *) strdup(argv[1]);
  assert(useq != NULL);

  ulen = strlen(argv[1]);
  assert(ulen > 0);

  vseq = (unsigned char *) strdup(argv[2]);
  assert(vseq != NULL);

  vlen = strlen(argv[2]);
  assert(vlen > 0);

  start1 = atoi(argv[3]);

  end1 = atoi(argv[4]);
  assert(end1 > 0);
  assert(start1 < end1);

  start2 = atoi(argv[5]);

  end2 = atoi(argv[6]);
  assert(end2 > 0);
  assert(start2 < end2);

  delta = atoi(argv[7]);
  assert(delta > 0);

  //tp_data->start1 = start1;  
  //printf("START.11 %s\n",tp_data->start1);
  /*
  tp_data->useq = useq;
  tp_data->ulen = ulen;
  tp_data->vseq = vseq;
  tp_data->vlen = vlen;
  tp_data->start1 = start1;
  tp_data->end1 = end1;
  tp_data->start2 = start2;
  tp_data->end2 = end2;
  tp_data->delta = delta;
  */

  //GtUword TP;
  //TP = encode(tp_data);
  printf("encode\n");
  //encode(tp_data);

  gt_tracepoint_data_delete(tp_data);

  return 0;
}
