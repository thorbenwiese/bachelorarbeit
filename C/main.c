#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
// ctype.h for isdigit()
#include <ctype.h>

#include "traceback-code/gt-defs.h"
#include "traceback-code/gt-alloc.h"
#include "traceback-code/front-with-trace.h"

#include "TracePoint.h"

int main(int argc, char *argv[]){

  unsigned char *useq, *vseq;//, *cig;
  GtUword ulen, vlen, edist, start1, end1, start2, end2, delta;//, ciglen;
  FrontEdistTrace *fet;
  GtEoplist *eoplist;
  GtEoplistReader *eoplist_reader;
  GtCigarOp co;

  if (argc != 8){
    fprintf(stderr,"Usage: %s <sequence1> <sequence2> <start1> <end1> <start2> <end2> <delta>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  TracePointData *tp_data;
  tp_data = tracepoint_data_new();

  useq = (unsigned char *) strdup(argv[1]);
  ulen = strlen(argv[1]);
  vseq = (unsigned char *) strdup(argv[2]);
  vlen = strlen(argv[2]);
  start1 = atoi(argv[3]);
  end1 = atoi(argv[4]);
  start2 = atoi(argv[5]);
  end2 = atoi(argv[6]);
  delta = atoi(argv[7]);

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
  GtUword *TP = encode(tp_data);
  printf("TP\n");
  printf("Trace Point 1: %lu\n",  TP[0]);

  //cig = (unsigned char *) strdup(argv[8]);
  //ciglen = strlen(argv[8]);

  // create CIGAR String
  fet = front_edist_trace_new();
  eoplist = gt_eoplist_new();
  edist = front_edist_trace_eoplist(eoplist,
                                    fet,
                                    useq,
                                    ulen,
                                    vseq,
                                    vlen,
                                    false);
  assert(edist == gt_eoplist_unit_cost(eoplist));
  eoplist_reader = gt_eoplist_reader_new(eoplist);

  while (gt_eoplist_reader_next_cigar(&co, eoplist_reader))
  {
    printf("%lu%c",co.iteration, gt_eoplist_pretty_print(co.eoptype, false));
  }
  printf("\n");
  /* 
  // struct to store input
  TracePointAlignment tp_aln;

  tp_aln.useq = useq;
  tp_aln.ulen = ulen;
  tp_aln.vseq = vseq;
  tp_aln.vlen = vlen;
  tp_aln.start1 = start1;
  tp_aln.end1 = end1;
  tp_aln.start2 = start2;
  tp_aln.end2 = end2;
  tp_aln.delta = delta;
  tp_aln.cigar = cig;
  tp_aln.ciglen = ciglen;

  // print struct
  printf("\nSequenz 1: \t\t%s\n", tp_aln.useq);
  printf("Länge Sequenz 1: \t%d\n", tp_aln.ulen);
  printf("Sequenz 2: \t\t%s\n", tp_aln.vseq);
  printf("Länge Sequenz 2: \t%d\n", tp_aln.vlen);
  printf("Start 1: \t\t%d\n", tp_aln.start1);
  printf("Ende 1: \t\t%d\n", tp_aln.end1);
  printf("Start 2: \t\t%d\n", tp_aln.start2);
  printf("Ende 2: \t\t%d\n", tp_aln.end2);
  printf("Delta: \t\t\t%d\n", tp_aln.delta);
  printf("CIGAR: \t\t\t%s\n", tp_aln.cigar);
  printf("Länge CIGAR: \t\t%d\n", tp_aln.ciglen);
  printf("\n");
  
  encode(tp_aln);
  */
  gt_tracepoint_data_delete(tp_data);
  return 0;
}
