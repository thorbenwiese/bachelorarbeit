#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "eoplist.h"
#include "TracePoint.h"
#include "gt-alloc.h"

int main(int argc, char *argv[])
{

  GtUchar *long_useq, *long_vseq, *useq = NULL, *vseq = NULL;
  GtUword start1, end1, start2, end2, delta;
  TracePointList *tp_list = NULL;
  GtEoplist *eoplist = NULL;
  GtEoplistReader *eoplist_reader = NULL;
  GtCigarOp co;
  long readstart1, readend1, readstart2, readend2, readdelta;
  bool haserr = false;

  if (argc != 8)
  {
    fprintf(stderr,"Usage: %s <sequence1> <sequence2> <start1> <end1> "
                             "<start2> <end2> <delta>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  /* read and check parameters */
  //TODO check if characters
  long_useq = (unsigned char *) strdup(argv[1]);
  long_vseq = (unsigned char *) strdup(argv[2]);

  if (sscanf(argv[3],"%ld",&readstart1) != 1 || readstart1 < 0)
  {
    fprintf(stderr,"%s: <start1> must be non-negative number\n",argv[0]);
    haserr = true;
  }
  else if (sscanf(argv[4],"%ld",&readend1) != 1 || readend1 < readstart1)
  {
    fprintf(stderr,"%s: <end1> must be greater than <start1>\n",argv[0]);
    haserr = true;
  }
  else if (sscanf(argv[5],"%ld",&readstart2) != 1 || readstart2 < 0)
  {
    fprintf(stderr,"%s: <start2> must be non-negative number\n",argv[0]);
    haserr = true;
  }
  else if (sscanf(argv[6],"%ld",&readend2) != 1 || readend2 < readstart2)
  {
    fprintf(stderr,"%s: <end2> must be greater than <start2>\n",argv[0]);
    haserr = true;
  }
  else if (sscanf(argv[7],"%ld",&readdelta) != 1 || readdelta < 0)
  {
    fprintf(stderr,"%s: <delta> must be non-negative number\n",argv[0]);
    haserr = true;
  }
  if (!haserr)
  {
    start1 = (GtUword) readstart1;
    end1 = (GtUword) readend1;
    start2 = (GtUword) readstart2;
    end2 = (GtUword) readend2;
    delta = (GtUword) readdelta;

    /* get substring with regard to start of sequence */
    useq = long_useq + start1;
    vseq = long_vseq + start2;

    /* create TracePointData */
    tp_list = gt_tracepoint_list_new();
    gt_tracepoint_list_set(tp_list, useq, vseq, start1, end1, 
                           start2, end2, delta);

    /* encode list to Trace Point Array */
    printf("ENCODE\n");
    eoplist = gt_eoplist_new();
    gt_tracepoint_encode(tp_list, eoplist);

    /* print Trace Points */
    gt_print_tracepoint_list(tp_list);

    /* decode TracePoint Array and TracePointData to GtEoplist */
    printf("DECODE\n");
    gt_eoplist_reset(eoplist);
    eoplist = gt_tracepoint_decode(tp_list);
    while (gt_eoplist_reader_next_cigar(&co, eoplist_reader))                    
    {                                                                            
      printf("%lu%c",co.iteration, gt_eoplist_pretty_print(co.eoptype, false));  
    }      
  }

  gt_tracepoint_list_delete(tp_list);
  gt_eoplist_delete(eoplist);
  gt_free(long_useq);
  gt_free(long_vseq);

  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}
