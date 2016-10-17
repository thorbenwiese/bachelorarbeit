#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "eoplist.h"
#include "TracePoint.h"
#include "gt-alloc.h"
#include "front-with-trace.h"

int main(int argc, char *argv[])
{

  GtUchar *long_useq, *long_vseq;
  GtUword start1, end1, start2, end2, delta, edist;

  long readstart1, readend1, readstart2, readend2, readdelta;
  bool haserr = false;

  bool decode = true;
  bool store = false;

  if (argc != 8)
  {
    fprintf(stderr,"Usage: %s <sequence1> <sequence2> <start1> <end1> "
                             "<start2> <end2> <delta>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  /* read and check parameters */
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
    GtEoplist *eoplist;
    FrontEdistTrace *fet = NULL;
    TracePointList *tp_list = NULL;
    GtUchar *useq = NULL, *vseq = NULL;
    char *enc_cigar = NULL, *dec_cigar = NULL;
    GtUword unitcost;

    struct timeval start_time;
    struct timeval comp_time;

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
    eoplist = gt_eoplist_new();
    fet = front_edist_trace_new();
    edist = front_edist_trace_eoplist(eoplist,
                                      fet,
                                      useq,
                                      end1 - start1 + 1,
                                      vseq,
                                      end2 - start2 + 1,
                                      true);
    front_edist_trace_delete(fet);
    unitcost = gt_eoplist_unit_cost(eoplist);
    assert(edist == unitcost);

    gettimeofday(&start_time, NULL);
    gt_tracepoint_encode(tp_list, eoplist);
    gettimeofday(&comp_time, NULL); 

    enc_cigar = gt_eoplist2cigar_string(eoplist,false);
    //printf("CIGAR Encode: %s\n", enc_cigar);
    gt_free(enc_cigar);
    
    double enc_time = (comp_time.tv_sec - start_time.tv_sec) + 
                      (comp_time.tv_usec - start_time.tv_usec) * 1e-6; 
    printf("Berechnungszeit Encode: %.6f\n", enc_time);

    if(store)
    {
      FILE *f1 = fopen("enc_time.txt", "a");
      if (f1 == NULL)
      {
        printf("Error opening file!\n");
        exit(EXIT_FAILURE);
      }

      fprintf(f1, "%f\n", enc_time);

      fclose(f1);
    }

    /* decode TracePoint Array and TracePointData to GtEoplist */
    if(decode)
    {
      GtEoplist *eoplist_tp = gt_tracepoint_decode(tp_list);
      GtUword unitcost_dc;

      gettimeofday(&start_time, NULL);
      dec_cigar = gt_eoplist2cigar_string(eoplist_tp,false);
      gettimeofday(&comp_time, NULL); 

      //printf("CIGAR Decode: %s\n", dec_cigar);
      gt_free(dec_cigar);
      unitcost_dc = gt_eoplist_unit_cost(eoplist_tp);
      if (unitcost_dc > unitcost)
      {
        fprintf(stderr,"unitcost_decode = %lu > %lu = unitcost_encode\n", 
                 unitcost_dc,unitcost);
        exit(EXIT_FAILURE);
      }
      gt_eoplist_delete(eoplist_tp);
      double dec_time = (comp_time.tv_sec - start_time.tv_sec) + 
                        (comp_time.tv_usec - start_time.tv_usec) * 1e-6;
      printf("Berechnungszeit Decode: %.6f\n", dec_time);

      if(store)
      {
        FILE *f2 = fopen("dec_time.txt", "a");
        if (f2 == NULL)
        {
          printf("Error opening file!\n");
          exit(EXIT_FAILURE);
        }

        fprintf(f2, "%f\n", dec_time);

        fclose(f2);
      }
    }
    gt_eoplist_delete(eoplist);
    gt_tracepoint_list_delete(tp_list);
  }

  gt_free(long_useq);
  gt_free(long_vseq);

  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}
