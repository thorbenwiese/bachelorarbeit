#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "eoplist.h"
#include "TracePoint.h"
#include "gt-alloc.h"
#include "front-with-trace.h"

#define OPTIONS "hfpdaex"

/* struct to store input parameters: 
   input sequence file, start/end positions, delta*/
typedef struct
{
  char *inputfile;
  GtUword start1, end1, start2, end2, delta, amount;
  bool encode, decode;
} Options;

static void parse_options(Options *options, int argc, char * const argv[])
{
  long readstart1, readend1, readstart2, readend2, readdelta, readamount;
  options->inputfile = NULL;
  options->start1 = 0;
  options->end1 = 0;
  options->start2 = 0;
  options->end2 = 0;
  options->delta = 0;
  options->amount = 0;
  options->encode = false;
  options->decode = false;

  if (argc == 1)
  {
    fprintf(stderr,"Usage: %s <input file> <start1> <end1> "
                             "<start2> <end2> <delta>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  while (true)
  {
    int opt = getopt(argc, argv, OPTIONS);
    if (opt == -1)
    {
      break;
    }
    switch((char) opt)
    {
      case 'h':
        printf("Usage:\n"
               "-f <inputfile>\n"
               "-p <positions of sequences: start1 end1 start2 end2>\n"
               "-d <delta value>\n"
               "-e <encode Trace Points>\n"
               "-x <decode Trace Points>\n");
      case 'f':
        options->inputfile = argv[optind];
        optind++;
        break;
      case 'p':
        if (optind + 3 >= argc)
        {
          fprintf(stderr,"%s: missing argument to option -%c\n",argv[0],
                  (char) opt);
          exit(EXIT_FAILURE);
        }
        if (sscanf(argv[optind],"%ld",&readstart1) != 1 || readstart1 < 0)
        {
          fprintf(stderr,"%s: <start1> must be non-negative number\n",argv[0]);
          exit(EXIT_FAILURE);
        }
        else if (sscanf(argv[optind + 1],"%ld",&readend1) != 1 || readend1 < readstart1)
        {
          fprintf(stderr,"%s: <end1> must be greater than <start1>\n",argv[0]);
          exit(EXIT_FAILURE);
        }
        else if (sscanf(argv[optind + 2],"%ld",&readstart2) != 1 || readstart2 < 0)
        {
          fprintf(stderr,"%s: <start2> must be non-negative number\n",argv[0]);
          exit(EXIT_FAILURE);
        }
        else if (sscanf(argv[optind + 3],"%ld",&readend2) != 1 || readend2 < readstart2)
        {
          fprintf(stderr,"%s: <end2> must be greater than <start2>\n",argv[0]);
          exit(EXIT_FAILURE);
        }
        options->start1 = (GtUword) readstart1;
        options->end1 = (GtUword) readend1;
        options->start2 = (GtUword) readstart2;
        options->end2 = (GtUword) readend2;
        optind += 4;
        break;
      case 'd':
        if (sscanf(argv[optind],"%ld",&readdelta) != 1 || readdelta < 0)
        {
          fprintf(stderr,"%s: <delta> must a be non-negative number\n",argv[0]);
          exit(EXIT_FAILURE);
        }
        options->delta = (GtUword) readdelta;
        optind++;
        break;
      case 'a':
        if (sscanf(argv[optind],"%ld",&readamount) != 1 || readamount < 1)
        {
          fprintf(stderr,"%s: <amount> must be greater than 0\n",argv[0]);
          exit(EXIT_FAILURE);
        }
        options->amount = (GtUword) readamount;
        optind++;
        break;
      case 'e':
        options->encode = true;
        break;
      case 'x':
        options->decode = true;
        break;
      default:
        printf("Usage:\n"
               "-f <inputfile>\n"
               "-p <positions of sequences: start1 end1 start2 end2>\n"
               "-d <delta value>\n"
               "-e <encode Trace Points>\n"
               "-x <decode Trace Points>\n");
    }
  }
  
}

int main(int argc, char *argv[])
{
  Options options;
  GtEoplist *eoplist;
  FrontEdistTrace *fet = NULL;
  TracePointList *tp_list = NULL;
  GtUchar *useq = NULL, *vseq = NULL;
  //char *enc_cigar = NULL, *dec_cigar = NULL;
  GtUword unitcost, edist;
  FILE *inf;//, *outf;

  struct timeval start_enc_time;
  struct timeval comp_enc_time;

  struct timeval start_dec_time;
  struct timeval comp_dec_time;

  double total_enc_time = 0.0, total_dec_time = 0.0;
  
  parse_options(&options, argc, argv);

  if(options.inputfile != NULL)
  {
    GtUword len = 10000, i;
    char *line1 = malloc(len);
    char *line2 = malloc(len);
    inf = fopen(options.inputfile, "r");
    if(inf == NULL)
    {
      fprintf(stderr,"%s: cannot read file \"%s\"\n",argv[0],options.inputfile);
      exit(EXIT_FAILURE);
    }
    for(i = 1; i < options.amount * 2; i++)
    {
      fgets(line1, len, inf);
      useq = (unsigned char *) line1 + options.start1;
      fgets(line2, len, inf);
      vseq = (unsigned char *) line2 + options.start2;

      // create TracePointData
      tp_list = gt_tracepoint_list_new();
      gt_tracepoint_list_set(tp_list,
                             useq, 
                             vseq, 
                             options.start1,
                             options.end1,
                             options.start2,
                             options.end2,
                             options.delta);

      // encode list to Trace Point Array
      eoplist = gt_eoplist_new();
      fet = front_edist_trace_new();
      edist = front_edist_trace_eoplist(eoplist,
                                        fet,
                                        useq,
                                        options.end1 - options.start1 + 1,
                                        vseq,
                                        options.end2 - options.start2 + 1,
                                        true);
      front_edist_trace_delete(fet);
      unitcost = gt_eoplist_unit_cost(eoplist);
      assert(edist == unitcost);

      gettimeofday(&start_enc_time, NULL);
      gt_tracepoint_encode(tp_list, eoplist);
      gettimeofday(&comp_enc_time, NULL); 
      //store to TP File

      //enc_cigar = gt_eoplist2cigar_string(eoplist,false);
      //printf("ENCODE CIGAR: %s\n",enc_cigar);
      //gt_free(enc_cigar);
      gt_eoplist_delete(eoplist);
    
      double enc_time = (comp_enc_time.tv_sec - start_enc_time.tv_sec) + 
                        (comp_enc_time.tv_usec - start_enc_time.tv_usec) * 1e-6; 

      total_enc_time += enc_time;
      //printf("Berechnungszeit Encode: %.6f\n", enc_time);

      if(options.decode)
      {
        gettimeofday(&start_dec_time, NULL);
        GtEoplist *eoplist_tp = gt_tracepoint_decode(tp_list);
        gettimeofday(&comp_dec_time, NULL); 
        GtUword unitcost_dc;

        //dec_cigar = gt_eoplist2cigar_string(eoplist_tp,false);
        //printf("DECODE CIGAR: %s\n",dec_cigar);

        //gt_free(dec_cigar);
        unitcost_dc = gt_eoplist_unit_cost(eoplist_tp);
        if (unitcost_dc > unitcost)
        {
          fprintf(stderr,"unitcost_decode = %lu > %lu = unitcost_encode\n", 
                  unitcost_dc,unitcost);
          exit(EXIT_FAILURE);
        }
        gt_eoplist_delete(eoplist_tp);

        double dec_time = (comp_dec_time.tv_sec - start_dec_time.tv_sec) + 
                          (comp_dec_time.tv_usec - start_dec_time.tv_usec)*1e-6;

        total_dec_time += dec_time;
        //printf("Berechnungszeit Decode: %.6f\n", dec_time);
      }
    }
    //printf("Gesamte Berechnungszeit Encode: %.6f\n", total_enc_time);
    printf("%.8f,", total_dec_time);
    //printf("%.6f\n", (float) total_dec_time/options.amount);
    gt_free(line1);
    gt_free(line2);
    fclose(inf);
    gt_tracepoint_list_delete(tp_list);
  }
  return EXIT_SUCCESS;
}
