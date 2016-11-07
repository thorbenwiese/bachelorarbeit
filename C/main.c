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

#define OPTIONS "hfpdax"

/* struct to store input parameters: 
   input sequence file, start/end positions, delta*/
typedef struct
{
  char *inputfile;
  GtUword start1, end1, start2, end2, delta, amount;
  bool decode;
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
  options->decode = false;

  if (argc == 1)
  {
    printf("Usage:\n"
           "\t-f <inputfile>\n"
           "\t-a <amount of sequence pairs read from file>\n"
           "\t-p <positions of sequences: start1 end1 start2 end2>\n"
           "\t-d <delta value>\n"
           "\t-x <decode Trace Points>\n");
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
               "\t-f <inputfile>\n"
               "\t-a <amount of sequence pairs read from file>\n"
               "\t-p <positions of sequences: start1 end1 start2 end2>\n"
               "\t-d <delta value>\n"
               "\t-x <decode Trace Points>\n");
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
        else if (sscanf(argv[optind + 1],"%ld",&readend1) != 1 || 
                 readend1 < readstart1)
        {
          fprintf(stderr,"%s: <end1> must be greater than <start1>\n",argv[0]);
          exit(EXIT_FAILURE);
        }
        else if (sscanf(argv[optind + 2],"%ld",&readstart2) != 1 || 
                 readstart2 < 0)
        {
          fprintf(stderr,"%s: <start2> must be non-negative number\n",argv[0]);
          exit(EXIT_FAILURE);
        }
        else if (sscanf(argv[optind + 3],"%ld",&readend2) != 1 || 
                 readend2 < readstart2)
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
      case 'x':
        options->decode = true;
        break;
      default:
        printf("Usage:\n"
               "-f <inputfile>\n"
               "-a <amount of sequence pairs read from file>\n"
               "-p <positions of sequences: start1 end1 start2 end2>\n"
               "-d <delta value>\n"
               "-x <decode Trace Points>\n");
    }
  }
}

int main(int argc, char *argv[])
{
  Options options;

  parse_options(&options, argc, argv);

  if(options.inputfile != NULL)
  {
    GtEoplist *eoplist;
    FrontEdistTrace *fet = NULL;
    TracePointList *tp_list = NULL;
    GtUchar *useq = NULL, *vseq = NULL;
    GtUword unitcost, edist;
    FILE *inf;

    GtUword len = 10000, i;
    char *line1 = malloc(len);
    char *line2 = malloc(len);

    inf = fopen(options.inputfile, "r");
    if(inf == NULL)
    {
      fprintf(stderr,"%s: cannot read file \"%s\"\n",argv[0],options.inputfile);
      exit(EXIT_FAILURE);
    }

    for(i = 0; i < options.amount; i++)
    {
      fgets(line1, len, inf);
      useq = (unsigned char *) line1 + options.start1;
      fgets(line2, len, inf);
      vseq = (unsigned char *) line2 + options.start2;

      //printf("%s\n\n%s\n",useq,vseq);
      /* create TracePointData */
      tp_list = gt_tracepoint_list_new();
      gt_tracepoint_list_set(tp_list,
                             useq, 
                             vseq, 
                             options.start1,
                             options.end1,
                             options.start2,
                             options.end2,
                             options.delta);

      /* encode list to Trace Point Array */
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

      gt_tracepoint_encode(tp_list, eoplist);

      gt_eoplist_delete(eoplist);
    
      if(options.decode)
      {
        GtEoplist *eoplist_tp = gt_tracepoint_decode(tp_list);
          
        GtUword unitcost_dc;
        unitcost_dc = gt_eoplist_unit_cost(eoplist_tp);

        gt_eoplist_delete(eoplist_tp);

        if (unitcost_dc > unitcost)
        {
          fprintf(stderr,"unitcost_decode = %lu > %lu = unitcost_encode\n", 
                  unitcost_dc,unitcost);
          exit(EXIT_FAILURE);
        }
      }
    }
    
    gt_tracepoint_list_delete(tp_list);

    gt_free(line1);
    gt_free(line2);
    fclose(inf);
  }
  return EXIT_SUCCESS;
}
