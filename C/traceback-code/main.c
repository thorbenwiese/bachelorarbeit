#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "gt-defs.h"
#include "gt-alloc.h"
#include "front-with-trace.h"

int main(int argc,char *argv[])
{
  unsigned char *useq, *vseq;
  GtUword ulen, vlen, edist;
  FrontEdistTrace *fet;
  GtEoplist *eoplist, *eoplist2;
  GtEoplistReader *eoplist_reader;
  GtCigarOp co;
  char *cigar_string, *cigar_string2;
  bool distinguish_mismatch_match = false;
  int argc_start_seq;

  if (argc != 3 && argc != 4)
  {
    fprintf(stderr,"Usage: %s [-d] <sequence1> <sequence2>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  if (argc == 3)
  {
    argc_start_seq = 1;
  } else
  {
    if (strcmp(argv[1],"-d") == 0)
    {
      distinguish_mismatch_match = true;
    } else
    {
      fprintf(stderr,"Usage: %s [-d] <sequence1> <sequence2>\n",argv[0]);
      exit(EXIT_FAILURE);
    }
    argc_start_seq = 2;
  }
  useq = (unsigned char *) strdup(argv[argc_start_seq]);
  ulen = strlen(argv[argc_start_seq]);
  vseq = (unsigned char *) strdup(argv[argc_start_seq+1]);
  vlen = strlen(argv[argc_start_seq+1]);
  fet = front_edist_trace_new();
  eoplist = gt_eoplist_new();
  edist = front_edist_trace_eoplist(eoplist,
                                    fet,
                                    useq,
                                    ulen,
                                    vseq,
                                    vlen,
                                    true);
  assert(edist == gt_eoplist_unit_cost(eoplist));
  eoplist_reader = gt_eoplist_reader_new(eoplist);
  if (distinguish_mismatch_match)
  {
    gt_eoplist_reader_distinguish_mismatch_match(eoplist_reader);
  }
  while (gt_eoplist_reader_next_cigar(&co,eoplist_reader))
  {
    printf("%lu%c",co.iteration,
                   gt_eoplist_pretty_print(co.eoptype,
                                           distinguish_mismatch_match));
  }
  printf("\n");
  cigar_string = gt_eoplist2cigar_string(eoplist,distinguish_mismatch_match);
  printf("%s\n",cigar_string);
  eoplist2 = gt_eoplist_new_from_cigar(cigar_string,strlen(cigar_string));
  cigar_string2 = gt_eoplist2cigar_string(eoplist2,distinguish_mismatch_match);
  assert(strcmp(cigar_string,cigar_string2) == 0);
  front_edist_trace_delete(fet);
  gt_eoplist_reader_reset(eoplist_reader,eoplist,70);
  gt_eoplist_format_generic(stdout, eoplist_reader, useq, ulen, vseq, vlen,
                            distinguish_mismatch_match);
  gt_eoplist_reader_delete(eoplist_reader);
  gt_eoplist_delete(eoplist);
  gt_free(cigar_string);
  gt_free(cigar_string2);
  gt_eoplist_delete(eoplist2);
  gt_free(useq);
  gt_free(vseq);
  return EXIT_SUCCESS;
}
