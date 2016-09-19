#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
// ctype.h for isdigit()
#include <ctype.h>

#include "TracePoint.h"
#include "Alignment.h"
#include "Cigar_Pattern.h"

int main(int argc, char *argv[]){

  unsigned char *useq, *vseq, *cig;

  if (argc != 9){
    fprintf(stderr,"Usage: %s <sequence1> <sequence2> <start1> <end1> <start2> <end2> <delta> <cigar>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  useq = (unsigned char *) strdup(argv[1]);
  int ulen = strlen(argv[1]);
  vseq = (unsigned char *) strdup(argv[2]);
  int vlen = strlen(argv[2]);
  int start1 = atoi(argv[3]);
  int end1 = atoi(argv[4]);
  int start2 = atoi(argv[5]);
  int end2 = atoi(argv[6]);
  int delta = atoi(argv[7]);
  cig = (unsigned char *) strdup(argv[8]);
  int ciglen = strlen(argv[8]);
  
  // struct to store input
  struct TracePointAlignment tp_aln;

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
  
  return 0;
}
