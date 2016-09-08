#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <math.h>
// ctype.h for isdigit()
#include <ctype.h>

#include "TracePoint.h"
#include "Alignment.h"
#include "Cigar_Pattern.h"

# define LEN 100000

char seq1_arr[LEN];
char seq2_arr[LEN];
char cigar_arr[LEN];

int main(){

  // struct to store input
  struct TracePointAlignment tp_aln;
  printf("Enter Seq1, Seq2, Start1, End1, Start2 ,End2, Delta, Cigar:\n");

  // read input and store integers in struct and strings in char arrays
  scanf("%s %s %d %d %d %d %d %s", seq1_arr, seq2_arr, 
                                   &tp_aln.start_seq1, &tp_aln.end_seq1,
                                   &tp_aln.start_seq2, &tp_aln.end_seq2,
                                   &tp_aln.delta, cigar_arr);

  // check input
  assert(seq1_arr[0] != 0);
  assert(seq2_arr[0] != 0);
  assert(tp_aln.start_seq1 >= 0);
  assert(tp_aln.end_seq1 >= 1);
  assert(tp_aln.start_seq2 >= 0);
  assert(tp_aln.end_seq2 >= 1);
  assert(tp_aln.delta > 0);
  assert(cigar_arr[0] != 0);

  // pointer to malloc storage for seq1, seq2, cigar
  char *seq1;
  char *seq2;
  char *cigar;

  // len of seq1, seq2, cigar + null character
  int seq1len = strlen(seq1_arr) + 1;
  int seq2len = strlen(seq2_arr) + 1;
  int cigarlen = strlen(cigar_arr) + 1;

  // copy seq1, seq2, cigar into dynamically allocated storage
  seq1 = (char*)malloc(seq1len * sizeof(char));
  memcpy(seq1, &seq1_arr[tp_aln.start_seq1], 
         tp_aln.end_seq1 - tp_aln.start_seq1 + 1);
  seq2 = (char*)malloc(seq2len * sizeof(char));
  memcpy(seq2, &seq2_arr[tp_aln.start_seq2], 
         tp_aln.end_seq2 - tp_aln.start_seq2 + 1);
  cigar = (char*)malloc(cigarlen * sizeof(char));
  strcpy(cigar, cigar_arr);

  // store in struct
  tp_aln.seq1 = seq1;
  tp_aln.seq2 = seq2;
  tp_aln.cigar = cigar;

  // print struct
  printf("\nSequenz 1: \t%s\n", tp_aln.seq1);
  printf("Sequenz 2: \t%s\n", tp_aln.seq2);
  printf("Start 1: \t%d\n", tp_aln.start_seq1);
  printf("Ende 1: \t%d\n", tp_aln.end_seq1);
  printf("Start 2: \t%d\n", tp_aln.start_seq2);
  printf("Ende 2: \t%d\n", tp_aln.end_seq2);
  printf("Delta: \t\t%d\n", tp_aln.delta);
  printf("CIGAR: \t\t%s\n", tp_aln.cigar);
  printf("\n");

  // encode alignment
  encode(tp_aln);
  
  return 0;
}
