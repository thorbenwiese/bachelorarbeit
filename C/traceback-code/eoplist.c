#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
/*#include "core/assert_api.h"*/
#include "gt-defs.h"
#include "gt-alloc.h"
#include "eoplist.h"

#define DELETION_CHAR    'D'
#define INSERTION_CHAR   'I'
#define MATCH_CHAR       '='
#define MISMATCH_CHAR    'X'
#define REPLACEMENT_CHAR 'M'

char gt_eoplist_pretty_print(GtEopType eoptype,bool distinguish_mismatch_match)
{
  switch(eoptype)
  {
    case GtDeletionOp:
      return DELETION_CHAR;
    case GtInsertionOp:
      return INSERTION_CHAR;
    case GtMismatchOp:
      return distinguish_mismatch_match ? MISMATCH_CHAR : REPLACEMENT_CHAR;
    case GtMatchOp:
      return distinguish_mismatch_match ? MATCH_CHAR : REPLACEMENT_CHAR;
    default:
      gt_assert(false);
  }
}

struct GtEoplist
{
  GtUword nextfreeuint8_t, allocateduint8_t, countmismatches, countmatches,
                                             countdeletions, countinsertions;
  uint8_t *spaceuint8_t;
};

void gt_eoplist_reset(GtEoplist *eoplist)
{
  if(eoplist != NULL)
  {
    eoplist->nextfreeuint8_t = 0;
    eoplist->countmismatches = 0;
    eoplist->countmatches = 0;
    eoplist->countdeletions = 0;
    eoplist->countinsertions = 0;
  }
}

GtEoplist *gt_eoplist_new(void)
{
  GtEoplist *eoplist = gt_malloc(sizeof *eoplist);

  gt_assert(eoplist != NULL);
  eoplist->allocateduint8_t = 0;
  eoplist->spaceuint8_t = NULL;
  gt_eoplist_reset(eoplist);
  return eoplist;
}

typedef struct
{
  char *space;
  size_t nextfree, allocated;
} GtStringBuffer;

static void stringbuffer_append_cigar(GtStringBuffer *sbuf,
                                      const GtCigarOp *co,
                                      bool distinguish_mismatch_match)
{
  const size_t gt_uword_maxwidth = sizeof("18446744073709551615");

  if (sbuf->nextfree + gt_uword_maxwidth + 1 + 1 >= sbuf->allocated)
  {
    sbuf->allocated = sbuf->allocated * 1.2 + gt_uword_maxwidth + 1 + 1 + 1;
    sbuf->space = gt_realloc(sbuf->space,sizeof *sbuf->space * sbuf->allocated);
  }
  sbuf->nextfree +=
    sprintf(sbuf->space + sbuf->nextfree,"%lu%c",
            co->iteration,gt_eoplist_pretty_print(co->eoptype,
                                                  distinguish_mismatch_match));
}

char *gt_eoplist2cigar_string(const GtEoplist *eoplist,
                              bool distinguish_mismatch_match)
{
  GtStringBuffer sbuf = {NULL,0,0};
  GtCigarOp co;
  GtEoplistReader *eoplist_reader = gt_eoplist_reader_new(eoplist);

  if (distinguish_mismatch_match)
  {
    gt_eoplist_reader_distinguish_mismatch_match(eoplist_reader);
  }
  while (gt_eoplist_reader_next_cigar(&co,eoplist_reader))
  {
    stringbuffer_append_cigar(&sbuf,&co,distinguish_mismatch_match);
  }
  gt_eoplist_reader_delete(eoplist_reader);
  return sbuf.space;
}

GtEoplist *gt_eoplist_new_from_cigar(const char *cigarstring,GtUword length)
{
  const char *cptr;
  GtUword iteration = 0;
  GtEoplist *eoplist = gt_eoplist_new();

  for (cptr = cigarstring; cptr < cigarstring + length; cptr++)
  {
    if (isdigit(*cptr))
    {
      iteration = iteration * 10 + (GtUword) (*cptr - '0');
    } else
    {
      unsigned long idx;
      switch (*cptr)
      {
        case DELETION_CHAR:
          for (idx = 0; idx < iteration; idx++)
          {
            gt_eoplist_deletion_add(eoplist);
          }
          break;
        case INSERTION_CHAR:
          for (idx = 0; idx < iteration; idx++)
          {
            gt_eoplist_insertion_add(eoplist);
          }
          break;
        case MATCH_CHAR:
        case REPLACEMENT_CHAR:
          gt_eoplist_match_add(eoplist,iteration);
          break;
        case MISMATCH_CHAR:
          for (idx = 0; idx < iteration; idx++)
          {
            gt_eoplist_mismatch_add(eoplist);
          }
          break;
        default:
          fprintf(stderr,"Illegal symbol %c in Cigar string\n",*cptr);
          exit(EXIT_FAILURE);
      }
      iteration = 0;
    }
  }
  return eoplist;
}

void gt_eoplist_delete(GtEoplist *eoplist)
{
  if (eoplist != NULL)
  {
    gt_free(eoplist->spaceuint8_t);
    gt_free(eoplist);
  }
}

#define FT_EOPCODE_MAXMATCHES     253
#define FT_EOPCODE_MISMATCH       253
#define FT_EOPCODE_DELETION       254
#define FT_EOPCODE_INSERTION      255

#define GT_CHECKARRAYSPACE(A,TYPE,L)\
        do {\
          if ((A)->nextfree##TYPE >= (A)->allocated##TYPE)\
          {\
            (A)->allocated##TYPE += L;\
            (A)->space##TYPE = (TYPE *) gt_realloc((A)->space##TYPE,\
                                                   sizeof (TYPE) *\
                                                   (A)->allocated##TYPE);\
          }\
          gt_assert((A)->space##TYPE != NULL);\
        } while (false)

#define GT_STOREINARRAY(A,TYPE,L,VAL)\
        do {\
          GT_CHECKARRAYSPACE(A,TYPE,L);\
          (A)->space##TYPE[(A)->nextfree##TYPE++] = VAL;\
        } while (false)

#define GT_EOPLIST_PUSH(EOPLIST,EOP)\
        gt_assert((EOPLIST) != NULL);\
        {\
          const GtUword addamount = (EOPLIST)->allocateduint8_t * 0.2 + 128;\
          GT_STOREINARRAY(EOPLIST,uint8_t,addamount,(uint8_t) (EOP));\
        }

void gt_eoplist_match_add(GtEoplist *eoplist,GtUword length)
{
  gt_assert(eoplist != NULL && length > 0);
  eoplist->countmatches += length;
  while (true)
  {
    if (length <= FT_EOPCODE_MAXMATCHES)
    {
      gt_assert(length > 0);
      GT_EOPLIST_PUSH(eoplist,(uint8_t) (length - 1)); /* R length */
      break;
    }
    GT_EOPLIST_PUSH(eoplist,FT_EOPCODE_MAXMATCHES - 1); /* R max */
    length -= FT_EOPCODE_MAXMATCHES;
  }
}

void gt_eoplist_mismatch_add(GtEoplist *eoplist)
{
  GT_EOPLIST_PUSH(eoplist,FT_EOPCODE_MISMATCH); /* R 1 */
  eoplist->countmismatches++;
}

void gt_eoplist_deletion_add(GtEoplist *eoplist)
{
  GT_EOPLIST_PUSH(eoplist,FT_EOPCODE_DELETION);
  eoplist->countdeletions++;
}

void gt_eoplist_insertion_add(GtEoplist *eoplist)
{
  gt_assert(eoplist != NULL);
  GT_EOPLIST_PUSH(eoplist,FT_EOPCODE_INSERTION);
  eoplist->countinsertions++;
}

GtUword gt_eoplist_length(const GtEoplist *eoplist)
{
  if (eoplist == NULL)
  {
    return 0;
  }
  return eoplist->nextfreeuint8_t;
}

void gt_eoplist_reverse_end(GtEoplist *eoplist,GtUword firstindex)
{
  uint8_t *fwd, *bck;

  gt_assert(eoplist != NULL);
  if (firstindex + 1 >= eoplist->nextfreeuint8_t)
  {
    return;
  }
  for (fwd = eoplist->spaceuint8_t + firstindex,
       bck = eoplist->spaceuint8_t + eoplist->nextfreeuint8_t - 1; fwd < bck;
       fwd++, bck--)
  {
    uint8_t tmp = *fwd;
    *fwd = *bck;
    *bck = tmp;
  }
}

GtUword gt_eoplist_matches_count(const GtEoplist *eoplist)
{
  gt_assert(eoplist != NULL);
  return eoplist->countmatches;
}

GtUword gt_eoplist_mismatches_count(const GtEoplist *eoplist)
{
  gt_assert(eoplist != NULL);
  return eoplist->countmismatches;
}

GtUword gt_eoplist_deletions_count(const GtEoplist *eoplist)
{
  gt_assert(eoplist != NULL);
  return eoplist->countdeletions;
}

GtUword gt_eoplist_insertions_count(const GtEoplist *eoplist)
{
  gt_assert(eoplist != NULL);
  return eoplist->countinsertions;
}

GtUword gt_eoplist_unit_cost(const GtEoplist *eoplist)
{
  return gt_eoplist_insertions_count(eoplist) +
         gt_eoplist_deletions_count(eoplist) +
         gt_eoplist_mismatches_count(eoplist);
}

struct GtEoplistReader
{
  const uint8_t *endeoplist,
                *currenteop;
  bool distinguish_mismatch_match;
  GtUchar *outbuffer;
  unsigned int width;
  GtUword aligned_u, aligned_v, repcount;
};

void gt_eoplist_reader_reset(GtEoplistReader *eoplist_reader,
                             const GtEoplist *eoplist,unsigned int width)
{
  gt_assert(eoplist != NULL && eoplist_reader != NULL);
  if (eoplist->spaceuint8_t == NULL || eoplist->nextfreeuint8_t == 0)
  {
    eoplist_reader->currenteop = NULL;
    eoplist_reader->endeoplist = NULL;
  } else
  {
    eoplist_reader->currenteop = eoplist->spaceuint8_t;
    eoplist_reader->endeoplist = eoplist->spaceuint8_t +
                                 eoplist->nextfreeuint8_t;
  }
  if (eoplist_reader->width < width)
  {
    eoplist_reader->width = width;
    eoplist_reader->outbuffer = gt_realloc(eoplist_reader->outbuffer,
                                           sizeof *eoplist_reader->outbuffer *
                                           3 * (width+1));
    gt_assert(eoplist_reader->outbuffer != NULL);
  }
  eoplist_reader->aligned_u = eoplist_reader->aligned_v
                            = eoplist_reader->repcount = 0;
}

GtEoplistReader *gt_eoplist_reader_new(const GtEoplist *eoplist)
{
  GtEoplistReader *eoplist_reader;

  gt_assert(eoplist != NULL);
  eoplist_reader = gt_malloc(sizeof *eoplist_reader);
  gt_assert(eoplist_reader != NULL);
  eoplist_reader->width = 0;
  eoplist_reader->outbuffer = NULL;
  eoplist_reader->distinguish_mismatch_match = false;
  gt_eoplist_reader_reset(eoplist_reader,eoplist,0);
  return eoplist_reader;
}

void gt_eoplist_reader_distinguish_mismatch_match(
      GtEoplistReader *eoplist_reader)
{
  eoplist_reader->distinguish_mismatch_match = true;
}

void gt_eoplist_reader_delete(GtEoplistReader *eoplist_reader)
{
  if (eoplist_reader != NULL)
  {
    if (eoplist_reader->outbuffer != NULL)
    {
      gt_free(eoplist_reader->outbuffer);
    }
    gt_free(eoplist_reader);
  }
}

bool gt_eoplist_reader_next_cigar(GtCigarOp *cigar_op,
                                  GtEoplistReader *eoplist_reader)
{
  if (eoplist_reader->currenteop == NULL ||
      eoplist_reader->currenteop >= eoplist_reader->endeoplist)
  {
    return false;
  }
  cigar_op->eoptype = GtUndefinedOp;
  cigar_op->iteration = 0;
  while (true)
  {
    if (cigar_op->iteration > 0)
    {
      switch (*eoplist_reader->currenteop)
      {
        case FT_EOPCODE_DELETION:
          gt_assert(cigar_op->eoptype != GtUndefinedOp);
          if (cigar_op->eoptype == GtDeletionOp)
          {
            cigar_op->iteration++;
            eoplist_reader->currenteop++;
            break;
          }
          return true;
        case FT_EOPCODE_INSERTION:
          gt_assert(cigar_op->eoptype != GtUndefinedOp);
          if (cigar_op->eoptype == GtInsertionOp)
          {
            cigar_op->iteration++;
            eoplist_reader->currenteop++;
            break;
          }
          return true;
       case FT_EOPCODE_MISMATCH:
          gt_assert(cigar_op->eoptype != GtUndefinedOp);
          if (eoplist_reader->distinguish_mismatch_match)
          {
            if (cigar_op->eoptype == GtMismatchOp)
            {
              cigar_op->iteration++;
              eoplist_reader->currenteop++;
              break;
            } else
            {
              return true;
            }
          } else
          {
            if (cigar_op->eoptype == GtMatchOp)
            {
              cigar_op->iteration++;
              eoplist_reader->currenteop++;
              break;
            } else
            {
              return true;
            }
          }
          gt_assert(false);
       default:
          if (cigar_op->eoptype == GtMatchOp)
          {
            gt_assert(*eoplist_reader->currenteop < FT_EOPCODE_MAXMATCHES);
            cigar_op->iteration += (1UL + *eoplist_reader->currenteop);
            eoplist_reader->currenteop++;
          } else
          {
            return true;
          }
      }
    } else
    {
      switch(*eoplist_reader->currenteop)
      {
        case FT_EOPCODE_DELETION:
          cigar_op->eoptype = GtDeletionOp;
          cigar_op->iteration = 1UL;
          break;
        case FT_EOPCODE_INSERTION:
          cigar_op->eoptype = GtInsertionOp;
          cigar_op->iteration = 1UL;
          break;
        case FT_EOPCODE_MISMATCH:
          cigar_op->eoptype
            = eoplist_reader->distinguish_mismatch_match ? GtMismatchOp
                                                         : GtMatchOp;
          cigar_op->iteration = 1UL;
          break;
        default:
          cigar_op->eoptype = GtMatchOp;
          cigar_op->iteration = (1UL + *eoplist_reader->currenteop);
          break;
      }
      eoplist_reader->currenteop++;
    }
    if (eoplist_reader->currenteop >= eoplist_reader->endeoplist)
    {
      return true;
    }
  }
}

GtUword gt_eoplist_num_segments(const GtEoplist *eoplist,GtUword delta)
{
  GtUword length_u = gt_eoplist_matches_count(eoplist) +
                     gt_eoplist_mismatches_count(eoplist) +
                     gt_eoplist_deletions_count(eoplist);

  return length_u/delta + 1;
}

bool gt_eoplist_reader_next_segment(GtEoplistSegment *segment,
                                    GtEoplistReader *eoplist_reader,
                                    GtUword delta)
{
  while (true)
  {
    if (eoplist_reader->repcount > 0)
    {
      eoplist_reader->aligned_u++;
      eoplist_reader->aligned_v++;
      eoplist_reader->repcount--;
    } else
    {
      if (eoplist_reader->currenteop >= eoplist_reader->endeoplist)
      {
        break;
      }
      switch (*eoplist_reader->currenteop)
      {
        case FT_EOPCODE_DELETION:
          eoplist_reader->aligned_u++;
          break;
        case FT_EOPCODE_INSERTION:
          eoplist_reader->aligned_v++;
          break;
        case FT_EOPCODE_MISMATCH:
          eoplist_reader->aligned_u++;
          eoplist_reader->aligned_v++;
          break;
        default:
          eoplist_reader->aligned_u++;
          eoplist_reader->aligned_v++;
          eoplist_reader->repcount = (GtUword) *eoplist_reader->currenteop;
      }
      eoplist_reader->currenteop++;
    }
    if (eoplist_reader->aligned_u == delta)
    {
      segment->aligned_u = eoplist_reader->aligned_u;
      segment->aligned_v = eoplist_reader->aligned_v;
      eoplist_reader->aligned_u = eoplist_reader->aligned_v = 0;
      return true;
    }
  }
  if (eoplist_reader->aligned_v > 0 || eoplist_reader->aligned_u > 0)
  {
    gt_assert(eoplist_reader->repcount == 0);
    segment->aligned_u = eoplist_reader->aligned_u;
    segment->aligned_v = eoplist_reader->aligned_v;
    eoplist_reader->aligned_v = 0;
    eoplist_reader->aligned_u = 0;
    gt_assert(eoplist_reader->repcount == 0 &&
           eoplist_reader->currenteop >= eoplist_reader->endeoplist);
    return true;
  }
  return false;
}

double gt_eoplist_segments_entropy(const GtEoplist *eoplist,GtUword delta)
{
  GtEoplistReader *eoplist_reader = gt_eoplist_reader_new(eoplist);
  GtEoplistSegment segment;
  const GtUword max_value = 2 * delta + 1;
  GtUword segment_count = 0, idx,
          *segment_dist = gt_calloc(max_value + 1,sizeof *segment_dist);
  double entropy = 0.0;

  while (gt_eoplist_reader_next_segment(&segment,eoplist_reader,delta))
  {
    gt_assert(segment.aligned_v <= max_value);
    segment_dist[segment.aligned_v]++;
    segment_count++;
  }
  gt_assert(segment_count <= gt_eoplist_num_segments(eoplist,delta));
  gt_eoplist_reader_delete(eoplist_reader);
  for (idx = 0; idx <= max_value; idx++)
  {
    if (segment_dist[idx] > 0)
    {
      double prob = (double) segment_dist[idx]/segment_count;
      entropy += prob * log2(prob);
    }
  }
  gt_free(segment_dist);
  return entropy == 0.0 ? 0.0 : -entropy;
}

void gt_eoplist_show_plain(const GtEoplist *eoplist)
{
  GtUword idx;

  for (idx = 0; idx < eoplist->nextfreeuint8_t; idx++)
  {
    if (eoplist->spaceuint8_t[idx] == FT_EOPCODE_DELETION)
    {
      printf("D\n");
    } else
    {
      if (eoplist->spaceuint8_t[idx] == FT_EOPCODE_INSERTION)
      {
        printf("I\n");
      } else
      {
        if (eoplist->spaceuint8_t[idx] == FT_EOPCODE_MISMATCH)
        {
          printf("M\n");
        } else
        {
          printf("%d\n",eoplist->spaceuint8_t[idx]);
        }
      }
    }
  }
}

void gt_eoplist_reader_verify(const GtEoplist *eoplist,
                              const GtUchar *useq,
                              GtUword ulen,
                              const GtUchar *vseq,
                              GtUword vlen,
                              GtUword edist,
                              bool distinguish_mismatch_match)
{
  GtEoplistReader *eoplist_reader = gt_eoplist_reader_new(eoplist);
  GtCigarOp co;
  GtUword sumulen = 0, sumvlen = 0, sumdist = 0;

  if (distinguish_mismatch_match)
  {
    gt_eoplist_reader_distinguish_mismatch_match(eoplist_reader);
  }
  while (gt_eoplist_reader_next_cigar(&co,eoplist_reader))
  {
    if (co.eoptype == GtDeletionOp)
    {
      sumulen += co.iteration;
      sumdist += co.iteration;
    } else
    {
      if (co.eoptype == GtInsertionOp)
      {
        sumvlen += co.iteration;
        sumdist += co.iteration;
      } else
      {
        GtUword idx;
        if (co.eoptype == GtMismatchOp)
        {
          gt_assert(eoplist_reader->distinguish_mismatch_match);
          sumdist += co.iteration;
        }
        for (idx = 0; idx < co.iteration; idx++)
        {
          GtUchar a = useq[sumulen+idx],
                  b = vseq[sumvlen+idx];
          if (a == b)
          {
            gt_assert(co.eoptype == GtMatchOp);
          } else
          {
            gt_assert(!eoplist_reader->distinguish_mismatch_match ||
                   co.eoptype == GtMismatchOp);
          }
          if (!eoplist_reader->distinguish_mismatch_match)
          {
            if (a != b)
            {
              sumdist++;
            }
          }
        }
        sumulen += co.iteration;
        sumvlen += co.iteration;
      }
    }
  }
  if (ulen != sumulen)
  {
    fprintf(stderr,"ulen = %lu != %lu = sumulen\n",ulen,sumulen);
    exit(EXIT_FAILURE);
  }
  if (vlen != sumvlen)
  {
    fprintf(stderr,"vlen = %lu != %lu = sumvlen\n",vlen,sumvlen);
    exit(EXIT_FAILURE);
  }
  if (edist != sumdist)
  {
    fprintf(stderr,"edist = %lu != %lu = sumdist\n",edist,sumdist);
    exit(EXIT_FAILURE);
  }
  gt_eoplist_reader_delete(eoplist_reader);
}

static unsigned int gt_eoplist_show_advance(unsigned int pos,
                                            unsigned int width,
                                            const GtUchar *topbuf,
                                            FILE *fp)
{
  gt_assert(width > 0);
  if (pos < width - 1)
  {
    return pos + 1;
  }
  gt_assert(pos == width - 1);
  fwrite(topbuf,sizeof *topbuf,3 * (width+1),fp);
  return 0;
}

void gt_eoplist_format_generic(FILE *fp,
                               GtEoplistReader *eoplist_reader,
                               const GtUchar *useq,
                               GtUword ulen,
                               const GtUchar *vseq,
                               GtUword vlen,
                               bool distinguish_mismatch_match)
{
  GtCigarOp co;
  GtUword idx_u = 0, idx_v = 0;
  unsigned int pos = 0;
  GtUchar *topbuf = eoplist_reader->outbuffer, *midbuf = NULL, *lowbuf = NULL;
  const GtUchar matchsymbol = '|',
                mismatchsymbol = ' ',
                gapsymbol = '-';

  gt_assert(eoplist_reader != NULL);
  if (distinguish_mismatch_match)
  {
    gt_eoplist_reader_distinguish_mismatch_match(eoplist_reader);
  }
  topbuf[eoplist_reader->width] = '\n';
  midbuf = topbuf + eoplist_reader->width + 1;
  midbuf[eoplist_reader->width] = '\n';
  lowbuf = midbuf + eoplist_reader->width + 1;
  lowbuf[eoplist_reader->width] = '\n';
  while (gt_eoplist_reader_next_cigar(&co,eoplist_reader))
  {
    switch (co.eoptype)
    {
      GtUword j;

      case GtMatchOp:
      case GtMismatchOp:
        for (j = 0; j < co.iteration && idx_u < ulen && idx_v < vlen; j++)
        {
          topbuf[pos] = useq[idx_u++];
          lowbuf[pos] = vseq[idx_v++];
          if (eoplist_reader->distinguish_mismatch_match)
          {
            midbuf[pos] = co.eoptype == GtMatchOp ? matchsymbol
                                                  : mismatchsymbol;
          } else
          {
            midbuf[pos] = topbuf[pos] == lowbuf[pos] ? matchsymbol
                                                     : mismatchsymbol;
          }
          pos = gt_eoplist_show_advance(pos,eoplist_reader->width,topbuf,fp);
        }
        break;
      case GtDeletionOp:
        for (j = 0; j < co.iteration && idx_u < ulen; j++)
        {
          topbuf[pos] = useq[idx_u++];
          midbuf[pos] = mismatchsymbol;
          lowbuf[pos] = gapsymbol;
          pos = gt_eoplist_show_advance(pos,eoplist_reader->width,topbuf,fp);
        }
        break;
      case GtInsertionOp:
        for (j = 0; j < co.iteration && idx_v < vlen; j++)
        {
          topbuf[pos] = gapsymbol;
          midbuf[pos] = mismatchsymbol;
          lowbuf[pos] = vseq[idx_v++];
          pos = gt_eoplist_show_advance(pos,eoplist_reader->width,topbuf,fp);
        }
        break;
      default:
        gt_assert(false);
    }
  }
  if (pos > 0)
  {
    topbuf[pos] = '\n';
    fwrite(topbuf,sizeof *topbuf,pos+1,fp);
    midbuf[pos] = '\n';
    fwrite(midbuf,sizeof *midbuf,pos+1,fp);
    lowbuf[pos] = '\n';
    fwrite(lowbuf,sizeof *lowbuf,pos+1,fp);
  }
}
