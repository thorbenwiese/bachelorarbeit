#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include "gt-defs.h"
#include "gt-alloc.h"
#include "tracer.h"
#include "front-with-trace.h"

typedef unsigned int Rowvaluetype;
typedef uint8_t Backreferencetype;

typedef struct
{
  Rowvaluetype row,
               localmatch_count;
  Backreferencetype backreference;
} Frontvalue;

#define FRONT_DIAGONAL(FRONTPTR) (GtWord) ((FRONTPTR) - midfront)

static void inline front_prune_add_matches(Frontvalue *midfront,
                                           Frontvalue *fv,
                                           const GtUchar *useq,
                                           GtUword ulen,
                                           const GtUchar *vseq,
                                           GtUword vlen)
{
  GtUword upos, vpos;

  for (upos = fv->row, vpos = fv->row + FRONT_DIAGONAL(fv);
       upos < ulen && vpos < vlen &&
       useq[upos] == vseq[vpos];
       upos++, vpos++)
       /* Nothing */;
  fv->localmatch_count = upos - fv->row;
  fv->row = upos;
}

static void front_next_inplace(Frontvalue *midfront,
                               Frontvalue *lowfront,
                               Frontvalue *highfront,
                               const GtUchar *useq,
                               GtUword ulen,
                               const GtUchar *vseq,
                               GtUword vlen)
{
  Frontvalue bestfront, insertion_value, replacement_value, *frontptr;

  insertion_value = *lowfront; /* from previous diag -(d-1) => -d => DELETION */
  bestfront = insertion_value;
  bestfront.row++;
  *lowfront = bestfront;
  lowfront->backreference = FT_EOP_DELETION;
  front_prune_add_matches(midfront,lowfront,useq,ulen,vseq,vlen);

  replacement_value = *(lowfront+1);
  if (bestfront.row < replacement_value.row + 1)
  {
    bestfront = replacement_value;
    bestfront.backreference = FT_EOP_DELETION;
    bestfront.row++;
  } else
  {
    bestfront.backreference = FT_EOP_MISMATCH;
    if (bestfront.row == replacement_value.row + 1)
    {
      bestfront.backreference |= FT_EOP_DELETION;
    }
  }
  *(lowfront+1) = bestfront;
  front_prune_add_matches(midfront,lowfront + 1,useq,ulen,vseq,vlen);
  for (frontptr = lowfront+2; frontptr <= highfront; frontptr++)
  {
    bestfront = insertion_value;
    bestfront.backreference = FT_EOP_INSERTION;
    if (frontptr <= highfront - 1)
    {
      if (bestfront.row < replacement_value.row + 1)
      {
        bestfront = replacement_value;
        bestfront.backreference = FT_EOP_MISMATCH;
        bestfront.row++;
      } else
      {
        if (bestfront.row == replacement_value.row + 1)
        {
          bestfront.backreference |= FT_EOP_MISMATCH;
        }
      }
    }
    if (frontptr <= highfront - 2)
    {
      if (bestfront.row < frontptr->row + 1)
      {
        bestfront = *frontptr;
        bestfront.backreference = FT_EOP_DELETION;
        bestfront.row++;
      } else
      {
        if (bestfront.row == frontptr->row + 1)
        {
          bestfront.backreference |= FT_EOP_DELETION;
        }
      }
    }
    if (frontptr < highfront)
    {
      insertion_value = replacement_value;
      replacement_value = *frontptr;
    }
    *frontptr = bestfront;
    front_prune_add_matches(midfront,frontptr,useq,ulen,vseq,vlen);
  }
}

static void front_second_inplace(Frontvalue *midfront,
                                 Frontvalue *lowfront,
                                 const GtUchar *useq,
                                 GtUword ulen,
                                 const GtUchar *vseq,
                                 GtUword vlen)
{
  *(lowfront+1) = *(lowfront+2) = *lowfront;
  lowfront->row++;
  lowfront->backreference = FT_EOP_DELETION;
  front_prune_add_matches(midfront,lowfront,useq,ulen,vseq,vlen);

  (lowfront+1)->row++;
  (lowfront+1)->backreference = FT_EOP_MISMATCH;
  front_prune_add_matches(midfront,lowfront + 1,useq,ulen,vseq,vlen);

  (lowfront+2)->backreference = FT_EOP_INSERTION;
  front_prune_add_matches(midfront,lowfront + 2,useq,ulen,vseq,vlen);
}

struct FrontEdistTrace
{
  Frontvalue *spaceFrontvalue;
  GtUword allocatedFrontvalue;
  GtFrontTrace *front_trace;
};

FrontEdistTrace *front_edist_trace_new(void)
{
  FrontEdistTrace *fet = gt_malloc(sizeof *fet);

  gt_assert(fet != NULL);
  fet->spaceFrontvalue = NULL;
  fet->allocatedFrontvalue = 0;
  fet->front_trace = front_trace_new();
  return fet;
}

void front_edist_trace_delete(FrontEdistTrace *fet)
{
  if (fet != NULL)
  {
    gt_free(fet->spaceFrontvalue);
    front_trace_delete(fet->front_trace);
    gt_free(fet);
  }
}

static void front_trace_add_gen(GtFrontTrace *front_trace,
                                const Frontvalue *lowfront,
                                const Frontvalue *highfront)
{
  const Frontvalue *fv;

  for (fv = lowfront; fv <= highfront; fv++)
  {
    front_trace_add_trace(front_trace,fv->backreference,fv->localmatch_count);
  }
}

GtUword front_edist_trace_distance(FrontEdistTrace *fet,
                                   const GtUchar *useq,
                                   GtUword ulen,
                                   const GtUchar *vseq,
                                   GtUword vlen)
{
  const GtUword sumseqlength = ulen + vlen;
  GtUword distance;

  front_trace_reset(fet->front_trace);
  for (distance = 0; distance <= sumseqlength; distance++)
  {
    Frontvalue *basefront;
    if (2 * distance >= fet->allocatedFrontvalue)
    {
      fet->allocatedFrontvalue = fet->allocatedFrontvalue * 1.2 + 32;
      fet->spaceFrontvalue = gt_realloc(fet->spaceFrontvalue,
                                        sizeof *fet->spaceFrontvalue *
                                        fet->allocatedFrontvalue);
      gt_assert(fet->spaceFrontvalue != NULL);
    }
    basefront = fet->spaceFrontvalue;
    if (distance == 0)
    {
      basefront->row = 0;
      basefront->backreference = 0; /* No back reference */
      front_prune_add_matches(basefront,basefront,useq,ulen,vseq,vlen);
    } else
    {
      if (distance == 1)
      {
        front_second_inplace(basefront + distance,basefront,
                             useq,ulen,vseq,vlen);
      } else
      {
        front_next_inplace(basefront + distance,basefront,
                           basefront + 2 * distance,
                           useq,ulen,vseq,vlen);
      }
    }
    front_trace_add_gen(fet->front_trace,basefront,basefront + 2 * distance);
    if ((vlen > ulen && vlen - ulen <= distance) ||
        (vlen <= ulen && ulen - vlen <= distance))
    {
      if (basefront[distance + vlen - ulen].row == ulen)
      {
        break;
      }
    }
  }
  gt_assert(distance <= sumseqlength);
  return distance;
}

GtUword front_edist_trace_distance_void(void *info,
                                        const GtUchar *useq,
                                        GtUword ulen,
                                        const GtUchar *vseq,
                                        GtUword vlen)
{
  return front_edist_trace_distance((FrontEdistTrace *) info,
                                    useq,
                                    ulen,
                                    vseq,
                                    vlen);
}

GtUword front_edist_trace_eoplist(GtEoplist *eoplist,
                                  FrontEdistTrace *fet,
                                  const GtUchar *useq,
                                  GtUword ulen,
                                  const GtUchar *vseq,
                                  GtUword vlen,
                                  bool testeoplist)
{
  GtUword distance = front_edist_trace_distance(fet, useq, ulen, vseq, vlen);

  gt_assert(eoplist != NULL);
  if (testeoplist)
  {
    front_trace2eoplist_directed(NULL,
                                 fet->front_trace,
                                 distance,
                                 useq,
                                 ulen,
                                 vseq,
                                 vlen);
    gt_eoplist_reset(eoplist);
  }
  front_trace2eoplist_directed(eoplist,
                               fet->front_trace,
                               distance,
                               useq,
                               ulen,
                               vseq,
                               vlen);
  if (testeoplist)
  {
    gt_eoplist_reader_verify(eoplist,
                             useq,
                             ulen,
                             vseq,
                             vlen,
                             distance,
                             true);
    gt_eoplist_reader_verify(eoplist,
                             useq,
                             ulen,
                             vseq,
                             vlen,
                             distance,
                             false);
  }
  return distance;
}

GtUword front_edist_trace_eoplist_void(void *info,
                                       const GtUchar *useq,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vlen)
{
  FrontEdistTracewitheoplist *fet_eoplist = (FrontEdistTracewitheoplist *) info;

  return front_edist_trace_eoplist(fet_eoplist->eoplist,
                                   fet_eoplist->fet,
                                   useq,
                                   ulen,
                                   vseq,
                                   vlen,
                                   true);
}

void front_edist_segments(GtEoplist *eoplist,
                          FrontEdistTrace *fet,
                          const GtUchar *useq,
                          GtUword ulen,
                          const GtUchar *vseq,
                          GtUword vlen,
                          GtUword delta)
{
  GtUword distance = front_edist_trace_distance(fet, useq, ulen, vseq, vlen),
          unit_cost, offset_u = 0, offset_v = 0, sum_segment_dist = 0;
  GtEoplistReader *eoplist_reader;
  GtEoplistSegment segment;
  GtEoplist *eoplist_concat = gt_eoplist_new(), *eoplist2;
  char *cigar_string, *cigar_string2;

  gt_assert(eoplist != NULL);
  gt_eoplist_reset(eoplist);
  front_trace2eoplist_directed(eoplist,
                               fet->front_trace,
                               distance,
                               useq,
                               ulen,
                               vseq,
                               vlen);
  eoplist_reader = gt_eoplist_reader_new(eoplist);
  while (gt_eoplist_reader_next_segment(&segment,eoplist_reader,delta))
  {
    FrontEdistTrace *fet_segment = front_edist_trace_new();
    GtUword this_distance = front_edist_trace_distance(fet_segment,
                                                       useq + offset_u,
                                                       segment.aligned_u,
                                                       vseq + offset_v,
                                                       segment.aligned_v);
    sum_segment_dist += this_distance;
    front_trace2eoplist_directed(eoplist_concat,
                                 fet_segment->front_trace,
                                 this_distance,
                                 useq + offset_u,
                                 segment.aligned_u,
                                 vseq + offset_v,
                                 segment.aligned_v);
    offset_u += segment.aligned_u;
    offset_v += segment.aligned_v;
    front_edist_trace_delete(fet_segment);
  }
  gt_assert(offset_u == ulen);
  gt_assert(offset_v == vlen);
  gt_eoplist_reader_delete(eoplist_reader);
  unit_cost = gt_eoplist_unit_cost(eoplist);
  gt_assert(sum_segment_dist <= unit_cost);
  gt_assert(sum_segment_dist == gt_eoplist_unit_cost(eoplist_concat));
  gt_eoplist_reader_verify(eoplist_concat,
                           useq,
                           ulen,
                           vseq,
                           vlen,
                           sum_segment_dist,
                           true);
  cigar_string = gt_eoplist2cigar_string(eoplist_concat,true);
  eoplist2 = gt_eoplist_new_from_cigar(cigar_string,strlen(cigar_string));
  cigar_string2 = gt_eoplist2cigar_string(eoplist2,true);
  gt_assert(strcmp(cigar_string,cigar_string2) == 0);
  gt_free(cigar_string);
  gt_free(cigar_string2);
  gt_eoplist_delete(eoplist2);
  gt_eoplist_delete(eoplist_concat);
}

GtUword front_edist_segments_void(void *info,
                                  const GtUchar *useq,
                                  GtUword ulen,
                                  const GtUchar *vseq,
                                  GtUword vlen)
{
  FrontEdistTracewitheoplist *fet_eoplist = (FrontEdistTracewitheoplist *) info;
  const GtUword delta = (GtUword) 40;

  front_edist_segments(fet_eoplist->eoplist,
                       fet_eoplist->fet,
                       useq,
                       ulen,
                       vseq,
                       vlen,
                       delta);
  return gt_eoplist_unit_cost(fet_eoplist->eoplist);
}
