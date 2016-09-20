#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include "gt-defs.h"
#include "gt-alloc.h"
#include "tracer.h"

#define BACKTRACEBITS 3

typedef struct
{
  uint32_t bits:BACKTRACEBITS,     /* combination of FT_EOP_MISMATCH
                                                     FT_EOP_INSERTION
                                                     FT_EOP_DELETION */
           lcs:(32-BACKTRACEBITS); /* longest common suffix */
} GtBackreftable;

struct GtFrontTrace
{
  GtBackreftable *backref_table;
  GtUword maxlcs,
          backref_nextfree,
          backref_allocated;
};

GtFrontTrace *front_trace_new(void)
{
  GtFrontTrace *front_trace = gt_malloc(sizeof *front_trace);

  gt_assert(front_trace != NULL);
  front_trace->maxlcs = (1 << (32-BACKTRACEBITS)) - 1;
  front_trace->backref_table = NULL;
  front_trace->backref_nextfree = 0;
  front_trace->backref_allocated = 0;
  return front_trace;
}

void front_trace_delete(GtFrontTrace *front_trace)
{
  if (front_trace != NULL)
  {
    gt_free(front_trace->backref_table);
    gt_free(front_trace);
  }
}

void front_trace_reset(GtFrontTrace *front_trace)
{
  gt_assert (front_trace != NULL);
  front_trace->backref_nextfree = 0;
}

void front_trace_add_trace(GtFrontTrace *front_trace,
                           uint8_t backreference,
                           uint32_t localmatch_count)
{
  gt_assert (front_trace != NULL);
  if (front_trace->backref_nextfree >= front_trace->backref_allocated)
  {
    front_trace->backref_allocated
      = front_trace->backref_allocated * 1.2 + 128UL;
    front_trace->backref_table
      = gt_realloc(front_trace->backref_table,sizeof *front_trace->backref_table
                                              * front_trace->backref_allocated);
    gt_assert(front_trace->backref_table != NULL);
  }
  gt_assert(front_trace->backref_nextfree < front_trace->backref_allocated);
  front_trace->backref_table[front_trace->backref_nextfree].bits
    = backreference;
  gt_assert(localmatch_count <= front_trace->maxlcs);
  front_trace->backref_table[front_trace->backref_nextfree++].lcs
    = localmatch_count;
}

static void gt_check_diagonal_run(const GtUchar *useq,
                                  const GtUchar *vseq,
                                  GtWord diagonal,
                                  unsigned int firstrow,
                                  unsigned int nextrow)
{
  GtUword idx;

  gt_assert(useq != NULL && vseq != NULL && firstrow <= nextrow);
  for (idx = firstrow; idx < nextrow; idx++)
  {
    gt_assert (useq[idx] == vseq[idx+diagonal]);
  }
}

void front_trace2eoplist_directed(GtEoplist *eoplist,
                                  const GtFrontTrace *front_trace,
                                  GtUword distance,
                                  const GtUchar *useq,
                                  GtUword ulen,
                                  const GtUchar *vseq,
                                  GtUword vlen)
{
  const GtBackreftable *basefront, *current;
  GtUword firstindex;
  GtWord diagonal = (GtWord) vlen - (GtWord) ulen;
  uint32_t row;
  uint8_t preferred_eop = FT_EOP_MISMATCH;

  gt_assert(front_trace != NULL &&
         front_trace->backref_nextfree >= 2 * distance + 1);
  basefront = front_trace->backref_table + front_trace->backref_nextfree
                                       - (2 * distance + 1);
  current = basefront + distance + diagonal;
  firstindex = gt_eoplist_length(eoplist);
  row = ulen;
  while (distance > 0)
  {
    GtUword nextrowadd;

    if (eoplist != NULL)
    {
      if (current->lcs > 0)
      {
        gt_eoplist_match_add(eoplist,current->lcs);
      }
    } else
    {
      gt_check_diagonal_run(useq, vseq, diagonal, row - current->lcs, row);
    }
    if (current->bits & preferred_eop)
    {
      if (preferred_eop == FT_EOP_MISMATCH)
      {
        nextrowadd = 1;
      } else
      {
        if (preferred_eop == FT_EOP_INSERTION)
        {
          gt_assert(-(GtWord) ulen < diagonal);
          diagonal--;
          nextrowadd = 0;
        } else
        {
          gt_assert(preferred_eop == FT_EOP_DELETION);
          gt_assert(diagonal < (GtWord) vlen);
          diagonal++;
          nextrowadd = 1;
        }
      }
    } else
    {
      if (current->bits & FT_EOP_MISMATCH)
      {
        preferred_eop = FT_EOP_MISMATCH;
        nextrowadd = 1;
      } else
      {
        if (current->bits & FT_EOP_INSERTION)
        {
          gt_assert(-(GtWord) ulen < diagonal);
          diagonal--;
          preferred_eop = FT_EOP_INSERTION;
          nextrowadd = 0;
        } else
        {
          gt_assert(current->bits & FT_EOP_DELETION);
          gt_assert(diagonal < (GtWord) vlen);
          diagonal++;
          preferred_eop = FT_EOP_DELETION;
          nextrowadd = 1;
        }
      }
    }
    if (eoplist != NULL)
    {
      if (preferred_eop == FT_EOP_DELETION)
      {
        gt_eoplist_deletion_add(eoplist);
      } else
      {
        if (preferred_eop == FT_EOP_INSERTION)
        {
          gt_eoplist_insertion_add(eoplist);
        } else
        {
          gt_eoplist_mismatch_add(eoplist);
        }
      }
    }
    distance--;
    basefront -= (2 * distance + 1);
    gt_assert(basefront >= front_trace->backref_table);
    gt_assert(row >= current->lcs + nextrowadd);
    row -= current->lcs + nextrowadd;
    current = basefront + distance + diagonal;
  }
  gt_assert(basefront == front_trace->backref_table && current->bits == 0);
  if (eoplist != NULL)
  {
    if (current->lcs > 0)
    {
      gt_eoplist_match_add(eoplist,current->lcs);
    }
    gt_eoplist_reverse_end(eoplist,firstindex);
  }
}
