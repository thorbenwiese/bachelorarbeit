#ifndef FRONT_WITH_TRACE_H
#define FRONT_WITH_TRACE_H
#include <stdbool.h>
#include "gt-defs.h"
#include "eoplist.h"

/* For computing the front and keeping track of the relevant traceback
   information, we use the following opaque type */

typedef struct FrontEdistTrace FrontEdistTrace;

/* Here is the corresponding constructor */

FrontEdistTrace *front_edist_trace_new(void);

/* and here is the corresponding destructor */

void front_edist_trace_delete(FrontEdistTrace *fet);

/* The following function returns the unit edit distance of the given
   sequences <useq> and <vseq> of length <ulen> and <vlen>, respectively. In
   the object pointed to by <fet>, the bracktracing information is
   stored. Note that in case of performing many alignments, one only needs
   to create a single FrontEdistTrace-Object, as this is internally
   reset in each call of the following function. */

GtUword front_edist_trace_distance(FrontEdistTrace *fet,
                                   const GtUchar *useq,
                                   GtUword ulen,
                                   const GtUchar *vseq,
                                   GtUword vlen);

/* In our test framework we need a function with is like the prevous
   but requires a void pointer. So we added the following function. */

GtUword front_edist_trace_distance_void(void *info,
                                        const GtUchar *useq,
                                        GtUword ulen,
                                        const GtUchar *vseq,
                                        GtUword vlen);

/* The following function returns the unit edit distance of the given
   sequences <useq> and <vseq> of length <ulen> and <vlen>, respectively. In
   the object pointed to by <fet>, the bracktracing information is
   stored. Moreover, the backtracing is performed and the constructed
   list of edit operations is appended to the eoplist. So the method also
   works for the case that the eoplist already contains
   edit operations. Note that in case of performing many alignments, one
   only needs to create a single FrontEdistTrace-Object, as this is
   internally reset in each call of the following function. */

GtUword front_edist_trace_eoplist(GtEoplist *eoplist,
                                  FrontEdistTrace *fet,
                                  const GtUchar *useq,
                                  GtUword ulen,
                                  const GtUchar *vseq,
                                  GtUword vlen,
                                  bool testeoplist);

typedef struct
{
  FrontEdistTrace *fet;
  GtEoplist *eoplist;
} FrontEdistTracewitheoplist;

/* In our test framework we need a function with is like the previous
   but requires a void pointer, referring to a structure of type
   <FrontEdistTracewitheoplist>, which itself contains references to the
   first two arguments of front_edist_trace_eoplist. */

GtUword front_edist_trace_eoplist_void(void *info,
                                       const GtUchar *useq,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vlen);

GtUword front_edist_segments_void(void *info,
                                  const GtUchar *useq,
                                  GtUword ulen,
                                  const GtUchar *vseq,
                                  GtUword vlen);
#endif
