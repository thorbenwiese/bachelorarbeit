#ifndef SUBSTRING_H
#define SUBSTRING_H

#include <stdlib.h>
#include "gt-alloc.h"
#include "gt-defs.h"


/* function to get a pointer to a substring*/
GtUchar *substring(GtUchar *string, GtUword position, GtUword length);

#endif
