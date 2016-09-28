#include "substring.h"

GtUchar *substring(GtUchar *string, GtUword position, GtUword length) 
{
   GtUchar *pointer;
   GtUword i;
 
   pointer = gt_malloc(length++);
 
   for (i = 0 ; i < position - 1; i++) 
      string++; 
 
   for (i = 0 ; i < length ; i++)
   {
      *(pointer + i) = *string;      
      string++;   
   }
 
   *(pointer + i) = '\0';
 
   return pointer;
}
