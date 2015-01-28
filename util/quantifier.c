#include <stdlib.h>
#include "quantifier.h"
#include "all.h"

/* Comment in .h file. */
int i_test_equals(int a, void* b)
{
  return (NULL == b) ? FALSE : (a == *(int*)b);
}

/* Comment in .h file. */
int i_exists(int* array, int size, int(*test)(int, void*), void* test_params)
{
  int ii; // Loop counter.
  int result = FALSE;

  if (NULL != array && NULL != test)
    {
      for(ii = 1; !result && ii <= size; ii++)
        {
          result = test(array[ii], test_params);
        }
    }

  return result;
}

/* Comment in .h file. */
int i_forall(int* array, int size, int(*test)(int, void*), void* test_params)
{
  int ii; // Loop counter.
  int result = TRUE;

  if (NULL != array && NULL != test)
    {
      for(ii = 1; result && ii <= size; ii++)
        {
          result = test(array[ii], test_params);
        }
    }

  return result;
}
