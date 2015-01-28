#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "memfunc.h"
#include "all.h"

#define THREAD_SAFE // Leave this defined to have the code use mutexes to be thread safe.

#ifdef THREAD_SAFE
#include <pthread.h>

static pthread_mutex_t allocated_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif // THREAD_SAFE

// A singly linked list that records what memory was allocated.
typedef struct allocation_record allocation_record;
struct allocation_record
{
  void*              location; // The pointer to the allocated memory.
  int                bytes;    // The number of bytes allocated.
  allocation_record* next;     // The next record in the list.
};

static allocation_record* allocated = NULL;

/* Comment in .h file. */
int v_alloc(void** ptr, int bytes)
{
  int                error = FALSE; // Error flag.
  allocation_record* record;        // For adding a record of the allocation to allocated.
  
#if (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)
  if (NULL == ptr)
    {
      fprintf(stderr, "ERROR: ptr must not be NULL.");
      error = TRUE;
    }
  else
    {
      *ptr = NULL;
    }
  
  if (0 >= bytes)
    {
      fprintf(stderr, "ERROR: bytes must be greater than zero.");
      error = TRUE;
    }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)

  if (!error)
    {
      // Allocate memory.
      *ptr = malloc(bytes);

#if (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
      if (NULL == *ptr)
        {
          fprintf(stderr, "ERROR: Failed to allocate memory.");
          error = TRUE;
        }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
    }

  if (!error)
    {
      // Initialize memory to zeros.
      memset(*ptr, 0, bytes);
      
      // Record allocation to make sure it gets freed.
      record = malloc(sizeof(allocation_record));

#if (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
      if (NULL == record)
        {
          free(*ptr);
          *ptr = NULL;
          fprintf(stderr, "ERROR: Failed to allocate memory.");
          error = TRUE;
        }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
    }
  
  if (!error)
    {
      record->location = *ptr;
      record->bytes = bytes;

#ifdef THREAD_SAFE
#if (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
      error =
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
          pthread_mutex_lock(&allocated_mutex);
      
#if (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
      if (error)
        {
          free(record);
          free(*ptr);
          *ptr = NULL;
          fprintf(stderr, "ERROR: Failed to lock mutex.");
        }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
#endif // THREAD_SAFE
    }
  
  if (!error)
    {
      record->next = allocated;
      allocated = record;

#ifdef THREAD_SAFE
#if (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
      error =
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
          pthread_mutex_unlock(&allocated_mutex);
      
#if (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
      if (error)
        {
          allocated = record->next;
          free(record);
          free(*ptr);
          *ptr = NULL;
          fprintf(stderr, "ERROR: Failed to unlock mutex.");
        }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
#endif // THREAD_SAFE
    }

  return error;
}

/* Comment in .h file. */
int v_dealloc(void** ptr, int bytes)
{
  int error = FALSE; // Error flag.

#if (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)
  if (NULL == ptr || NULL == *ptr)
    {
      fprintf(stderr, "ERROR: ptr must not be NULL.");
      error = TRUE;
    }
  
  if (0 >= bytes)
    {
      fprintf(stderr, "ERROR: bytes must be greater than zero.");
      error = TRUE;
    }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)

  if (NULL != ptr && NULL != *ptr)
    {
#ifdef THREAD_SAFE
#if (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
      int error2 =
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
          pthread_mutex_lock(&allocated_mutex);
#if (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
      if (error2)
        {
          error = TRUE;
          fprintf(stderr, "ERROR: Failed to lock mutex.");
        }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
#endif // THREAD_SAFE
      
      // Find the record of the allocation.
      allocation_record** record = &allocated;

      /* FIXME this search is really slow.  For now we are commenting it out.  The code will just take the first record.  It will still verify that there are
       * the same number of deallocations as allocations, but it won't verify that the pointers and number of bytes match.
      while (NULL != *record && (*record)->location != *ptr)
        {
          record = &(*record)->next;
        }
      */

      // Free memory
      free(*ptr);
      *ptr = NULL;

      // Record that the allocation has been freed.
      if (NULL != *record)
        {
          allocation_record* temp_record = *record;

#if (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)
          /* FIXME Since we don't have the correct record (see the FIXME above) don't check that bytes is correct.
          if (bytes != (*record)->bytes)
            {
              fprintf(stderr, "WARNING: Deallocating memory with a different value for bytes than was used for allocation.\n");
            }
          */
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)
          
          *record = (*record)->next;
          free(temp_record);
        }
#if (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)
      else
        {
          fprintf(stderr, "WARNING: Deallocating memory with no record of the location having been allocated.\n");
        }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)
      
#ifdef THREAD_SAFE
#if (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
      error2 =
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
          pthread_mutex_unlock(&allocated_mutex);
#if (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
      if (error2)
        {
          error = TRUE;
          fprintf(stderr, "ERROR: Failed to unlock mutex.");
        }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
#endif // THREAD_SAFE
    }
  
  return error;
}

/* Comment in .h file. */
int i_alloc(int** array, int size)
{
  return v_alloc((void**)array, (size + 1) * sizeof(int));
}

/* Comment in .h file. */
int i_dealloc(int** array, int size)
{
  return v_dealloc((void**)array, (size + 1) * sizeof(int));
}

/* Comment in .h file. */
int d_alloc(double** array, int size)
{
  return v_alloc((void**)array, (size + 1) * sizeof(double));
}

/* Comment in .h file. */
int d_dealloc(double** array, int size)
{
  return v_dealloc((void**)array, (size + 1) * sizeof(double));
}

/* Comment in .h file. */
int vtwo_alloc(void*** array, int rows, int bytes)
{
  int error = FALSE; // Error flag.
  int ii;            // Loop Counter.

#if (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)
  if (NULL == array)
    {
      fprintf(stderr, "ERROR: array must not be NULL.");
      error = TRUE;
    }
  else
    {
      *array = NULL;
    }
  
  if (0 >= rows)
    {
      fprintf(stderr, "ERROR: rows must be greater than zero.");
      error = TRUE;
    }
  
  if (0 >= bytes)
    {
      fprintf(stderr, "ERROR: bytes must be greater than zero.");
      error = TRUE;
    }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)

  if (!error)
    {
      error = v_alloc((void**)array, rows * sizeof(void*));
    }

  if (!error)
    {
      // Allocate all of the memory for the rows.  Assign (*array)[0] to point to the zeroth row.
      error = v_alloc(*array, rows * bytes);

      if (error)
        {
          v_dealloc((void**)array, rows * sizeof(void*));
        }
    }

  if (!error)
    {
      for (ii = 1; ii < rows; ii++)
        {
          // Assign each row pointer.
          (*array)[ii] = (*array)[0] + (ii * bytes);
        }
    }

  return (error);
}

/* Comment in .h file. */
int vtwo_dealloc(void*** array, int rows, int bytes)
{
  int error = FALSE; // Error flag.

#if (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)
  if (NULL == array || NULL == *array || NULL == (*array)[0])
    {
      fprintf(stderr, "ERROR: array must not be NULL.");
      error = TRUE;
    }
  
  if (0 >= rows)
    {
      fprintf(stderr, "ERROR: rows must be greater than zero.");
      error = TRUE;
    }
  
  if (0 >= bytes)
    {
      fprintf(stderr, "ERROR: bytes must be greater than zero.");
      error = TRUE;
    }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)

  if (NULL != array && NULL != *array)
    {
      if (NULL != (*array)[0])
        {
          if (v_dealloc(*array, rows * bytes))
            {
              error = TRUE;
            }
        }

      if (v_dealloc((void**)array, rows * sizeof(void*)))
        {
          error = TRUE;
        }
    }
  
  return error;
}

/* Comment in .h file. */
int itwo_alloc(int*** array, int rows, int cols)
{
  return vtwo_alloc((void***)array, rows + 1, (cols + 1) * sizeof(int));
}

/* Comment in .h file. */
int itwo_dealloc(int*** array, int rows, int cols)
{
  return vtwo_dealloc((void***)array, rows + 1, (cols + 1) * sizeof(int));
}

/* Comment in .h file. */
int dtwo_alloc(double*** array, int rows, int cols)
{
  return vtwo_alloc((void***)array, rows + 1, (cols + 1) * sizeof(double));
}

/* Comment in .h file. */
int dtwo_dealloc(double*** array, int rows, int cols)
{
  return vtwo_dealloc((void***)array, rows + 1, (cols + 1) * sizeof(double));
}

/* Comment in .h file. */
int check_all_memory_freed(int cleanup)
{
  int error = FALSE; // Error flag.

#ifdef THREAD_SAFE
#if (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
  error =
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
      pthread_mutex_lock(&allocated_mutex);

#if (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
  if (error)
    {
      fprintf(stderr, "ERROR: Failed to lock mutex.");
    }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
#endif // THREAD_SAFE
  
  if (!error)
    {
      allocation_record* record = allocated;

      while (NULL != record)
        {
          fprintf(stderr, "WARNING: %d bytes of unfreed memory\n", record->bytes);
          record = record->next;

          if (cleanup)
            {
              free(allocated->location);
              free(allocated);
              allocated = record;
            }
        }

#if (DEBUG_LEVEL & DEBUG_LEVEL_INTERNAL_ASSERTIONS)
      assert(!cleanup || NULL == allocated); // Cleanup implies allocated is NULL.
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_INTERNAL_ASSERTIONS)

#ifdef THREAD_SAFE
#if (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
      error =
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
          pthread_mutex_unlock(&allocated_mutex);

#if (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
      if (error)
        {
          fprintf(stderr, "ERROR: Failed to unlock mutex.");
        }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_LIBRARY_ERRORS)
#endif // THREAD_SAFE
    }
  
  return error;
}
