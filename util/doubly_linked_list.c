#include <stdio.h>
#include <assert.h>
#include "doubly_linked_list.h"
#include "all.h"

/* Comment in .h file. */
int doubly_linked_list_check_invariant(doubly_linked_list_element* head, doubly_linked_list_element* tail, doubly_linked_list_element* element_in_list)
{
  int error                 = FALSE; // Error flag.
  int found_element_in_list = FALSE; // Flag to check if we have found element_in_list.
  
  if (NULL == head)
    {
      if (NULL != tail)
        {
          fprintf(stderr, "ERROR: head of list is NULL but tail is not NULL.\n");
          error = TRUE;
        }
    }
  else
    {
      if (NULL == tail || NULL != head->prev || NULL != tail->next)
        {
          fprintf(stderr, "ERROR: incorrect linkage of first or last element in list.\n");
          error = TRUE;
        }
      
      doubly_linked_list_element* temp_element    = head;
      doubly_linked_list_element* laggard_element = head;  // Used to check for cycles.
      int                         advance_laggard = FALSE; // Toggle to advance laggard_element every other time through loop.
      
      while (NULL != temp_element)
        {
          if (element_in_list == temp_element)
            {
              found_element_in_list = TRUE;
            }
          
          if (NULL != temp_element->next)
            {
              if (temp_element->next->prev != temp_element)
                {
                  fprintf(stderr, "ERROR: missing reciprocal linkage in list.\n");
                  error = TRUE;
                }
            }
          else
            {
              if (tail != temp_element)
                {
                  fprintf(stderr, "ERROR: following the list from head does not arrive at tail.\n");
                  error = TRUE;
                }
            }
          
          temp_element = temp_element->next;
          
          if (advance_laggard)
            {
              laggard_element = laggard_element->next;
              advance_laggard = FALSE;
            }
          else
            {
              advance_laggard = TRUE;
            }
          
          if (laggard_element == temp_element)
            {
              fprintf(stderr, "ERROR: cycle in list.\n");
              error = TRUE;
            }
        }
    }
  
  if (NULL != element_in_list && !found_element_in_list)
    {
      fprintf(stderr, "ERROR: an element required to be in the list is missing.\n");
      error = TRUE;
    }
  
  return error;
}

/* Comment in .h file. */
int doubly_linked_list_insert_after(doubly_linked_list_element** head, doubly_linked_list_element** tail,
                                    doubly_linked_list_element* prev_element, doubly_linked_list_element* element_to_insert)
{
  int error = FALSE; // Error flag.
  
#if (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)
  if (NULL == head)
    {
      fprintf(stderr, "ERROR: head must not be NULL.\n");
      error = TRUE;
    }
  
  if (NULL == tail)
    {
      fprintf(stderr, "ERROR: tail must not be NULL.\n");
      error = TRUE;
    }
  
  if (NULL == element_to_insert)
    {
      fprintf(stderr, "ERROR: element_to_insert must not be NULL.\n");
      error = TRUE;
    }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)
  
#if (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_INVARIANTS)
  if (!error)
    {
      error = doubly_linked_list_check_invariant(*head, *tail, prev_element);
    }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_INVARIANTS)

  if (!error)
    {
      element_to_insert->prev = prev_element; // If prev_element is NULL element_to_insert->prev should be NULL.

      if (NULL == element_to_insert->prev)
        {
          // Insert at the head of the list.
          element_to_insert->next = *head;
          *head                   = element_to_insert;
        }
      else
        {
          // Insert not at the head of the list.
          element_to_insert->next       = element_to_insert->prev->next;
          element_to_insert->prev->next = element_to_insert;
        }

      if (NULL == element_to_insert->next)
        {
#if (DEBUG_LEVEL & DEBUG_LEVEL_INTERNAL_ASSERTIONS)
          assert(*tail == prev_element); // tail should be prev_element or NULL if prev_element is NULL.
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_INTERNAL_ASSERTIONS)

          // Insert at the tail of the list.
          *tail = element_to_insert;
        }
      else
        {
          // Insert not at the tail of the list.
          element_to_insert->next->prev = element_to_insert;
        }
    }
  
  return error;
}

/* Comment in .h file. */
int doubly_linked_list_insert_before(doubly_linked_list_element** head, doubly_linked_list_element** tail,
                                     doubly_linked_list_element* next_element, doubly_linked_list_element* element_to_insert)
{
  doubly_linked_list_element* prev_element = NULL;
  
  if (NULL != next_element)
    {
      prev_element = next_element->prev; // If next_element is at the head of the list prev_element will get set to NULL
                                         // and doubly_linked_list_insert_after will insert at the head of the list.
    }
  else if (NULL != tail) // Protect against seg fault on *tail.  If false, doubly_linked_list_insert_after will throw the error.
    {
      prev_element = *tail; // Insert after the tail of the list.
    }
  
  return doubly_linked_list_insert_after(head, tail, prev_element, element_to_insert);
}

/* Comment in .h file. */
int doubly_linked_list_remove_element(doubly_linked_list_element** head, doubly_linked_list_element** tail, doubly_linked_list_element* element_to_remove)
{
  int error = FALSE; // Error flag.
  
#if (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)
  if (NULL == head)
    {
      fprintf(stderr, "ERROR: head must not be NULL.\n");
      error = TRUE;
    }
  
  if (NULL == tail)
    {
      fprintf(stderr, "ERROR: tail must not be NULL.\n");
      error = TRUE;
    }
  
  if (NULL == element_to_remove)
    {
      fprintf(stderr, "ERROR: element_to_remove must not be NULL.\n");
      error = TRUE;
    }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE)
  
#if (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_INVARIANTS)
  if (!error)
    {
      error = doubly_linked_list_check_invariant(*head, *tail, element_to_remove);
    }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_PUBLIC_FUNCTIONS_INVARIANTS)

  if (!error)
    {
      if (NULL == element_to_remove->prev)
        {
#if (DEBUG_LEVEL & DEBUG_LEVEL_INTERNAL_ASSERTIONS)
          assert(*head == element_to_remove); // Remove from the head of the list.
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_INTERNAL_ASSERTIONS)

          *head = element_to_remove->next;
        }
      else
        {
          // Remove not from the head of the list.
          element_to_remove->prev->next = element_to_remove->next;
        }

      if (NULL == element_to_remove->next)
        {
#if (DEBUG_LEVEL & DEBUG_LEVEL_INTERNAL_ASSERTIONS)
          assert(*tail == element_to_remove); // Remove from the tail of the list.
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_INTERNAL_ASSERTIONS)

          *tail = element_to_remove->prev;
        }
      else
        {
          // Remove not from the tail of the list.
          element_to_remove->next->prev = element_to_remove->prev;
        }
    }
  
  return error;
}
