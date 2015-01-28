#ifndef DOUBLY_LINKED_LIST_H
#define DOUBLY_LINKED_LIST_H

/* A doubly_linked_list_element struct represents an abstract list element.
 * To use this code, define your list element like doubly_linked_list_element
 * except with a data payload after next, like this:
 * 
 * typedef struct my_list_element my_list_element;
 * struct my_list_element
 * {
 *   my_list_element* prev;
 *   my_list_element* next;
 *   int              data1;
 *   double           data2;
 *   etc.
 * }
 * 
 * prev and next must come first in the definition to have the same location
 * in the struct as in a doubly_linked_list_element struct.  Then call the
 * doubly_linked_list functions and cast your pointers to
 * doubly_linked_list_element.
 */
typedef struct doubly_linked_list_element doubly_linked_list_element;
struct doubly_linked_list_element
{
  doubly_linked_list_element* prev; // The previous element in the list or NULL if this is the first element.
  doubly_linked_list_element* next; // The next     element in the list or NULL if this is the last  element.
                                    // Data payload goes here.
};

/* Return TRUE if the list invariant is violated, FALSE otherwise.
 * Also print an error message indicating the details of what is wrong.
 *
 * Parameters:
 *
 * head            - A pointer to the head of the list.
 * tail            - A pointer to the tail of the list.
 * element_in_list - An element that must be in the list.  If NULL this check
 *                   is ignored.
 */
int doubly_linked_list_check_invariant(doubly_linked_list_element* head, doubly_linked_list_element* tail, doubly_linked_list_element* element_in_list);

/* Add element_to_insert in the list after prev_element, or at the head of the
 * list if prev_element is NULL.
 * Return TRUE if there is an error, FALSE otherwise.
 * If there is an error nothing is inserted into the list.
 * 
 * Parameters:
 *
 * head              - A pointer to the head of the list passed by reference.
 *                     May be updated if element_to_insert is inserted at the
 *                     head of the list.
 * tail              - A pointer to the tail of the list passed by reference.
 *                     May be updated if element_to_insert is inserted at the
 *                     tail of the list.
 * prev_element      - The element to insert after, or NULL to insert at the
 *                     head of the list
 * element_to_insert - The element to insert.
 */
int doubly_linked_list_insert_after(doubly_linked_list_element** head, doubly_linked_list_element** tail,
                                    doubly_linked_list_element* prev_element, doubly_linked_list_element* element_to_insert);

/* Add element_to_insert in the list before next_element, or at the tail of the
 * list if next_element is NULL.
 * Return TRUE if there is an error, FALSE otherwise.
 * If there is an error nothing is inserted into the list.
 * 
 * Parameters:
 *
 * head              - A pointer to the head of the list passed by reference.
 *                     May be updated if element_to_insert is inserted at the
 *                     head of the list.
 * tail              - A pointer to the tail of the list passed by reference.
 *                     May be updated if element_to_insert is inserted at the
 *                     tail of the list.
 * prev_element      - The element to insert before, or NULL to insert at the
 *                     tail of the list
 * element_to_insert - The element to insert.
 */
int doubly_linked_list_insert_before(doubly_linked_list_element** head, doubly_linked_list_element** tail,
                                     doubly_linked_list_element* next_element, doubly_linked_list_element* element_to_insert);

/* Remove element_to_remove from the list.
 * Return TRUE if there is an error, FALSE otherwise.
 * If there is an error nothing is removed from the list.
 * 
 * Parameters:
 *
 * head              - A pointer to the head of the list passed by reference.
 *                     May be updated if element_to_remove is removed from the
 *                     head of the list.
 * tail              - A pointer to the tail of the list passed by reference.
 *                     May be updated if element_to_remove is removed from the
 *                     tail of the list.
 * element_to_remove - The element to remove.
 */
int doubly_linked_list_remove_element(doubly_linked_list_element** head, doubly_linked_list_element** tail, doubly_linked_list_element* element_to_remove);

#endif // DOUBLY_LINKED_LIST_H
