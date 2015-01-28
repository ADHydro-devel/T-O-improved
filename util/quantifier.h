#ifndef QUANTIFIER_H
#define QUANTIFIER_H

/* Return a == *b.
 * Return FALSE if b is NULL.
 *
 * Parameters:
 *
 * a - One of the integers.
 * b - A pointer to the other integer cast as a void pointer.
 */
int i_test_equals(int a, void* b);

/* Existential quantifier for an integer array.
 * Return TRUE if there exists an element of array for which
 * test(array[ii], test_params) returns TRUE.
 * Return FALSE if array has no elements (array == NULL || 1 > size) or there
 * is no test (NULL == test).
 *
 * Parameters:
 *
 * array       - The array of integers to test.
 * size        - The number of elements in array.  One based indexing is used.
 * test        - A function pointer to the test.
 * test_params - Passed in to test.  test must know what type of pointer
 *               this is and what to do with the value it points to.
 */
int i_exists(int* array, int size, int(*test)(int, void*), void* test_params);

/* Universal quantifier for an integer array.
 * Return TRUE if for all elements of array
 * test(array[ii], test_params) returns TRUE.
 * Return TRUE if array has no elements (array == NULL || 1 > size) or there
 * is no test (NULL == test).
 *
 * Parameters:
 *
 * array       - The array of integers to test.
 * size        - The number of elements in array.  One based indexing is used.
 * test        - A function pointer to the test.
 * test_params - Passed in to test.  test must know what type of pointer
 *               this is and what to do with the value it points to.
 */
int i_forall(int* array, int size, int(*test)(int, void*), void* test_params);

#endif // QUANTIFIER_H
