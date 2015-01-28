#ifndef MEMFUNC_H
#define MEMFUNC_H

/* Allocate memory and initialize all allocated bytes to zero.
 * Return TRUE if there is an error, FALSE otherwise.
 *
 * Parameters:
 *
 * ptr   - A pointer passed by reference which will be assigned to point
 *         to the newly allocated memory or NULL if there is an error.
 * bytes - The number of bytes to allocate.
 */
int v_alloc(void** ptr, int bytes);

/* Free memory allocated by v_alloc.
 * Return TRUE if there is an error, FALSE otherwise.
 * Even if there is an error make every effort to free as much as possible.
 *
 * Parameters:
 *
 * ptr   - A pointer to the memory to deallocate passed by reference.
 *         Will be set to NULL after the memory is deallocated.
 * bytes - The number of bytes in the allocated memory.
 */
int v_dealloc(void** ptr, int bytes);

/* Allocate a 1D array of integers and initialize all elements to zero.
 * Return TRUE if there is an error, FALSE otherwise.
 *
 * Parameters:
 *
 * array - A pointer passed by reference which will be assigned to point
 *         to the newly allocated array or NULL if there is an error.
 * size  - The number of elements in the array.
 *         One extra element will be allocated to support one based indexing.
 */
int i_alloc(int** array, int size);

/* Free memory allocated by i_alloc.
 * Return TRUE if there is an error, FALSE otherwise.
 * Even if there is an error make every effort to free as much as possible.
 *
 * Parameters:
 *
 * array - A pointer to the array to deallocate passed by reference.
 *         Will be set to NULL after the array is deallocated.
 * size  - The number of elements in the allocated array not including
 *         the extra element added to support one based indexing.
 */
int i_dealloc(int** array, int size);

/* Allocate a 1D array of doubles and initialize all elements to zero.
 * Return TRUE if there is an error, FALSE otherwise.
 *
 * Parameters:
 *
 * array - A pointer passed by reference which will be assigned to point
 *         to the newly allocated array or NULL if there is an error.
 * size  - The number of elements in the array.
 *         One extra element will be allocated to support one based indexing.
 */
int d_alloc(double** array, int size);

/* Free memory allocated by d_alloc.
 * Return TRUE if there is an error, FALSE otherwise.
 * Even if there is an error make every effort to free as much as possible.
 *
 * Parameters:
 *
 * array - A pointer to the array to deallocate passed by reference.
 *         Will be set to NULL after the array is deallocated.
 * size  - The number of elements in the allocated array not including
 *         the extra element to support one based indexing.
 */
int d_dealloc(double** array, int size);

/* Allocate an array of memory blocks as an array of pointers each pointing to
 * a newly allocated block of memory and initialize all bytes in each
 * allocated block to zero.
 * Return TRUE if there is an error, FALSE otherwise.
 *
 * Parameters:
 *
 * array - A pointer passed by reference which will be assigned to point
 *         to the newly allocated array or NULL if there is an error.
 * rows  - The number of rows to allocate in the array.  NOTE: THERE WILL
 *         NOT BE ONE EXTRA ROW ALLOCATED TO SUPPORT ONE BASED INDEXING.
 *         IT IS EXPECTED THAT THIS FUNCTION WILL BE CALLED BY
 *         A TYPE SPECIFIC FUNCTION THAT WILL ADD THE EXTRA ROW.
 * bytes - The number of bytes to allocate in each memory block in the array.
 */
int vtwo_alloc(void*** array, int rows, int bytes);

/* Free memory allocated by vtwo_alloc.
 * Return TRUE if there is an error, FALSE otherwise.
 * Even if there is an error make every effort to free as much as possible.
 *
 * Parameters:
 *
 * array - A pointer to the array to deallocate passed by reference.
 *         Will be set to NULL after the array is deallocated.
 * rows  - The number of rows in the allocated array including any extra row
 *         added to support one based indexing.
 * bytes - The number of bytes allocated in each memory block in the array.
 */
int vtwo_dealloc(void*** array, int rows, int bytes);

/* Allocate a 2D array of integers as an array of pointers to 1D arrays and
 * initialize all elements to zero.
 * Return TRUE if there is an error, FALSE otherwise.
 *
 * Parameters:
 *
 * array - A pointer passed by reference which will be assigned to point
 *         to the newly allocated array or NULL if there is an error.
 * rows  - The number of rows to allocate in the array.
 *         One extra row will be allocated to support one based indexing.
 * cols  - The number of columns to allocate in the array.
 *         One extra column will be allocated to support one based indexing.
 */
int itwo_alloc(int*** array, int rows, int cols);

/* Free memory allocated by itwo_alloc.
 * Return TRUE if there is an error, FALSE otherwise.
 * Even if there is an error make every effort to free as much as possible.
 *
 * Parameters:
 *
 * array - A pointer to the array to deallocate passed by reference.
 *         Will be set to NULL after the array is deallocated.
 * rows  - The number of rows in the allocated array not including
 *         the extra row added to support one based indexing.
 * cols  - The number of columns in the allocated array not including
 *         the extra column added to support one based indexing.
 */
int itwo_dealloc(int*** array, int rows, int cols);

/* Allocate a 2D array of doubles as an array of pointers to 1D arrays and
 * initialize all elements to zero.
 * Return TRUE if there is an error, FALSE otherwise.
 *
 * Parameters:
 *
 * array - A pointer passed by reference which will be assigned to point
 *         to the newly allocated array or NULL if there is an error.
 * rows  - The number of rows to allocate in the array.
 *         One extra row will be allocated to support one based indexing.
 * cols  - The number of columns to allocate in the array.
 *         One extra column will be allocated to support one based indexing.
 */
int dtwo_alloc(double*** array, int rows, int cols);

/* Free memory allocated by dtwo_alloc.
 * Return TRUE if there is an error, FALSE otherwise.
 * Even if there is an error make every effort to free as much as possible.
 *
 * Parameters:
 *
 * array - A pointer to the array to deallocate passed by reference.
 *         Will be set to NULL after the array is deallocated.
 * rows  - The number of rows in the allocated array not including
 *         the extra row added to support one based indexing.
 * cols  - The number of columns in the allocated array not including
 *         the extra column added to support one based indexing.
 */
int dtwo_dealloc(double*** array, int rows, int cols);

/* Print a warning if any allocated memory has not been freed and optionally
 * free all allocated memory.
 * Return TRUE if there is an error, FALSE otherwise.
 * If there is an error before freeing memory do not free memory.
 * Call this at the end of your program after you think you have deallocated
 * everything.  Do not use this as a substitute for properly freeing memory.
 *
 * Parameters:
 *
 * cleanup - If TRUE, free all allocated memory.
 */
int check_all_memory_freed(int cleanup);

#endif // MEMFUNC_H
