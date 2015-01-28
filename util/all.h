#ifndef ALL_H
#define ALL_H

/* Error checking, assertions, and invariants
 *
 * Error checking is the process of checking for invalid input.  Even bug free
 * code needs to perform error checking because the code cannot control what it
 * is given for input.  Code should handle invalid input gracefully.  Ideally
 * it should not crash or seg fault on any input and should print a helpful
 * error message.
 * 
 * Assertions check within the code for things that the code should be
 * guaranteeing to be true.  By definition, bug free code will never fire an
 * assertion on any input.  This code does not satisfy that condition if any
 * error checking is turned off.  Assertions are only needed for finding bugs
 * and can be compiled out of production code.
 * 
 * Invariants are conditions that should always be true.  Invariants are
 * usually written on complex data structures like a chanel network.
 * Invariants can be used for both error checking and assertions.  For example,
 * when reading in a channel network the code may check the invariant and throw
 * an error if the input mesh does not satisfy the invariant.  Thereafter, the
 * code may continue to check the invariant periodically as an assertion.
 * 
 * After error checking and assertions the rest of the code is written assuming
 * that those checks pass.  For example, if statement cases may be omitted if
 * they would be forbidden by an assertion.
 */
#define DEBUG_LEVEL_LIBRARY_ERRORS              (1 << 0) // Error checking of return values of libraries that this code calls.
#define DEBUG_LEVEL_USER_INPUT_SIMPLE           (1 << 1) // Error checking of simple conditions on user inputs such as file reads.
#define DEBUG_LEVEL_USER_INPUT_INVARIANTS       (1 << 2) // Error checking of invariants        on user inputs such as file reads.
#define DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE     (1 << 3) // Error checking of simple conditions on inputs to public functions.
#define DEBUG_LEVEL_PUBLIC_FUNCTIONS_INVARIANTS (1 << 4) // Error checking of invariants        on inputs to public functions.
#define DEBUG_LEVEL_ALL_FUNCTIONS_SIMPLE        (1 << 5) // Assertion      of simple conditions on inputs to all functions.
#define DEBUG_LEVEL_ALL_FUNCTIONS_INVARIANTS    (1 << 6) // Assertion      of invariants        on inputs to all functions.
#define DEBUG_LEVEL_INTERNAL_ASSERTIONS         (1 << 7) // Assertion      of conditions within the code of  all functions.

// No error checking or assertions.  Not recommended.
#define DEBUG_LEVEL_NONE (0)

// This is an appropriate debug level for production code released as an executable program.
#define DEBUG_LEVEL_PRODUCTION_EXECUTABLE (DEBUG_LEVEL_LIBRARY_ERRORS | DEBUG_LEVEL_USER_INPUT_SIMPLE | DEBUG_LEVEL_USER_INPUT_INVARIANTS)

// This is an appropriate debug level for production code called from other code as a library.
#define DEBUG_LEVEL_PRODUCTION_LIBRARY (DEBUG_LEVEL_LIBRARY_ERRORS | DEBUG_LEVEL_USER_INPUT_SIMPLE | DEBUG_LEVEL_USER_INPUT_INVARIANTS | \
                                        DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE | DEBUG_LEVEL_PUBLIC_FUNCTIONS_INVARIANTS)

// This is an appropriate debug level for development code.
#define DEBUG_LEVEL_DEVELOPMENT (DEBUG_LEVEL_LIBRARY_ERRORS | DEBUG_LEVEL_USER_INPUT_SIMPLE | DEBUG_LEVEL_USER_INPUT_INVARIANTS |                   \
                                 DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE | DEBUG_LEVEL_PUBLIC_FUNCTIONS_INVARIANTS | DEBUG_LEVEL_ALL_FUNCTIONS_SIMPLE | \
                                 DEBUG_LEVEL_INTERNAL_ASSERTIONS)

// This debug level is paranoid and expensive, but could help localize bugs.
#define DEBUG_LEVEL_ALL (DEBUG_LEVEL_LIBRARY_ERRORS | DEBUG_LEVEL_USER_INPUT_SIMPLE | DEBUG_LEVEL_USER_INPUT_INVARIANTS |                   \
                         DEBUG_LEVEL_PUBLIC_FUNCTIONS_SIMPLE | DEBUG_LEVEL_PUBLIC_FUNCTIONS_INVARIANTS | DEBUG_LEVEL_ALL_FUNCTIONS_SIMPLE | \
                         DEBUG_LEVEL_ALL_FUNCTIONS_INVARIANTS | DEBUG_LEVEL_INTERNAL_ASSERTIONS)

// Set this macro to the debug level to use.
#define DEBUG_LEVEL (DEBUG_LEVEL_DEVELOPMENT)

#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))

#ifndef TRUE
#define TRUE (1)
#endif // TRUE

#ifndef FALSE
#define FALSE (0)
#endif // FALSE

#define GRAVITY (9.81)    // Meters per second squared.
#define NaN     (0.0/0.0) // Not a number.

/* Special cases of element boundaries.
 * Used for mesh edges and channel upstream/downstream links.
 */
typedef enum
{
  NOFLOW  = -1, // This must be -1 because Triangle and TauDEM both use -1 to indicate no neighbor.
  INFLOW  = -2,
  OUTFLOW = -3
} boundary_condition_type;

#define IS_BOUNDARY(x) (NOFLOW == (x) || INFLOW == (x) || OUTFLOW == (x))

/* Used for the type of channel links.
 */
typedef enum
{
  UNUSED                             = 0, // Not all link numbers are used.
  STREAM                             = 1, // Stream link.
  WATERBODY                          = 2, // Waterbody link.
  ICEMASS                            = 3, // Equivalent to WATERBODY, but made of ice.
  FAKE                               = 4, // Used to create 3-way branches in a binary tree.
  INTERNAL_BOUNDARY_CONDITION        = 5, // For defining arbitrary flow constraints.
  PRUNED_STREAM                      = 6, // A pruned link was once used but is no longer needed.  Type is set to one of the pruned types,
  PRUNED_FAKE                        = 7, // and upstream and downstream links are broken, but other structure values are preserved.
  PRUNED_INTERNAL_BOUNDARY_CONDITION = 8
} link_type;

#endif // ALL_H
