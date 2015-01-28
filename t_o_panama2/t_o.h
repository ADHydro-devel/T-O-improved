#ifndef T_O_H
#define T_O_H

#include <pthread.h>

/* A t_o_parameters struct stores constant soil parameters for a Talbot-Ogden
 * domain.  It is pulled out as a separate struct from t_o_domain because
 * multiple domains being simulated may share the same parameters.
 */
typedef struct
{
  int             num_bins;                    // The number of bins.
  double          vg_alpha;                    // Van Genutchen parameter in one over meters.
  double          vg_n;                        // Van Genutchen parameter.
  double          bc_lambda;                   // Brook-Corey parameter.
  double          bc_psib;                     // Brook-Corey parameter in meters.
  double          delta_water_content;         // The water content width of each bin as a unitless fraction.
  double*         bin_water_content;           // 1D array containing the water content of each bin as a unitless fraction.
                                               // Varies between residual saturation and porosity.
  double*         cumulative_conductivity;     // 1D array. cumulative_conductivity[ii] contains the total conductivity of all bins up to and including bin ii
                                               // in units of meters of water per second.
  double          effective_capillary_suction; // The minimum value to use for bin_capilary_suction[last_bin] in the infiltrate_distance calculation.
  double*         bin_capillary_suction;       // The capillary suction head of each bin in meters.
  double          dry_depth_dt;                // The timestep in seconds used to calculate the values in bin_dry_depth.
  double*         bin_dry_depth;               // The dry depth of each bin in meters for a timestep of dry_depth_dt. This can be used as an exact value of the dry
                                               // depth for that timestep and an upper bound of the dry depth for smaller timesteps.
  pthread_mutex_t dry_depth_mutex;             // For thread-safe access to dry depth in shared parameters structures.
                                               // If THREAD_SAFE is not defined then this is uninitialized and unused.
  int             dry_depth_mutex_initialized; // Flag so that we know whether to destroy the mutex.
} t_o_parameters;

// FIXLATER possible optimization: Rather than searching for slugs to the left or right in contact with a given slug store pointers to them.
/* A slug struct represents a single slug of water in a single bin.
 * All of the slugs in a bin are stored in a doubly linked list.
 */
typedef struct slug slug;
struct slug
{
  slug*  prev; // The slug next closer to the surface or NULL if this is the top slug.
  slug*  next; // The slug next closer to the bottom  or NULL if this is the bottom slug.
  double top;  // The depth of the top of the slug in meters.
  double bot;  // The depth of the bottom of the slug in meters.
};

/* A t_o_domain struct stores all of the state of a single Talbot-Ogden domain.
 * This struct and the functions in this header should be taken together
 * like the member data and methods of a C++ object.
 */
typedef struct
{
  t_o_parameters* parameters;            // Constant soil parameters.
  double          layer_top_depth;       // The depth of the top of the t_o_domain in meters.
  double          layer_bottom_depth;    // The depth of the bottom of the t_o_domain in meters.
  double*         surface_front;         // 1D array containing the depth of the bottom of the surface front water in each bin in meters.
  slug**          top_slug;              // 1D array of pointers to the top    slug in each bin or NULL if the bin has no slugs.
  slug**          bot_slug;              // 1D array of pointers to the bottom slug in each bin or NULL if the bin has no slugs.
  int             yes_groundwater;       // Whether to simulate groundwater. If FALSE, groundwater_front is NULL.
  double*         groundwater_front;     // 1D array containing the depth of the top of the groundwater in meters in each bin.
                                         // Only used if yes_groundwater is TRUE.
  double          initial_water_content; // Bins with water content less than or equal to this are in contact with groundwater.
                                         // Only used if yes_groundwater is FALSE.
} t_o_domain;

/* Create a t_o_parameters struct and initialize it.
 * Return TRUE if there is an error, FALSE otherwise.
 *
 * Parameters:
 *
 * parameters          - A pointer passed by reference which will be assigned
 *                       to point to the newly allocated struct
 *                       or NULL if there is an error.
 * num_bins            - The number of bins.  One based indexing is used.
 * conductivity        - Hydrologic conductivity of the entire soil column
 *                       in meters of water per second [m/s].
 * porosity            - Soil porosity as a unitless fraction.
 * residual saturation - Soil residual saturation as a unitless fraction.
 * van_genutchen       - Flag to indicate if the Van Genutchen function
 *                       should be used to calculate the bin properties.
 *                       If FALSE, the Brook-Corey function is used instead.
 * vg_alpha            - Van Genutchen parameter in one over meters.
 * vg_n                - Van Genutchen parameter, unitless.
 * bc_lambda           - Brook-Corey parameter, unitless.
 * bc_psib             - Brook-Corey parameter in meters.
 */
int t_o_parameters_alloc(t_o_parameters** parameters, int num_bins, double conductivity, double porosity, double residual_saturation,
                         int van_genutchen, double vg_alpha, double vg_n, double bc_lambda, double bc_psib);

/* Free memory allocated by t_o_parameters_alloc.
 *
 * Parameters:
 *
 * parameters - A pointer to the t_o_parameters struct passed by reference.
 *              Will be set to NULL after the memory is deallocated.
 */
void t_o_parameters_dealloc(t_o_parameters** parameters);

/* Create a t_o_domain struct and initialize it.
 * Return TRUE if there is an error, FALSE otherwise.
 * If yes_groundwater is FALSE, then the domain is initialized to have no
 * water other than initial_water_content.  If yes_groundwater is TRUE and
 * initialize_to_hydrostatic is FALSE, then the domain is initialized to have
 * no water other than bin 1 and residual saturation.  If yes_groundwater is
 * TRUE and initialize_to_hydrostatic is TRUE, then the domain is initialized
 * to have no surface front water or slugs and groundwater is initialized to
 * hydrostatic equilibrium with the given water table.
 *
 * Parameters:
 *
 * domain                    - A pointer passed by reference which will be
 *                             assigned to point to the newly allocated struct
 *                             or NULL if there is an error.
 * parameters                - A pointer to the t_o_parameters struct.
 * layer_top_depth           - The depth of the top    of the t_o_domain in
 *                             meters.
 * layer_bottom_depth        - The depth of the bottom of the t_o_domain in
 *                             meters.
 * yes_groundwater           - Whether to simulate groundwater.
 * initial_water_content     - If yes_groundwater is FALSE, bins less than or
 *                             equal to initial_water_content are in contact
 *                             with deep groundwater, and bins greater than
 *                             initial_water_content have a free drainage
 *                             boundary condition at layer_depth.
 *                             Ignored if yes_groundwater is TRUE.
 * initialize_to_hydrostatic - Whether to initialize the groundwater front to
 *                             hydrostatic equilibrium with the given water
 *                             table.  Ignored if yes_groundwater is FALSE.
 * water_table               - The depth in meters of the water table to use
 *                             to initialize the groundwater front.
 *                             Ignored if yes_groundwater or
 *                             initialize_to_hydrostatic are FALSE.
 */
int t_o_domain_alloc(t_o_domain** domain, t_o_parameters* parameters, double layer_top_depth, double layer_bottom_depth, int yes_groundwater,
                     double initial_water_content, int initialize_to_hydrostatic, double water_table);

/* Free memory allocated for the Talbot-Ogden domain including memory
 * allocated by t_o_domain_alloc and memory subsequently allocated for slugs,
 * but excluding memory allocated by t_o_parameters_alloc because
 * the t_o_parameters struct might be shared.  You need to call
 * t_o_parameters_dealloc separately.
 *
 * Parameters:
 *
 * domain - A pointer to the t_o_domain struct passed by reference.
 *          Will be set to NULL after the memory is deallocated.
 */
void t_o_domain_dealloc(t_o_domain** domain);

/* Assert if any Talbot-Ogden domain invariant is violated.  In each bin
 * water must be non-overlapping and monotonically increasing in depth.
 * At every depth there can be no wet bin to the right of a dry bin.
 * This invariant may be violated temporarily between the steps of
 * a single timestep, but it will always hold after the redistribute step.
 *
 * Parameters:
 *
 * domain - A pointer to the t_o_domain struct.
 */
void t_o_check_invariant(t_o_domain* domain);

/* Return the total water in the Talbot-Ogden domain in meters of water.
 *
 * Parameters:
 *
 * domain - A pointer to the t_o_domain struct.
 */
double t_o_total_water_in_domain(t_o_domain* domain);

/* Step the Talbot-Ogden simulation forward one timestep.
 * Return TRUE if there is an error, FALSE otherwise.
 *
 * Parameters:
 *
 * domain               - A pointer to the t_o_domain struct.
 * dt                   - The duration of the timestep in seconds.
 * surfacewater_head    - The pressure head in meters of the surface water.
 *                        Normally, you will pass the same value as
 *                        surfacewater_depth, but this can be used to simulate
 *                        various laboratory test conditions.
 * surfacewater_depth   - A scalar passed by reference containing the depth
 *                        in meters of the surface water.  Will be updated for
 *                        the amount of infiltration.
 * water_table          - The depth in meters of the water table.  This can
 *                        also be used to simulate pressure head boundary
 *                        conditions at the botom of the domain, in which case
 *                        it is the negative of the pressure head.
 * groundwater_recharge - A scalar passed by reference containing any
 *                        previously accumulated groundwater recharge in meters
 *                        of water.  Will be updated for the amount of water
 *                        that flowed between the Talbot-Ogden domain and
 *                        groundwater.  Positive means water flowed down out of
 *                        the Talbot-Ogden domain.  Negative means water flowed
 *                        up in to the Talbot-Ogden domain.
 */
int t_o_timestep(t_o_domain* domain, double dt, double surfacewater_head, double* surfacewater_depth, double water_table, double* groundwater_recharge);

/* Arbitrarily add water to the groundwater front of a Talbot-Ogden domain.
 * This function is used to couple the domain to a separate groundwater
 * simulation.  Some groundwater simulations work by assuming that the
 * groundwater is effectively in a set of open topped buckets.  If water
 * moves in to a bucket you can put that water on top of the water already
 * in the bucket.
 *
 * However, when you couple that kind of groundwater simulation to a Talbot-
 * Ogden domain, the groundwater cannot be considered to be in an open topped
 * bucket.  The vadose zone water in the Talbot-Ogden domain is on top of the
 * groundwater, and there cannot be a gap or overlap between the two.  In that
 * case, when the groundwater simulation says that water moves in to a bucket,
 * you have to put that water on top of the Talbot-Ogden domain groundwater
 * front, and you want to do it in a way that does not grossly violate the
 * behavior of the Talbot-Ogden domain.
 *
 * t_o_add_groundwater adds water starting with the lowest point of the
 * groundwater front.  It fills in the lowest empty space creating a flat top
 * across all bins where it adds water.  It finds the depth needed to add the
 * requested quantity of water filling in all of the empty space below that
 * depth.
 *
 * groundwater_recharge is passed by reference because it is an in/out
 * parameter.  You pass in the amount of water you want to add and it gets
 * updated to the amount that was not able to be added.  If the updated value
 * of groundwater_recharge is not epsilon equal to zero then the entire domain
 * is full and you should move the water table to the surface and add the
 * remaining water to the surface water.
 * 
 * If you call t_o_add_groundwater on a domain with yes_groundwater FALSE
 * nothing is done.
 *
 * Parameters:
 *
 * domain               - A pointer to the t_o_domain struct.
 * groundwater_recharge - A scalar passed by reference containing the amount of
 *                        water to add to the domain in meters of water.  Will
 *                        be updated to the amount of water that was not able
 *                        to be added to the domain.  If the updated value is
 *                        epsilon equal to zero all of the water was added
 *                        successfully.  If the updated value is not epsilon
 *                        equal to zero the caller is responsible for putting
 *                        that much water somewhere else to maintain mass
 *                        conservation.  If it failed to add all of the water
 *                        that means the domain is entirely full and you can
 *                        move the water table to the surface and put the rest
 *                        of the water in the surface water.
 */
void t_o_add_groundwater(t_o_domain* domain, double* groundwater_recharge);

/* Arbitrarily remove water from the groundwater front of a Talbot-Ogden
 * domain.  This function is used to couple the domain to a separate
 * groundwater simulation.  See t_o_add_groundwater for more information.
 *
 * t_o_take_groundwater takes water starting with the rightmost bin.  It lowers
 * the groundwater front as low as 10% of capilary head above the water table.
 * If it can't get enough water from the rightmost bin it proceeds to bins to
 * the left until it can get all of the water or every bin is at 10% of
 * capillary head above the water table.
 * 
 * groundwater_recharge is passed by reference because it is an in/out
 * parameter.  You pass in the amount of water you want to take as a negative
 * number and it gets updated to the amount that was not able to be taken.
 * 
 * If you call t_o_take_groundwater on a domain with yes_groundwater FALSE
 * nothing is done.
 *
 * Parameters:
 *
 * domain               - A pointer to the t_o_domain struct.
 * water_table          - The depth in meters of the water table.
 * groundwater_recharge - A scalar passed by reference containing the amount of
 *                        water to take from the domain in meters of water as a
 *                        negative number.  Will be updated to the amount of
 *                        water that was not able to be taken from the domain.
 */
void t_o_take_groundwater(t_o_domain* domain, double water_table, double* groundwater_recharge);

/* Return an estimate of the specific yield.
 * 
 * Parameters:
 * 
 * domain      - A pointer to the t_o_domain struct.
 * water_table - The depth in meters of the water table.
 */
double t_o_specific_yield(t_o_domain* domain, double water_table);

// FIXME, wencong, add ET, Dec. 10, 2014. A very simple ET function.
/*
   domain      - A pointer to the t_o_domain struct.
   dt          - Time step size in seconds.
   root_depth  - Root depth in meters, if it's zero, it's bare soil evaporation, no transpiration.
   PET         - Potential evapotranspiration in meters per second.
   field_capacity         - Field capacity, unitless.
   wilting_point          - Wilting point, unitless.
   use_feddes             - whether to use the linear function of pressure head as Feddes, or use linear function of (field_capacity, wilting_point) 
                            to caculate actual ET.
   field_capacity_suction - Use if use_feddes is TRUE, absolute pressure coresponding to field_capacity. 
   wilting_point_suction  - Use if use_feddes is TRUE, absolute pressure coresponding to wilting_point. 
   surfacewater_depth     - A scalar passed by reference contains surfacewater_depth in meters of water.
   evaporated_water       - A scalar passed by reference contains evaporated water   in meters of water.
 */
int t_o_ET(t_o_domain* domain, double dt, double root_depth, double PET, double field_capacity, double wilting_point, 
           int use_feddes, double field_capacity_suction, double wilting_point_suction, double* surfacewater_depth, double* evaporated_water);

#endif // T_O_H
