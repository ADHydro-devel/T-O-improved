#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "t_o.h"
#include "doubly_linked_list.h"
#include "epsilon.h"
#include "memfunc.h"
#include "all.h"

#define SLIVER_SLUG_SIZE (0.001) // Meters.

#define THREAD_SAFE // Leave this defined to have the code use mutexes to be thread safe.

#ifdef THREAD_SAFE
static pthread_mutex_t slug_pool_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif // THREAD_SAFE

static slug* slug_pool = NULL; // A linked list of unused slug structs so that we don't have to allocate and deallocate every time.
                               // The list is singly linked.  Only the next pointers are used.

/* Documented assumptions:
 *
 * This code implements a single layer Talbot-Ogden domain. At the top of the
 * domain, the head boundary condition and water availability are both
 * determined by the ponded depth of surface water.  At the bottom of the
 * domain, the head boundary condition must be given and water availability is
 * assumed to be unlimited.  That is, any amount of water can flow out of or in
 * to the domain without changing the head boundary condition.
 *
 * If yes_groundwater is TRUE then there is groundwater at the bottom of every
 * bin.  The groundwater is not allowed to drop below layer_bottom_depth.  If
 * the equations would normally make groundwater drop below layer_bottom_depth,
 * the groundwater level is pinned at layer_bottom_depth, and a warning is
 * printed that the simulation will be inaccurate unless layer_depth is
 * increased.  As a consequence, every bin is always wet at layer_bottom_depth.
 *
 * If yes_groundwater is FALSE then the bins with water content less than or
 * equal to initial_water_content are in contact with deep groundwater, and all
 * other bins have a free drainage condition at layer_bottom_depth.  Surface
 * front water and slugs can reach layer_bottom_depth and continue to flow
 * downward.  The surface front water or slug is chopped off at
 * layer_bottom_depth, and the water removed from the domain is accounted for.
 * Surface front water or slugs that extend to layer_bottom_depth are not in
 * contact with groundwater.
 *
 * Bin 1 is forced to be completely full of water.  In any realistic simulation
 * this will be the case.  If it isn't then residual saturation is set too
 * high.  To achieve this, if yes_groundwater is TRUE then in the groundwater
 * simulation step if the groundwater in bin 1 wants to drop below the surface
 * it will be pinned to the surface and a warning will be printed.  If
 * yes_groundwater is FALSE and initial_water_content is less than the water
 * content of bin 1 it is set equal to the water content of bin 1 and a warning
 * is printed.  In any case, first_bin is forced to be at least 2, and last_bin
 * is forced to be at least 1.
 */

// FIXLATER possible optimization eliminate tiny slugs.

/* Comment in .h file. */
int t_o_parameters_alloc(t_o_parameters** parameters, int num_bins, double conductivity, double porosity, double residual_saturation,
                         int van_genutchen, double vg_alpha, double vg_n, double bc_lambda, double bc_psib)
{
  int    error = FALSE; // Error flag.
  int    ii;            // Loop counter.
  double m, p;          // Derived parameters.

  if (NULL == parameters)
    {
      fprintf(stderr, "ERROR: parameters must not be NULL\n");
      error = TRUE;
    }
  else
    {
      *parameters = NULL; // Prevent deallocating a random pointer.
    }

  if (0 >= num_bins)
    {
      fprintf(stderr, "ERROR: num_bins must be greater than 0\n");
      error = TRUE;
    }

  if (0.0 >= conductivity)
    {
      fprintf(stderr, "ERROR: conductivity must be greater than 0\n");
      error = TRUE;
    }

  if (0.0 >= residual_saturation)
    {
      fprintf(stderr, "ERROR: residual_saturation must be greater than 0\n");
      error = TRUE;
    }

  if (residual_saturation >= porosity)
    {
      fprintf(stderr, "ERROR: porosity must be greater than residual_saturation\n");
      error = TRUE;
    }

  if (van_genutchen)
    {
      // If given Van Genutchen parameters calculate bc_lambda and bc_psib.
      if (1.0 < vg_n)
        {
          m = 1.0 - 1.0 / vg_n;
          p = 1.0 + 2.0 / m;
          bc_lambda = 2 / (p - 3.0);
          bc_psib = (p + 3.0) * (147.8 + 8.1 * p + 0.092 * p * p) / (2.0 * vg_alpha * p * (p - 1.0) * (55.6 + 7.4 * p + p * p));
          assert(0.0 < bc_psib);
        }
      else
        {
          fprintf(stderr, "ERROR: Van Genutchen parameter n must be greater than 1\n");
          error = TRUE;
        }
    }
  else
    {
      // If given Brook-Corey parameters calculate vg_alpha and vg_n.
      if (0.0 < bc_psib)
        {
          p = 3.0 + 2.0 / bc_lambda;
          m = 2.0 / (p - 1.0);
          vg_alpha = (p + 3.0) * (147.8 + 8.1 * p + 0.092 * p * p) / (2.0 * bc_psib * p * (p - 1.0) * (55.6 + 7.4 * p + p * p));
          vg_n = 1.0 / (1.0 - m);
        }
      else
        {
          fprintf(stderr, "ERROR: Brook-Corey parameter psib must be greater than 0\n");
          error = TRUE;
        }
    }

  // Allocate the t_o_parameters struct.
  if (!error)
    {
      error = v_alloc((void**)parameters, sizeof(t_o_parameters));
    }

  // Initialize the t_o_parameters struct.
  if (!error)
    {
      (*parameters)->num_bins                = num_bins;
      (*parameters)->vg_alpha                = vg_alpha;
      (*parameters)->vg_n                    = vg_n;
      (*parameters)->bc_lambda               = bc_lambda;
      (*parameters)->bc_psib                 = bc_psib;
      (*parameters)->delta_water_content     = (porosity - residual_saturation) / num_bins;
      (*parameters)->bin_water_content       = NULL;
      (*parameters)->cumulative_conductivity = NULL;

      if (van_genutchen)
        {
          (*parameters)->effective_capillary_suction = (1.0 / vg_alpha) * (0.046 * m + 2.07 * m * m + 19.5 * m * m * m) / (1 + 4.7 * m + 16.0 * m * m);
        }
      else
        {
          (*parameters)->effective_capillary_suction = bc_psib * (2.0 + 3.0 * bc_lambda) / (1.0 + 3.0 * bc_lambda);
        }

      (*parameters)->bin_capillary_suction       = NULL;
      (*parameters)->dry_depth_dt                = 0.0;
      (*parameters)->bin_dry_depth               = NULL;
      (*parameters)->dry_depth_mutex_initialized = FALSE;
    }

  // Allocate bin_water_content.
  if (!error)
    {
      error = d_alloc(&(*parameters)->bin_water_content, num_bins);
    }

  // Initialize bin_water_content.
  if (!error)
    {
      for (ii = 1; ii <= num_bins; ii++)
        {
          (*parameters)->bin_water_content[ii] = residual_saturation + ii * (*parameters)->delta_water_content;
        }
        
      assert(epsilon_equal(porosity, (*parameters)->bin_water_content[num_bins]));
    }

  // 1D array containing the relative saturation of the center of each bin as a unitless fraction of the
  // available pore space between residual saturation and porosity.  Varies between 0 and 1.
  double bin_relative_saturation[num_bins + 1];

  // Calculate bin_relative_saturation.
  if (!error)
    {
      for (ii = 1; ii <= num_bins; ii++)
        {
          bin_relative_saturation[ii] = ((*parameters)->bin_water_content[ii] - residual_saturation) / (porosity - residual_saturation);
        }

      assert(1.0 == bin_relative_saturation[num_bins]);
    }

  // Allocate cumulative_conductivity.
  if (!error)
    {
      error = d_alloc(&(*parameters)->cumulative_conductivity, num_bins);
    }

  // Initialize cumulative_conductivity.
  if (!error)
    {
      for (ii = 1; ii <= num_bins; ii++)
        {
          if (van_genutchen)
            {
              (*parameters)->cumulative_conductivity[ii] = conductivity * pow(bin_relative_saturation[ii], 0.5) *
                  pow(1.0 - pow(1.0 - pow(bin_relative_saturation[ii], 1.0 / m), m), 2.0);
            }
          else
            {
              (*parameters)->cumulative_conductivity[ii] = conductivity * pow(bin_relative_saturation[ii], 3.0 + 2.0 / bc_lambda);
            }
        }

      assert(conductivity == (*parameters)->cumulative_conductivity[num_bins]);
    }

  // Allocate bin_capillary_suction.
  if (!error)
    {
      error = d_alloc(&(*parameters)->bin_capillary_suction, num_bins);
    }

  // Initialize bin_capillary_suction.
  if (!error)
    {
      for (ii = 1; ii <= num_bins; ii++)
        {
          if (van_genutchen)
            {
              double relative_saturation = bin_relative_saturation[ii];
              
              // Van Genutchen capillary suction goes to zero at 100% relative saturation so ise the mid-bin value for the last bin.
              if (ii == num_bins)
                {
                  relative_saturation = (bin_relative_saturation[ii] + bin_relative_saturation[ii - 1]) * 0.5;
                }
              
              (*parameters)->bin_capillary_suction[ii] = (1.0 / vg_alpha) * pow(pow(1.0 / relative_saturation, 1.0 / m) - 1.0, 1.0 / vg_n);
            }
          else
            {
              (*parameters)->bin_capillary_suction[ii] = bc_psib * pow(bin_relative_saturation[ii], -1.0 / bc_lambda);
            }
        }

      assert(0.0 < (*parameters)->bin_capillary_suction[num_bins]);
    }

  // Allocate bin_dry_depth.
  if (!error)
    {
      error = d_alloc(&(*parameters)->bin_dry_depth, num_bins);
    }

  // Initialize bin_dry_depth.
  if (!error)
    {
      for (ii = 1; ii <= num_bins; ii++)
        {
          (*parameters)->bin_dry_depth[ii] = 0.0;
        }
    }

  // Initialize dry_depth_mutex
#ifdef THREAD_SAFE
  if (!error)
    {
      error = pthread_mutex_init(&(*parameters)->dry_depth_mutex, NULL);
      
      if (!error)
        {
          (*parameters)->dry_depth_mutex_initialized = TRUE;
        }
    }
#endif // THREAD_SAFE

  if (error)
    {
      t_o_parameters_dealloc(parameters);
    }

  return error;
}

/* Comment in .h file. */
void t_o_parameters_dealloc(t_o_parameters** parameters)
{
  assert(NULL != parameters);

  if (NULL != parameters && NULL != *parameters)
    {
      // Deallocate bin_water_content.
      if (NULL != (*parameters)->bin_water_content)
        {
          d_dealloc(&(*parameters)->bin_water_content, (*parameters)->num_bins);
        }

      // Deallocate cumulative_conductivity.
      if (NULL != (*parameters)->cumulative_conductivity)
        {
          d_dealloc(&(*parameters)->cumulative_conductivity, (*parameters)->num_bins);
        }

      // Deallocate bin_capillary_suction.
      if (NULL != (*parameters)->bin_capillary_suction)
        {
          d_dealloc(&(*parameters)->bin_capillary_suction, (*parameters)->num_bins);
        }

      // Deallocate bin_dry_depth.
      if (NULL != (*parameters)->bin_dry_depth)
        {
          d_dealloc(&(*parameters)->bin_dry_depth, (*parameters)->num_bins);
        }

#ifdef THREAD_SAFE
      if ((*parameters)->dry_depth_mutex_initialized)
        {
          pthread_mutex_destroy(&(*parameters)->dry_depth_mutex);
        }
#endif // THREAD_SAFE
      
      // Deallocate the t_o_parameters struct.
      v_dealloc((void**)parameters, sizeof(t_o_parameters));
    }
}

/* Comment in .h file. */
int t_o_domain_alloc(t_o_domain** domain, t_o_parameters* parameters,double layer_top_depth, double layer_bottom_depth, int yes_groundwater,
                     double initial_water_content, int initialize_to_hydrostatic, double water_table)
{
  int error = FALSE; // Error flag.
  int ii;            // Loop counter.

  if (NULL == domain)
    {
      fprintf(stderr, "ERROR: domain must not be NULL\n");
      error = TRUE;
    }
  else
    {
      *domain = NULL; // Prevent deallocating a random pointer.
    }

  if (NULL == parameters)
    {
      fprintf(stderr, "ERROR: parameters must not be NULL\n");
      error = TRUE;
    }

  if (0.0 > layer_top_depth)
    {
      fprintf(stderr, "ERROR: layer_top_depth must be greater than or equal to 0\n");
      error = TRUE;
    }
    
  if (layer_top_depth >= layer_bottom_depth)
    {
      fprintf(stderr, "ERROR: layer_bottom_depth must be greater than layer_top_depth\n");
      error = TRUE;
    }

  if (!yes_groundwater && 0.0 >= initial_water_content)
    {
      fprintf(stderr, "ERROR: initial_water_content must be greater than 0\n");
      error = TRUE;
    }

  if (yes_groundwater && initialize_to_hydrostatic && 0.0 > water_table)
    {
      fprintf(stderr, "ERROR: water_table must be greater than or equal to 0\n");
      error = TRUE;
    }

  // Allocate the t_o_domain struct.
  if (!error)
    {
      error = v_alloc((void**)domain, sizeof(t_o_domain));
    }

  // Initialize the t_o_domain struct.
  if (!error)
    {
      (*domain)->parameters = parameters;
      (*domain)->layer_top_depth    = layer_top_depth;
      (*domain)->layer_bottom_depth = layer_bottom_depth;
      (*domain)->surface_front = NULL;
      (*domain)->top_slug = NULL;
      (*domain)->bot_slug = NULL;
      (*domain)->yes_groundwater = yes_groundwater;
      (*domain)->groundwater_front = NULL;
      if (!yes_groundwater)
        {
          if (initial_water_content >= parameters->bin_water_content[1])
            {
              (*domain)->initial_water_content = initial_water_content;
            }
          else
            {
              fprintf(stderr, "WARNING: initial_water_content less than the water content of bin 1, using the water content of bin 1 instead.  "
                              "You should set residual_saturation lower.\n");
              (*domain)->initial_water_content = parameters->bin_water_content[1];
            }
        }
    }

  // Allocate surface_front.
  if (!error)
    {
      error = d_alloc(&(*domain)->surface_front, parameters->num_bins);
    }

  // Initialize surface_front.
  if (!error)
    {
      for (ii = 1; ii <= parameters->num_bins; ii++)
        {
          (*domain)->surface_front[ii] = layer_top_depth;
        }
    }

  // Allocate top_slug.
  if (!error)
    {
      error = v_alloc((void**)&(*domain)->top_slug, (parameters->num_bins + 1) * sizeof(slug*));
    }

  // Initialize top_slug.
  if (!error)
    {
      for (ii = 1; ii <= parameters->num_bins; ii++)
        {
          (*domain)->top_slug[ii] = NULL;
        }
    }

  // Allocate bot_slug.
  if (!error)
    {
      error = v_alloc((void**)&(*domain)->bot_slug, (parameters->num_bins + 1) * sizeof(slug*));
    }

  // Initialize bot_slug.
  if (!error)
    {
      for (ii = 1; ii <= parameters->num_bins; ii++)
        {
          (*domain)->bot_slug[ii] = NULL;
        }
    }

  if (yes_groundwater)
    {
      // Allocate groundwater_front.
      if (!error)
        {
          error = d_alloc(&(*domain)->groundwater_front, parameters->num_bins);
        }

      // Initialize groundwater_front.
      if (!error)
        {
          // Force bin 1 to be full of water.
          (*domain)->groundwater_front[1] = layer_top_depth;

          for (ii = 2; ii <= parameters->num_bins; ii++)
            {
              if (initialize_to_hydrostatic)
                {
                  (*domain)->groundwater_front[ii] = water_table - parameters->bin_capillary_suction[ii];

                  if ((*domain)->groundwater_front[ii] < layer_top_depth)
                    {
                      (*domain)->groundwater_front[ii] = layer_top_depth;
                    }
                  else if ((*domain)->groundwater_front[ii] > layer_bottom_depth)
                    {
                      (*domain)->groundwater_front[ii] = layer_bottom_depth;
                    }
                }
              else
                {
                  (*domain)->groundwater_front[ii] = layer_bottom_depth;
                }
            }
        }
    }

  if (error)
    {
      t_o_domain_dealloc(domain);
    }

  return error;
}

/* Comment in .h file. */
void t_o_domain_dealloc(t_o_domain** domain)
{
  int ii; // Loop counter.

  assert(NULL != domain);

  if (NULL != domain && NULL != *domain)
    {
      // Deallocate surface_front.
      if (NULL != (*domain)->surface_front)
        {
          d_dealloc(&(*domain)->surface_front, (*domain)->parameters->num_bins);
        }

      // Deallocate top_slug.
      if (NULL != (*domain)->top_slug)
        {
          // Deallocate slugs.
          for(ii = 0; ii <= (*domain)->parameters->num_bins; ii++)
            {
              slug* temp_slug = (*domain)->top_slug[ii];

              while (NULL != temp_slug)
                {
                  slug* next_slug = temp_slug->next;

                  v_dealloc((void**)&temp_slug, sizeof(slug));
                  temp_slug = next_slug;
                }
            }

          v_dealloc((void**)&(*domain)->top_slug, ((*domain)->parameters->num_bins + 1) * sizeof(slug*));
        }

      // Deallocate bot_slug.
      if (NULL != (*domain)->bot_slug)
        {
          v_dealloc((void**)&(*domain)->bot_slug, ((*domain)->parameters->num_bins + 1) * sizeof(slug*));
        }

      // Deallocate groundwater_front.
      if (NULL != (*domain)->groundwater_front)
        {
          d_dealloc(&(*domain)->groundwater_front, (*domain)->parameters->num_bins);
        }

      // Deallocate the t_o_domain struct.
      v_dealloc((void**)domain, sizeof(t_o_domain));
    }

  // When you call this function deallocate all unused slug structs in the slug pool.
#ifdef THREAD_SAFE
  pthread_mutex_lock(&slug_pool_mutex);
#endif // THREAD_SAFE
  
  while (NULL != slug_pool)
    {
      slug* temp_slug = slug_pool;

      slug_pool = slug_pool->next;
      v_dealloc((void**)&temp_slug, sizeof(slug));
    }
  
#ifdef THREAD_SAFE
  pthread_mutex_unlock(&slug_pool_mutex);
#endif // THREAD_SAFE
}

/* Return TRUE if the given bin is completely wet from top to bot,
 * FALSE otherwise.
 *
 * Parameters:
 *
 * domain - A pointer to the t_o_domain struct.
 * bin    - Which bin to check for water.  One based indexing is used.
 * top    - The top of the region to check for water.
 * bot    - The bottom of the region to check for water.
 */
int has_water_at_depth(t_o_domain* domain, int bin, double top, double bot)
{
  int has_water = FALSE;
  
  assert(NULL != domain && 0 < bin && bin <= domain->parameters->num_bins && domain->layer_top_depth <= top && top <= bot && bot <= domain->layer_bottom_depth);

  if (!domain->yes_groundwater && domain->parameters->bin_water_content[bin] <= domain->initial_water_content)
    {
      // The bin is completely saturated.
      has_water = TRUE;
    }
  else if (domain->surface_front[bin] > domain->layer_top_depth && domain->surface_front[bin] >= bot)
    {  
      // There is surface front water from top to bot.
      has_water = TRUE;
    }
  else if (domain->yes_groundwater && domain->groundwater_front[bin] <= top)
    {
      // There is groundwater from top to bot.
      has_water = TRUE;
    }
  else
    {
      slug* temp_slug = domain->top_slug[bin];

      while (NULL != temp_slug && temp_slug->bot < top)
        {
          temp_slug = temp_slug->next;
        }

      if (NULL != temp_slug && temp_slug->top <= top && temp_slug->bot >= bot)
        {
          // There is water in a slug from top to bot.
          has_water = TRUE;
        }
    }

  return has_water;
}

/* Comment in .h file. */
void t_o_check_invariant(t_o_domain* domain)
{
#ifndef NDEBUG
  int ii; // Loop counter.

  assert(NULL != domain);

  if (NULL != domain)
    {
      // Process all bins.
      for (ii = 1; ii <= domain->parameters->num_bins; ii++)
        {
          if ((domain->yes_groundwater && domain->layer_top_depth == domain->groundwater_front[ii]) ||
              (!domain->yes_groundwater && domain->parameters->bin_water_content[ii] <= domain->initial_water_content))
            {
              // The bin is completely saturated.
              assert(domain->layer_top_depth == domain->surface_front[ii] && NULL == domain->top_slug[ii] && NULL == domain->bot_slug[ii]);

              // Water to the left.
              if (1 < ii)
                {
                  assert(has_water_at_depth(domain, ii - 1, domain->layer_top_depth, domain->layer_bottom_depth));
                }
            }
          else // The bin is not completely saturated.
            {
              // Surface front within domain.
              assert(domain->layer_top_depth <= domain->surface_front[ii]);

              // Water to the left of surface front water.
              if (1 < ii && domain->layer_top_depth < domain->surface_front[ii])
                {
                  assert(has_water_at_depth(domain, ii - 1, domain->layer_top_depth, domain->surface_front[ii]));
                }

              if (domain->yes_groundwater)
                {
                  // Surface front less than groundwater and groundwater within domain.
                  assert(domain->surface_front[ii] < domain->groundwater_front[ii] && domain->groundwater_front[ii] <= domain->layer_bottom_depth);

                  // Water to the left of groundwater.
                  if (1 < ii)
                    {
                      assert(has_water_at_depth(domain, ii - 1, domain->groundwater_front[ii], domain->layer_bottom_depth));
                    }
                }
              else
                {
                  // Surface front within domain.
                  assert(domain->surface_front[ii] <= domain->layer_bottom_depth);
                }

              // Check slugs.
              if (NULL == domain->top_slug[ii])
                {
                  assert(NULL == domain->bot_slug[ii]);
                }
              else
                {
                  // Linked list consistency.
                  assert(NULL != domain->bot_slug[ii] && NULL == domain->top_slug[ii]->prev && NULL == domain->bot_slug[ii]->next);

                  // Surface front less than slugs.
                  assert(domain->surface_front[ii] < domain->top_slug[ii]->top);

                  if (domain->yes_groundwater)
                    {
                      // Groundwater greater than slugs.
                      assert(domain->bot_slug[ii]->bot < domain->groundwater_front[ii]);
                    }
                  else
                    {
                      // Slugs within domain.
                      assert(domain->bot_slug[ii]->bot <= domain->layer_bottom_depth);
                    }

                  // Check each slug.
                  slug* temp_slug = domain->top_slug[ii];

                  while (NULL != temp_slug)
                    {
                      // Slug top less than bottom.
                      assert(temp_slug->top < temp_slug->bot);

                      // Water to the left of slug.
                      if (1 < ii)
                        {
                          assert(has_water_at_depth(domain, ii - 1, temp_slug->top, temp_slug->bot));
                        }

                      if (NULL != temp_slug->next)
                        {
                          // Linked list consistency.
                          assert(temp_slug->next->prev == temp_slug);

                          // Slugs non-overlapping.
                          assert(temp_slug->bot < temp_slug->next->top);
                        }
                      else
                        {
                          // Linked list consistency.
                          assert(domain->bot_slug[ii] == temp_slug);
                        }

                      temp_slug = temp_slug->next;
                    } // End while (NULL != temp_slug).
                } // End check slugs.
            } // End the bin is not completely saturated.
        } // End process all bins.
    } // End if (NULL != domain).
#endif // NDEBUG
}

/* Comment in .h file. */
double t_o_total_water_in_domain(t_o_domain* domain)
{
  int    ii;          // Loop counter.
  double water = 0.0; // Accumulator for water in meters of bin depth.

  assert(NULL != domain);

  if (NULL != domain)
    {
      // Process all bins.
      for (ii = 1; ii <= domain->parameters->num_bins; ii++)
        {
          if (!domain->yes_groundwater && domain->parameters->bin_water_content[ii] <= domain->initial_water_content)
            {
              // The bin is completely saturated.
              water += domain->layer_bottom_depth - domain->layer_top_depth;
            }
          else // The bin is not completely saturated.
            {
              // Add surface front water.
              water += domain->surface_front[ii] - domain->layer_top_depth;

              // Add slugs.
              slug* temp_slug = domain->top_slug[ii];

              while (NULL != temp_slug)
                {
                  water += temp_slug->bot - temp_slug->top;
                  temp_slug = temp_slug->next;
                }

              // Add groundwater.
              if (domain->yes_groundwater)
                {
                  water += domain->layer_bottom_depth - domain->groundwater_front[ii];
                }
            } // End the bin is not completely saturated.
        } // end process all bins.
    } // End if (NULL != domain).

  // Multiply by delta_water_content to convert from meters of bin depth to meters of water and add residual saturation.
  return (water * domain->parameters->delta_water_content) +
      ((domain->layer_bottom_depth - domain->layer_top_depth) * (domain->parameters->bin_water_content[1] - domain->parameters->delta_water_content));
}

/* Create a slug struct and initialize it.
 * Return TRUE if there is an error, FALSE otherwise.
 * top and bot are initialized to the passed parameters.  prev and next are
 * initialized to NULL.  Rather than allocating a slug struct each time it is
 * called, this function might get it from a pool of unused slug structs.
 *
 * Parameters:
 *
 * new_slug - A pointer passed by reference which will be assigned to point to
 *            the newly allocated struct or NULL if there is an error.
 * top      - The depth in meters of the top    of the new slug.
 * bot      - The depth in meters of the bottom of the new slug.
 */
int slug_alloc(slug** new_slug, double top, double bot)
{
  int error = FALSE; // Error flag.

  //assert(NULL != new_slug && 0.0 <= top && top < bot);
    assert(NULL != new_slug && 0.0 <= top && top <= bot); // FIXME, WENCONG, change top < bot to <=, so that it can create an empty slug.
#ifdef THREAD_SAFE
  pthread_mutex_lock(&slug_pool_mutex);
#endif // THREAD_SAFE

  if (NULL != slug_pool)
    {
      // Instead of allocating, get a slug struct from the slug pool.
      *new_slug = slug_pool;
      slug_pool = slug_pool->next;
      
#ifdef THREAD_SAFE
      pthread_mutex_unlock(&slug_pool_mutex);
#endif // THREAD_SAFE
    }
  else
    {
#ifdef THREAD_SAFE
      pthread_mutex_unlock(&slug_pool_mutex); // Unlock before allocating to reduce contention.
#endif // THREAD_SAFE

      // Allocate a new slug struct.
      error = v_alloc((void**)new_slug, sizeof(slug));
    }

  if (!error)
    {
      (*new_slug)->prev = NULL;
      (*new_slug)->next = NULL;
      (*new_slug)->top  = top;
      (*new_slug)->bot  = bot;
    }

  return error;
}

/* Free memory allocated by slug_alloc.
 * Rather than freeing a slug struct each time it is called, this function
 * puts it back in a pool of unused slug structs.
 * 
 * Parameters:
 *
 * slug_to_kill - A pointer to the slug struct passed by reference.
 *                Will be set to NULL after the memory is deallocated.
 */
void slug_dealloc(slug** slug_to_kill)
{
  assert(NULL != slug_to_kill && NULL != *slug_to_kill);
  
  // Instead of deallocating, add the slug struct to the slug pool.
  // v_dealloc((void**)slug_to_kill, sizeof(slug));
#ifdef THREAD_SAFE
  pthread_mutex_lock(&slug_pool_mutex);
#endif // THREAD_SAFE
  
  (*slug_to_kill)->next = slug_pool;
	  slug_pool = *slug_to_kill;
  
#ifdef THREAD_SAFE
  pthread_mutex_unlock(&slug_pool_mutex);
#endif // THREAD_SAFE
  
  *slug_to_kill = NULL;
}

/* Create a new slug in domain in the given bin number.
 * Return TRUE if there is an error, FALSE otherwise.
 * If there is an error no slug is created.
 * The slug is placed after prev_slug in the linked list or at the top of the
 * linked list if prev_slug is NULL.  This function assumes prev_slug is in the
 * given bin number and the linked list position is correct for the given
 * values of top and bot.
 *
 * Parameters:
 *
 * domain    - A pointer to the t_o_domain struct.
 * bin       - Which bin to add the slug to.  One based indexing is used.
 * prev_slug - The new slug will be placed after this slug in the linked list.
 *             Pass in NULL to place the new slug at the top of the linked list
 *             or if the slug list for the bin is empty.
 * top       - The depth in meters of the top    of the new slug.
 * bot       - The depth in meters of the bottom of the new slug.
 */
int create_slug_after(t_o_domain* domain, int bin, slug* prev_slug, double top, double bot)
{
  int   error = FALSE; // Error flag.
  slug* new_slug;      // To point to the created slug.

  // FIXME, wencong.
//assert(NULL != domain && 0 < bin && bin <= domain->parameters->num_bins && domain->layer_top_depth <= top && top < bot && bot <= domain->layer_bottom_depth);
  assert(NULL != domain && 0 < bin && bin <= domain->parameters->num_bins && domain->layer_top_depth <= top && top <=bot && bot <= domain->layer_bottom_depth); 
  // FIXME, WENCONG, top <bot

  // Allocate the new slug.
  error = slug_alloc(&new_slug, top, bot);

  // Place it in the linked list.
  if (!error)
    {
      doubly_linked_list_insert_after((doubly_linked_list_element**)&domain->top_slug[bin], (doubly_linked_list_element**)&domain->bot_slug[bin],
                                      (doubly_linked_list_element*)prev_slug, (doubly_linked_list_element*)new_slug);
    }

  return error;
}

/* Remove a slug from its linked list and deallocate its memory.
 *
 * Parameters:
 *
 * domain       - A pointer to the t_o_domain struct.
 * bin          - Which bin to remove the slug from.
 *                One based indexing is used.
 * slug_to_kill - The slug to remove.
 */
void kill_slug(t_o_domain* domain, int bin, slug* slug_to_kill)
{
  assert(NULL != domain && 0 < bin && bin <= domain->parameters->num_bins && NULL != slug_to_kill);

  // Remove from the list.
  doubly_linked_list_remove_element((doubly_linked_list_element**)&domain->top_slug[bin], (doubly_linked_list_element**)&domain->bot_slug[bin],
                                    (doubly_linked_list_element*)slug_to_kill);

  // Deallocate.
  slug_dealloc(&slug_to_kill);
}

/* Change surface water in one bin into a slug.
 * Return TRUE if there is an error, FALSE otherwise.
 * If there is an error the slug is not detached.
 * The water is removed from the surface front and placed in a new slug with
 * its top at zero and bottom at the location of the surface front.
 * There is no vertical movement of the water.  No slug will be created
 * if the surface front is at zero and thus has no water.
 *
 * Parameters:
 *
 * domain - A pointer to the t_o_domain struct.
 * bin    - Which bin to detach.  One based indexing is used.
 */
int detach_slug(t_o_domain* domain, int bin)
{
  int error = FALSE; // Error flag.

  assert(NULL != domain && 0 < bin && bin <= domain->parameters->num_bins);

  if (domain->layer_top_depth < domain->surface_front[bin])
    {
      error = create_slug_after(domain, bin, NULL, domain->layer_top_depth, domain->surface_front[bin]);

      if (!error)
        {
          domain->surface_front[bin] = domain->layer_top_depth;
        }
    }

  return error;
}

/* Change surface water in all bins into slugs.
 * Return TRUE if there is an error, FALSE otherwise.
 * If there is an error some but not all of the slugs might be detached.
 * The water is removed from the surface front and placed in new slugs with
 * their top at zero and bottom at the location of the surface front.
 * There is no vertical movement of the water.  No slug will be created
 * for a bin if the surface front is at zero and thus has no water.
 *
 * Parameters:
 *
 * domain - A pointer to the t_o_domain struct.
 */
int detach_all_slugs(t_o_domain* domain)
{
  int error = FALSE; // Error flag.
  int ii;            // Loop counter.

  assert(NULL != domain);

  // Process all bins.
  for (ii = 1; !error && ii <= domain->parameters->num_bins; ii++)
    {
      error = detach_slug(domain, ii);
    }

  return error;
}

// FIXLATER possible optimization: Calculate first_bin and last_bin once
// and pass the values around and update them only if they change.
// FIXLATER possible optimization binary search instead of linear.

/* Return the leftmost bin that is not completely full of water or num_bins + 1
 * if all bins are completely full of water.  If yes_groundwater is FALSE then
 * bins to the right of initial_water_content are not considered completely
 * full of water even if surface_front reaches to layer_depth because they are
 * not in contact with groundwater.  first_bin is forced to be at least 2
 * because bin 1 should always be completely full of water, and we often use
 * first_bin - 1 as an array index.  If bin 1 is ever not completely full of
 * water then residual_saturation is set too high, and the simulation will be
 * innacurate because it needs more bins to the left.
 *
 * Parameters:
 *
 * domain       - A pointer to the t_o_domain struct.
 * start_search - Only search from this bin to the right.  Certain situations
 *                have a lower bound for first_bin.  Passing in this lower
 *                bound can speed things up.  If you don't know a lower bound
 *                pass in 2.
 */
int find_first_bin(t_o_domain* domain, int start_search)
{
  assert(NULL != domain && 2 <= start_search && start_search <= domain->parameters->num_bins + 1);

  int first_bin = start_search; // The leftmost bin that is not completely full of water.

  // We cannot use has_water_at_depth(0.0 to domain->layer_depth) because if yes_groundwater is FALSE and surface_front reaches to layer_depth then
  // has_water_at_depth will return TRUE even though that bin should not be considered in contact with groundwater.
  while(first_bin <= domain->parameters->num_bins &&
        ((domain->yes_groundwater && domain->layer_top_depth == domain->groundwater_front[first_bin]) ||
         (!domain->yes_groundwater && domain->parameters->bin_water_content[first_bin] <= domain->initial_water_content)))
    {
      first_bin++;
    }

  return first_bin;
}

/* Return the rightmost bin that has surface front water or 1 if no bins have
 * surface front water.  last_bin is forced to be at least 1 because bin 1
 * should always be completely full of water.  See the comment of
 * find_first_bin for why.
 *
 * Parameters:
 *
 * domain - A pointer to the t_o_domain struct.
 */
int find_last_bin(t_o_domain* domain)
{
  assert(NULL != domain);

  int last_bin = domain->parameters->num_bins; // The rightmost bin that has surface front water.
  
  while (1 < last_bin && !has_water_at_depth(domain, last_bin, domain->layer_top_depth, domain->layer_top_depth))
    {
      last_bin--;
    }

  return last_bin;
}

/* Process infiltration into the completely saturated bins.
 * Return TRUE if there is an error, FALSE otherwise.
 * If there is an error some but not all of the surface front water
 * might be detached as slugs, and not all of the demand might be
 * satisfied even if there is still surface front water and/or slugs.
 * This function takes water from the surface first.  If there is not enough
 * water on the surface, it detaches all surface front water into slugs and
 * then takes water from slugs.
 *
 * Parameters:
 *
 * domain               - A pointer to the t_o_domain struct.
 * dt                   - The duration of the timestep in seconds.
 * first_bin            - The leftmost bin that is not completely full of
 *                        water.
 * surfacewater_depth   - A scalar passed by reference containing the depth
 *                        in meters of the surface water.  Will be updated for
 *                        the amount of infiltration.
 * ponded_water         - A scalar passed by reference.  Will be set to whether
 *                        the demand of the fully saturated bins was completely
 *                        satisfied from surface water.
 * groundwater_recharge - A scalar passed by reference containing any
 *                        previously accumulated groundwater recharge in meters
 *                        of water.  Will be updated for the amount of water
 *                        that flowed between the Talbot-Ogden domain and
 *                        groundwater.  Positive means water flowed down out of
 *                        the Talbot-Ogden domain.  Negative means water flowed
 *                        up in to the Talbot-Ogden domain.
 * water_table          - The depth in meters of the water table. (add 06/17/14)
 */
int t_o_satisfy_saturated_bins(t_o_domain* domain, double dt, int first_bin, double* surfacewater_depth, int* ponded_water, double* groundwater_recharge,
                               double water_table)
{
  assert(NULL != domain && 0.0 < dt && NULL != surfacewater_depth && 0.0 <= *surfacewater_depth && NULL != ponded_water && NULL != groundwater_recharge);

  int    error = FALSE; // Error flag.
  int    ii, jj;        // Loop counters.
  double demand;        // Total water in meters of water that can infiltrate this timestep.
  double demand2;       // Demand of firstbin - 2.  We want the groundwater in saturated bins to fall only if they cannot be satisfied from surface water.
                        // However, if it is raining less than the saturated conductivity it may be able to satisfy only part of the last saturated bin.
                        // This can lead to an oscilation where it can't satisfy all of the saturated bins so the groundwater in all of the saturated bins
                        // falls including those bins that are being satisfied from surface water.  The rainfall can no longer infiltrate through the saturated
                        // bins and it becomes surface front water, which collides with the groundwater front increasing the number of saturated bins until the
                        // bin that can only be partially satisfied with rainfall becomes saturated and the cycle repeats.  To avoid this oscilation we allow
                        // the groundwater to fall only if the next to last saturated bin cannot be satisfied from surface water.  demand2 is the demand of the
                        // next to last saturated bin.
  
  // Calculate the demand of all the bins up to but not including first_bin, if not saturated.
  // FIXME do we want to vary the demand with surface water pressure head boundary condition?
  
   if (domain->yes_groundwater && domain->groundwater_front[domain->parameters->num_bins] <= domain->layer_top_depth)
    { // FIXME, wecnong, add 06/17/14, if fully saturated, calculated demand using Darcy's law.
      assert(water_table >= domain->layer_top_depth);
      demand = dt * domain->parameters->cumulative_conductivity[domain->parameters->num_bins] * (water_table - domain->layer_top_depth) / 
               (domain->layer_bottom_depth - domain->layer_top_depth);
      if (*surfacewater_depth   >= dt * domain->parameters->cumulative_conductivity[domain->parameters->num_bins])
        {
          *ponded_water          = TRUE;
        }
      else
        {
          *ponded_water          = FALSE;
          error = detach_all_slugs(domain);
        }
        
      if (*surfacewater_depth >= demand)
        {
          *surfacewater_depth   -= demand;
          *groundwater_recharge += demand;
        }
      else
        {
          *groundwater_recharge += *surfacewater_depth;
          *surfacewater_depth    = 0.0; 
        }
     
      return error;
    }
   else
    {
      demand = domain->parameters->cumulative_conductivity[first_bin - 1] * dt;
      if (first_bin > 2)
        {
          demand2 = domain->parameters->cumulative_conductivity[first_bin - 2] * dt;
        }
      else
        {
          demand2 = 0.95 * demand;
        }
    }// End of calculated demand of saturated bins.

  // Satisfy as much of the demand as possible from surface water.
  if (*surfacewater_depth   >= demand)
    {
      *ponded_water          = TRUE;
      *surfacewater_depth   -= demand;
      *groundwater_recharge += demand;
      demand                 = 0.0;
    }
  else
    {
      if (*surfacewater_depth >= demand2)
        {
          *ponded_water = TRUE;
        }
      else
        {
          *ponded_water = FALSE;

          // If the demand wasn't entirely satisfied from surface water, all surface front water will detach as slugs.
          error = detach_all_slugs(domain);
        }
      
      if (0.0 < *surfacewater_depth)
        {
          *groundwater_recharge += *surfacewater_depth;
          demand                -= *surfacewater_depth;
          *surfacewater_depth    = 0.0;
        }
    }
  
  if (!error)
    {
      // Satisfy the rest of the demand from slugs.
      int has_slugs = TRUE;

      while (epsilon_less(0.0, demand) && has_slugs)
        {
          // Get the rightmost slug.
          ii = domain->parameters->num_bins;

          while (ii > first_bin && NULL == domain->top_slug[ii])
            {
              ii--;
            }

          if (ii <= first_bin)
            {
              has_slugs = FALSE;
            }
          else
            {
              // Get the water from all slugs connected to domain->top_slug[ii].
              // FIXME get water from all slugs?

              double demand_save = demand; // We need to remember the value of the demand at the beginning of the loop iteration.
              slug*  get_slugs[ii + 1];    // Pointers to the slug in each bin that we will get water from.
              // Some elements might be NULL if we will not get water from a slug in that bin.
              double weight[ii + 1];       // Take from each slug a fraction of the demand equal to its normalized weight.
              double total_weight = 0.0;   // For normalizing weights.

              // Determine the weights.
              for (jj = ii; jj >= first_bin; jj--)
                {
                  // Find the slug.
                  get_slugs[jj] = domain->top_slug[jj];

                  while (NULL != get_slugs[jj] && domain->top_slug[ii]->top >= get_slugs[jj]->bot)
                    {
                      get_slugs[jj] = get_slugs[jj]->next;
                    }

                  if (NULL != get_slugs[jj] && !(domain->top_slug[ii]->top >= get_slugs[jj]->top && domain->top_slug[ii]->bot <= get_slugs[jj]->bot))
                    {
                      get_slugs[jj] = NULL;
                    }

                  // At this point in the code if get_slugs[jj] is NULL then there is no slug in that bin at that height.
                  // If it is not NULL, get water from it.

                  // Determine the weight of this slug.
                  if (NULL != get_slugs[jj])
                    {
                      weight[jj]    = (domain->parameters->bin_capillary_suction[first_bin - 1] - domain->parameters->bin_capillary_suction[jj]);
                      total_weight += weight[jj];
                    }
                }

              // Get the water.
              for (jj = ii; jj >= first_bin; jj--)
                {
                  if (NULL != get_slugs[jj])
                    {
                      double water = demand_save * weight[jj] / total_weight;         // Water to take from this bin in meters of water.
                      double depth = water / domain->parameters->delta_water_content; // Water to take from this bin in meters of bin depth.

                      if (depth >= get_slugs[jj]->bot - get_slugs[jj]->top)
                        {
                          // Get the entire slug.
                          water = (get_slugs[jj]->bot - get_slugs[jj]->top) * domain->parameters->delta_water_content;
                          kill_slug(domain, jj, get_slugs[jj]);
                        }
                      else
                        {
                          // FIXLATER implement complicated function to allocate water gotten to top and bottom?
                          // Get the water evenly from the top and bottom of the slug.
                          get_slugs[jj]->top += depth / 2.0;
                          get_slugs[jj]->bot -= depth / 2.0;
                        }

                      *groundwater_recharge += water;
                      demand                -= water;
                    }
                }
            } // End get the water from all slugs connected to domain->top_slug[ii].
        } // End while (epsilon_less(0.0, demand) && has_slugs)
    } // End if (!error)
   
  return error;
}

/*  Find dry depth for Green-Ampt equaiton.
 *   The funciton use to define dry depth is:
 *         Z * porosity = Ks * t + Hc * porosity * ln(1 + Z  / Hc) 
 *   Parameters:
 *   Z        - dry depth to be found in meter [m];
 *   Ks       - saturated conductiviey in meter/second [m/s];
 *   Hc       - effecitve capillary drive(Geff) in meter [m];
 *   dt       - time step in second [s];
 *   porosity - porosity [-].
 */
double GA_drydepth(double Ks, double porosity, double Geff, double dt)
{
  double minimum_dry_depth     = 1.0e-4; // Meters.
  double convergence_tolerance = 1.0e-6; // Meters.
  double iteration_difference;           // Meters.
  double z_old;                          // Meters.
  double z_new;                          // Meters.
  double f;
  double f_prime;
  int    iteration_count = 0;
  int    iteration_limit = 10000;        // maximum number of iterations before giving up.
  
  assert(0 < Ks && 0 < porosity && 0 < Geff && 0.0 < dt);
  
  z_old = Geff;
  
  do
    {
      f       = Ks * dt + Geff * porosity * log(1.0 + z_old / Geff) - z_old * porosity;
      f_prime = Geff * porosity / (Geff + z_old) - porosity;
      if (0.0 == f_prime)
        {
          // Prevent divide by zero.
          z_new = Ks * dt / porosity + Geff * log(1.0 + Ks * dt / (porosity * Geff));
          iteration_difference = 0.0;
        }
      else
        {
          z_new = z_old - (f / f_prime);
          iteration_difference = z_new - z_old;
        }
      
      z_old = z_new;
      iteration_count++;
    }
  while (fabs(iteration_difference) > convergence_tolerance && iteration_count < iteration_limit);

  if (iteration_count >= iteration_limit)
    {
      fprintf(stderr, "WARNING: No convergence in GA_drydepth.\n");
      z_new = Ks * dt / porosity + Geff * log(1.0 + Ks * dt / (porosity * Geff));
    }
  else if (z_new < minimum_dry_depth)
    {
      z_new = minimum_dry_depth;
    }

  return z_new;
}

/* There are numerical problems calculating the distance that water will
 * infiltrate into a bin that has little or no surface front water.
 * This function is used instead in that case.
 * Dry depth of each bin is calculated according to Han's method.
 * Dry depth is the solution z of following equation:
 * K * dt + psi * delth * ln (1 + z / psi ) - z * delth = 0
 * K     = conductivity [m/s]
 * dt    = time step    [second]
 * delth = delta water content [-]
 * psi   = bin capillary suction [m]
 * z     = dry depth [m]
 * Netwon-Raphson iteration is used for the solution.
 * F(z)_prime = delth / (1 + z / psi) - delth
 */
double t_o_find_dry_depth(t_o_parameters* parameters, int bin, double dt)
{
  double dry_depth;                      // Meters
  double minimum_dry_depth = 1.0e-4;     // Meters.
  double maximum_dry_depth;              // Meters.
  double convergence_tolerance = 1.0e-6; // Meters.
  double iteration_difference;           // Meters.
  double z_old;                          // Meters.
  double z_new;                          // Meters.
  double f;
  double f_prime;
  int    iteration_count = 0;
  int    iteration_limit = 10000;        // maximum number of iterations before giving up.

  assert(NULL != parameters && 0 < bin && bin <= parameters->num_bins && 0.0 < dt);

  z_old = parameters->cumulative_conductivity[bin] * dt;

  do
    {
      f       = (parameters->cumulative_conductivity[bin] * dt) + (parameters->bin_capillary_suction[bin] * parameters->delta_water_content *
          log(1.0 + (z_old / parameters->bin_capillary_suction[bin]))) - (z_old * parameters->delta_water_content);
      f_prime = (parameters->delta_water_content / (1.0 + (z_old / parameters->bin_capillary_suction[bin]))) - parameters->delta_water_content;

      if (0.0 == f_prime)
        {
          // Prevent divide by zero.
          z_new = minimum_dry_depth;
          iteration_difference = 0.0;
        }
      else
        {
          z_new = z_old - (f / f_prime);
          iteration_difference = z_new - z_old;
        }

      z_old = z_new;
      iteration_count++;
    }
  while (fabs(iteration_difference) > convergence_tolerance && iteration_count < iteration_limit);

  // FIXME, later. If Brooks_Corey,maximum_dry_depth = GA_dry_depth is fine, but if using van_Genutchen, maximum_dry_depth = 10 * GA_dry_depth.
  // Since there is no Flag to indicate what parameters is used, just time 10.
  maximum_dry_depth = 10 * GA_drydepth(parameters->cumulative_conductivity[parameters->num_bins], 
                                  parameters->bin_water_content[parameters->num_bins],
                                  parameters->effective_capillary_suction, dt);

  if (iteration_count >= iteration_limit)
    {
      fprintf(stderr, "WARNING: No convergence in t_o_find_dry_depth in bin %d.  Using a value of %le meters\n", bin, minimum_dry_depth);
      dry_depth = minimum_dry_depth;
    }
  else if (z_new < minimum_dry_depth)
    {
      dry_depth = minimum_dry_depth;
    }
  else if (z_new > maximum_dry_depth)
    {
      dry_depth = maximum_dry_depth;
    }
  else
    {
      dry_depth = z_new;
    }
  
  return dry_depth;
}

/* Calculate the distance in meters that water will infiltrate into all of the
 * bins in one timestep.
 *
 * Parameters:
 *
 * domain            - A pointer to the t_o_domain struct.
 * dt                - The duration of the timestep in seconds.
 * first_bin         - The leftmost bin that is not completely full of water.
 * surfacewater_head - The pressure head in meters of the surface water.
 * distance          - A 1D array sized to hold domain->parameters->num_bins
 *                     elements with one based indexing.  distance[ii] is
 *                     filled in with the distance that water will infiltrate
 *                     into bin ii this timestep.
 */
void infiltrate_distance(t_o_domain* domain, double dt, int first_bin, double surfacewater_head, double* distance)
{
  assert(NULL != domain && 0.0 < dt && NULL != distance);

  int ii;                               // Loop counter.
  int last_bin = find_last_bin(domain); // The rightmost bin that has surface front water.
  
  while (last_bin >= first_bin && 0 >= domain->parameters->bin_capillary_suction[last_bin] + surfacewater_head)
    {
      last_bin--;
    }

  // There are numerical problems calculating the distance that water will
  // infiltrate into a bin that has little or no surface front water.
  // In this case a different calculation called dry depth is used to
  // figure out the infiltration distance. The dry depth should be used
  // any time the surface front is less than the dry depth.  So you need
  // to know the dry depth to know if you need to use the dry depth.
  // Calculating dry depths is somewhat expensive so we have developed
  // a caching scheme to improve performance.
  //
  // Inside the t_o_parameters structure we cache the dry depth of every bin
  // for a particular timestep duration.  If the current timestep is equal to
  // the cached timestep we just use those values.  If the current timestep is
  // greater than the cached timestep we update the cached dry depths at the
  // larger timestep and then use those values.  If the current timestep is
  // less than the cached timestep we do not update the cached values.
  // Instead, we use the cached values as upper bounds to exclude some bins
  // that don't need to use dry depth and calculate the exact dry depth for
  // bins that we cannot exclude with the upper bound.  This will result in
  // values being cached for the largest timestep ever used in the simulation,
  // which is an upper bound for all other timesteps.
#ifdef THREAD_SAFE
  pthread_mutex_lock(&domain->parameters->dry_depth_mutex);
#endif // THREAD_SAFE
  
  if (domain->parameters->dry_depth_dt < dt)
    {
      domain->parameters->dry_depth_dt = dt;

      for (ii = 1; ii <= domain->parameters->num_bins; ii++)
        {
          domain->parameters->bin_dry_depth[ii] = t_o_find_dry_depth(domain->parameters, ii, dt);
        }
    }

  assert(domain->parameters->dry_depth_dt >= dt);

#ifdef THREAD_SAFE
  pthread_mutex_unlock(&domain->parameters->dry_depth_mutex);
#endif // THREAD_SAFE
  
  for (ii = first_bin; ii <= domain->parameters->num_bins; ii++)
    {
      double dry_depth;

#ifdef THREAD_SAFE
      pthread_mutex_lock(&domain->parameters->dry_depth_mutex);
#endif // THREAD_SAFE
  
      if (domain->parameters->dry_depth_dt > dt && domain->parameters->bin_dry_depth[ii] + domain->layer_top_depth >= domain->surface_front[ii])
        {
#ifdef THREAD_SAFE
          pthread_mutex_unlock(&domain->parameters->dry_depth_mutex); // Unlock before calling to reduce contention.
#endif // THREAD_SAFE
          
          // Cannot exclude bin based on upper bound.  Must calculate exact dry depth.
          dry_depth = t_o_find_dry_depth(domain->parameters, ii, dt);
        }
      else
        {
          dry_depth = domain->parameters->bin_dry_depth[ii];

#ifdef THREAD_SAFE
  pthread_mutex_unlock(&domain->parameters->dry_depth_mutex);
#endif // THREAD_SAFE
        }

      if (0 >= domain->parameters->bin_capillary_suction[ii] + surfacewater_head)
        {
          // If the surfacewater head has more suction than the bin capillarity set the distance to zero.
          distance[ii] = 0.0;
        }
      else if (dry_depth + domain->layer_top_depth >= domain->surface_front[ii])
        {
          // If there is downward infiltration and surface_front is less than dry_depth set distance to dry_depth.
          distance[ii] = dry_depth;
        }
      else
        {
          // Clip capillary suction with effective capillary suction.
          double last_bin_capillary_suction = domain->parameters->bin_capillary_suction[last_bin];

          if (last_bin_capillary_suction < domain->parameters->effective_capillary_suction)
            {
              last_bin_capillary_suction = domain->parameters->effective_capillary_suction;
            }

          // If last_bin is equal to first_bin - 1 then all bins will be set to zero or dry_depth and this equation will not be evaluated.
          distance[ii] = ((domain->parameters->cumulative_conductivity[last_bin] - domain->parameters->cumulative_conductivity[first_bin - 1]) /
              (domain->parameters->bin_water_content[last_bin] - domain->parameters->bin_water_content[first_bin - 1])) *
              ((last_bin_capillary_suction + surfacewater_head) / domain->surface_front[ii] + 1) * dt;
          
          // Runge-Kutta 4 implementation.
          /*
          double k1, k2, k3, k4, k0;
          k0 = ((domain->parameters->cumulative_conductivity[last_bin] - domain->parameters->cumulative_conductivity[first_bin - 1]) /
              (domain->parameters->bin_water_content[last_bin] - domain->parameters->bin_water_content[first_bin - 1]));
          k1 = k0 * ((last_bin_capillary_suction + surfacewater_head) /  domain->surface_front[ii] + 1.0) * dt;
          k2 = k0 * ((last_bin_capillary_suction + surfacewater_head) / (domain->surface_front[ii] + 0.5 * k1) + 1.0) * dt;
          k3 = k0 * ((last_bin_capillary_suction + surfacewater_head) / (domain->surface_front[ii] + 0.5 * k2) + 1.0) * dt;
          k4 = k0 * ((last_bin_capillary_suction + surfacewater_head) / (domain->surface_front[ii] + k3) + 1.0) * dt;
          distance[ii] = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
          */
          
          // 1-GARTO type.
          /*distance[ii] = (domain->parameters->cumulative_conductivity[ii] - domain->parameters->cumulative_conductivity[ii - 1]) /
              (domain->parameters->delta_water_content) *
              ((last_bin_capillary_suction + surfacewater_head) / domain->surface_front[ii] + 1) * dt;
          */
          /*
          // 2-exact k'
          double saturation  = (domain->parameters->bin_water_content[ii] - 
                                  (domain->parameters->bin_water_content[1] - domain->parameters->delta_water_content)) /
                                (domain->parameters->bin_water_content[domain->parameters->num_bins] - 
                                  (domain->parameters->bin_water_content[1] - domain->parameters->delta_water_content)); 
              
          distance[ii] = domain->parameters->cumulative_conductivity[domain->parameters->num_bins] * (3.0 + 2.0 / domain->parameters->bc_lambda) * 
                         pow(saturation, 2.0 + 2.0 / domain->parameters->bc_lambda) / (domain->parameters->bin_water_content[domain->parameters->num_bins] - 
                                  (domain->parameters->bin_water_content[1] - domain->parameters->delta_water_content)) * 
                           ((last_bin_capillary_suction + surfacewater_head) / domain->surface_front[ii] + 1) * dt;
            */  
        } // End if (dry_depth >= domain->surface_front[ii]).
    } // End for (ii = first_bin; ii <= domain->parameters->num_bins; ii++).
}

/* Process infiltration into not completely saturated bins.
 * Return TRUE if there is an error, FALSE otherwise.
 * Actually always returns FALSE.  No conditions generate an error.
 * This function takes water from the surface and then from surface front
 * water in bins to the right if there is not enough surface water.
 * If domain->yes_groundwater is FALSE then surface front water might
 * fall out the bottom of the domain.  This water is accounted for in
 * groundwater_recharge
 *
 * Parameters:
 *
 * domain               - A pointer to the t_o_domain struct.
 * dt                   - The duration of the timestep in seconds.
 * first_bin            - A scalar passed by reference containing the leftmost
 *                        bin that is not completely full of water.  Will be
 *                        updated to the new value after infiltration.
 * surfacewater_head    - The pressure head in meters of the surface water.
 *                        Normally, you will pass the same value as
 *                        surfacewater_depth, but this can be used to simulate
 *                        various laboratory test conditions.
 * surfacewater_depth   - A scalar passed by reference containing the depth
 *                        in meters of the surface water.  Will be updated for
 *                        the amount of infiltration.
 * groundwater_recharge - A scalar passed by reference containing any
 *                        previously accumulated groundwater recharge in meters
 *                        of water.  Will be updated for the amount of water
 *                        that flowed between the Talbot-Ogden domain and
 *                        groundwater.  Positive means water flowed down out of
 *                        the Talbot-Ogden domain.  Negative means water flowed
 *                        up in to the Talbot-Ogden domain.
 * ponded_water         - If ponded water is TRUE infiltrate even if
 *                        surfacewater_depth is zero.
 */
int t_o_infiltrate(t_o_domain* domain, double dt, int* first_bin, double surfacewater_head, double* surfacewater_depth, double* groundwater_recharge, 
                   int ponded_water)
{
  int error = FALSE; // Error flag.
  int ii;            // Loop counter.

  assert(NULL != domain && 0.0 < dt && NULL != surfacewater_depth && 0.0 <= *surfacewater_depth && NULL != groundwater_recharge);

  if (0.0 < *surfacewater_depth || ponded_water)
    {
      double delta_z[domain->parameters->num_bins + 1]; // The distance that water can infiltrate into each bin this timestep.

      // Calculate the depth that water can infiltrate into each bin.
      infiltrate_distance(domain, dt, *first_bin, surfacewater_head, delta_z);

      // All bins to the left of firstbin have already had their demand satisfied by satisfy_saturated_bins so start processing at first_bin.
      for (ii = *first_bin; ii <= domain->parameters->num_bins; ii++)
        {
          double supplied_z = 0.0; // Depth actually infiltrated.

          // Infiltration might be limited by the surface water hitting a slug or groundwater.
          int hit_slug        = FALSE;
          int hit_groundwater = FALSE;

          if (NULL != domain->top_slug[ii])
            {
              double gap = (domain->top_slug[ii]->top - domain->surface_front[ii]);

              if (delta_z[ii] >= gap)
                {
                  delta_z[ii]  = gap;
                  hit_slug     = TRUE;
                }
            }
          else if (domain->yes_groundwater)
            {
              double gap = (domain->groundwater_front[ii] - domain->surface_front[ii]);

              if (delta_z[ii] >= gap)
                {
                  delta_z[ii]  = gap;
                  hit_groundwater = TRUE;
                }
            }

          // Satisfy as much of the demand as possible from surface water.
          if (*surfacewater_depth >= delta_z[ii] * domain->parameters->delta_water_content)
            {
              *surfacewater_depth -= delta_z[ii] * domain->parameters->delta_water_content;
              supplied_z          += delta_z[ii];
              delta_z[ii]          = 0.0;
            }
          else if (0.0 < *surfacewater_depth)
            {
              supplied_z          += *surfacewater_depth / domain->parameters->delta_water_content;
              delta_z[ii]         -= *surfacewater_depth / domain->parameters->delta_water_content;
              *surfacewater_depth  = 0.0;
            }

          // Satisfy the rest of the demand from the rightmost bin that has surface front water.

          int get_bin = domain->parameters->num_bins; // The bin to get water from.

          while (0.0 < delta_z[ii] && get_bin > ii)
            {
              if (domain->surface_front[get_bin] - domain->layer_top_depth >= delta_z[ii])
                {
                  // The bin has enough to completely satisfy remaining demand.
                  domain->surface_front[get_bin] -= delta_z[ii];
                  supplied_z                     += delta_z[ii];
                  delta_z[ii]                     = 0.0;
                }
              else
                {
                  // Get everything the bin has, if any, and go on to the next bin.
                  if (0.0 < domain->surface_front[get_bin] - domain->layer_top_depth)
                    {   
                      supplied_z                     += domain->surface_front[get_bin] - domain->layer_top_depth;
                      delta_z[ii]                    -= domain->surface_front[get_bin] - domain->layer_top_depth;
                      domain->surface_front[get_bin]  = domain->layer_top_depth;
                    }
                 
                  get_bin--;
                }
            }

          // Add the water to the surface front water in this bin.
          if (hit_slug && 0.0 == delta_z[ii])
            {
              // Surface water reaches the top slug.
              domain->surface_front[ii] = domain->top_slug[ii]->bot;
              kill_slug(domain, ii, domain->top_slug[ii]);
            }
          else if (hit_groundwater && 0.0 == delta_z[ii])
            {
              // Surface water reaches groundwater.
              domain->surface_front[ii]     = domain->layer_top_depth;
              domain->groundwater_front[ii] = domain->layer_top_depth;
            }
          else if (domain->surface_front[ii] + supplied_z > domain->layer_bottom_depth)
            {
              // Surface water reaches the bottom of the domain.
              *groundwater_recharge += (domain->surface_front[ii] + supplied_z - domain->layer_bottom_depth) * domain->parameters->delta_water_content;
              domain->surface_front[ii] = domain->layer_bottom_depth;
            }
          else
            {
              // Advance surface_front.
              domain->surface_front[ii] += supplied_z;
            }
        } // End loop over all bins starting at first_bin
      
      // first_bin can only change if there was infiltration, and it can only move to the right.
      *first_bin = find_first_bin(domain, *first_bin);
    } // End if (0.0 < *surfacewater_depth)

  return error;
}

/* Calculate the distance in meters that a slug will move in one timestep.
 * This distance can be different for the top and bottom of the slug.
 * Rather than returning one value it uses two output parameters.
 *
 * Parameters:
 *
 * domain       - A pointer to the t_o_domain struct.
 * bin          - Which bin falling_slug is in.
 * first_bin    - The leftmost bin that is not completely full of water.
 * dt           - The duration of the timestep in seconds.
 * falling_slug - Which slug to determine the distance for.
 * top_distance - A scalar passed by reference that gets set to the distance
 *                in meters that the top of the slug will move.
 *                Positive means down toward the bottom of the domain.
 * bot_distance - A scalar passed by reference that gets set to the distance
 *                in meters that the bottom of the slug will move.
 *                Positive means down toward the bottom of the domain.
 */
void slug_fall_distance(t_o_domain* domain, int bin, int first_bin, double dt, slug* falling_slug, double* top_distance, double* bot_distance)
{
  int   ii;        // Loop counter.
  slug* temp_slug; // Used to search for connected slugs.

  assert(NULL != domain && 0 < bin && bin <= domain->parameters->num_bins && 0 < first_bin && first_bin <= domain->parameters->num_bins &&
      0.0 < dt && NULL != falling_slug && NULL != top_distance && NULL != bot_distance);

  // find the capillary suction of the rightmost bin with a connected slug that this slug can steal water from.
  ii = domain->parameters->num_bins;

  do
    {
      temp_slug = domain->top_slug[ii];

      while (NULL != temp_slug && !(falling_slug->top <= temp_slug->bot && falling_slug->bot >= temp_slug->top))
        {
          if (falling_slug->bot < temp_slug->top)
            {
              temp_slug = NULL;
            }
          else
            {
              temp_slug = temp_slug->next;
            }
        }

      if (NULL == temp_slug)
        {
          ii--;
        }
    }
  while (ii > bin && NULL == temp_slug);

  double distance = ((domain->parameters->cumulative_conductivity[bin] - domain->parameters->cumulative_conductivity[bin - 1]) /
      domain->parameters->delta_water_content) * dt;

  double growth;

  if (ii <= bin)
    {
      growth = 0.0;
    }
  else
    {
      growth = ((domain->parameters->cumulative_conductivity[bin] - domain->parameters->cumulative_conductivity[first_bin - 1]) /
          (domain->parameters->bin_water_content[bin] - domain->parameters->bin_water_content[first_bin - 1])) *
              ((domain->parameters->bin_capillary_suction[bin] - domain->parameters->bin_capillary_suction[ii]) /
                  ((falling_slug->bot - falling_slug->top) / 2.0)) * dt;
    }

  *top_distance = distance - growth;
  *bot_distance = distance + growth;
}

/* Process the falling slugs step of the simulation.
 * Return TRUE if there is an error, FALSE otherwise.
 * Actually always return FALSE.  No conditions generate an error.
 * If domain->groundwater_yes is FALSE then slugs might fall out the bottom of
 * the domain.  This water is accounted for in groundwater_recharge
 *
 * Parameters:
 *
 * domain               - A pointer to the t_o_domain struct.
 * dt                   - The duration of the timestep in seconds.
 * first_bin            - The leftmost bin that is not completely full of
 *                        water.
 * groundwater_recharge - A scalar passed by reference containing any
 *                        previously accumulated groundwater recharge in meters
 *                        of water.  Will be updated for the amount of water
 *                        that flowed between the Talbot-Ogden domain and
 *                        groundwater.  Positive means water flowed down out of
 *                        the Talbot-Ogden domain.  Negative means water flowed
 *                        up in to the Talbot-Ogden domain.
 */
int t_o_falling_slugs(t_o_domain* domain, double dt, int first_bin, double* groundwater_recharge)
{
  assert(NULL != domain && 0.0 < dt && NULL != groundwater_recharge);

  int error = FALSE; // Error flag.
  int ii;            // Loop counter.

  // Process all bins except bin 1, which is guaranteed to be completely full of water and thus have no slugs.
  for (ii = 2; ii <= domain->parameters->num_bins; ii++)
    {
      // Process the slugs from bottom up.
      slug* temp_slug = domain->bot_slug[ii];

      while (NULL != temp_slug)
        {
          double top_delta_z; // The distance the top    of the slug will move.
          double bot_delta_z; // The distance the bottom of the slug will move.

          // Calculate the distance the slug will move this timestep.
          slug_fall_distance(domain, ii, first_bin, dt, temp_slug, &top_delta_z, &bot_delta_z);

          // FIXME deal with this
          // Prevent the slug from moving upward.
          if (0.0 >= top_delta_z)
            {
              top_delta_z = 0.000001;
            }

          if (0.0 >= bot_delta_z)
            {
              bot_delta_z = 0.000001;
            }

          // Stop the slug if it hits a lower slug or groundwater.
          if (NULL != temp_slug->next)
            {
              if (bot_delta_z > temp_slug->next->top - temp_slug->bot)
                {
                  bot_delta_z = temp_slug->next->top - temp_slug->bot;
                }
            }
          else if (domain->yes_groundwater)
            {
              if (bot_delta_z > domain->groundwater_front[ii] - temp_slug->bot)
                {
                  bot_delta_z = domain->groundwater_front[ii] - temp_slug->bot;
                }
            }

          // Prevent the slug from shrinking.
          if (top_delta_z > bot_delta_z)
            {
              top_delta_z = bot_delta_z;
            }

          // If the slug grows, steal the water from the rightmost overlapping slug, and the topmost if there are multiple rightmost.
          // FIXME steal from all connected slugs to the right with weights.
          double demand   = bot_delta_z - top_delta_z;    // Needed water in meters of bin depth.
          int    get_bin  = domain->parameters->num_bins; // Which bin  to try to get from.
          slug*  get_slug = domain->top_slug[get_bin];    // Which slug to try to get from.

          while (0.0 < demand && ii < get_bin)
            {
              while (0.0 < demand && NULL != get_slug && temp_slug->bot >= get_slug->top)
                {
                  // Save a pointer to the slug below get_slug in case we need to kill get_slug.
                  slug* next_slug = get_slug->next;

                  if (temp_slug->top <= get_slug->bot && temp_slug->bot >= get_slug->top)
                    {
                      // We can get water from get_slug.
                      if (get_slug->bot - get_slug->top >= demand)
                        {
                          // FIXLATER more complicated than taking equally from top and bot?
                          get_slug->top += demand / 2.0;
                          get_slug->bot -= demand / 2.0;
                          demand = 0.0;
                        }
                      else
                        {
                          demand -= get_slug->bot - get_slug->top;
                          kill_slug(domain, get_bin, get_slug);
                        }
                    }

                  get_slug = next_slug;
                }

              get_bin--;
              get_slug = domain->top_slug[get_bin];
            }

          // If we couldn't get all of the desired water take the deficit equally from the top and bottom of where the slug wants to be.
          // FIXLATER more complicated than taking equally from top and bot?
          if (0.0 < demand)
            {
              top_delta_z += demand / 2.0;
              bot_delta_z -= demand / 2.0;
            }

          // Save a pointer to the slug above temp_slug in case we need to kill temp_slug.
          slug* prev_slug = temp_slug->prev;

          if (NULL == temp_slug->next)
            {
              // The bottom slug might hit the groundwater or fall beyond layer_depth.
              if (domain->yes_groundwater)
                {
                  if (temp_slug->bot + bot_delta_z >= domain->groundwater_front[ii])
                    {
                      // The slug hits groundwater.
                      domain->groundwater_front[ii] -= (temp_slug->bot + bot_delta_z) - (temp_slug->top + top_delta_z);
                      kill_slug(domain, ii, temp_slug);
                    }
                  else
                    {
                      // Advance the slug.
                      temp_slug->top += top_delta_z;
                      temp_slug->bot += bot_delta_z;
                    }
                }
              else
                {
                  if (temp_slug->top + top_delta_z >= domain->layer_bottom_depth)
                    {
                      // The slug falls entirely past layer_depth.
                      *groundwater_recharge += ((temp_slug->bot + bot_delta_z) - (temp_slug->top + top_delta_z)) * domain->parameters->delta_water_content;
                      kill_slug(domain, ii, temp_slug);
                    }
                  else if (temp_slug->bot + bot_delta_z > domain->layer_bottom_depth)
                    {
                      // The slug falls partially past layer depth.
                      *groundwater_recharge += (temp_slug->bot + bot_delta_z - domain->layer_bottom_depth) * domain->parameters->delta_water_content;
                      temp_slug->top        += top_delta_z;
                      temp_slug->bot         = domain->layer_bottom_depth;
                    }
                  else
                    {
                      // Advance the slug.
                      temp_slug->top += top_delta_z;
                      temp_slug->bot += bot_delta_z;
                    }
                }
            }
          else
            {
              // A middle slug might hit the slug below it.
              if (temp_slug->bot + bot_delta_z >= temp_slug->next->top)
                {
                  // The slug hits the slug below it.
                  temp_slug->next->top -= (temp_slug->bot + bot_delta_z) - (temp_slug->top + top_delta_z);
                  kill_slug(domain, ii, temp_slug);
                }
              else
                {
                  // Advance the slug.
                  temp_slug->top += top_delta_z;
                  temp_slug->bot += bot_delta_z;
                }
            }

          temp_slug = prev_slug;
        } // End while (NULL != temp_slug).
    } // End for (ii = 2; ii <= domain->parameters->num_bins; ii++).

  return error;
}

/* Return the distance in meters that groundwater in a bin will move in one
 * timestep.  Positive means down toward the bottom of the domain.
 *
 * Parameters:
 *
 * domain      - A pointer to the t_o_domain struct.
 * bin         - Which bin to calculate distance for.
 * first_bin   - The leftmost bin that is not completely full of water.
 * dt          - The duration of the timestep in seconds.
 * water_table - The depth in meters of the water table.
 * inflow_rate - Flow rate through fully saturated bins in meters per second.
 */
double groundwater_distance(t_o_domain* domain, int bin, int first_bin, double dt, double water_table, double inflow_rate)
{
  double distance = 0.0;

  assert(NULL != domain && 0 < bin && bin <= domain->parameters->num_bins && 0 < first_bin && first_bin <= domain->parameters->num_bins &&
      0.0 < dt && 0.0 <= water_table);

  if (domain->yes_groundwater)
    {  
      // Original code;
      /*
      if (0.0 >= water_table && 0.0 == domain->groundwater_front[bin])
        {
          distance = 0.0;
        }
      else if (water_table <= domain->groundwater_front[bin])
        {
          // If groundwater falls below the water table there are numerical problems.  Move the groundwater to 10% of its capillary head above the water table.
          fprintf(stderr, "WARNING: Groundwater below water table moved up to 10 percent of capillary head above water table.\n");
          distance = (water_table - 0.1 * domain->parameters->bin_capillary_suction[bin]) - domain->groundwater_front[bin];
        }
      else
        {
          distance = ((domain->parameters->cumulative_conductivity[bin] - domain->parameters->cumulative_conductivity[first_bin - 1]) /
              (domain->parameters->bin_water_content[bin] - domain->parameters->bin_water_content[first_bin - 1])) *
                  ((-domain->parameters->bin_capillary_suction[bin] / (water_table - domain->groundwater_front[bin])) + 1) * dt;

          // Do not allow the groundwater to travel beyond hydrostatic.
          double distance_to_hydrostatic = (water_table - domain->parameters->bin_capillary_suction[bin]) - domain->groundwater_front[bin];

          if ((0.0 >= distance_to_hydrostatic && distance < distance_to_hydrostatic) || (0.0 <= distance_to_hydrostatic && distance > distance_to_hydrostatic))
            {
              distance = distance_to_hydrostatic;
            }
        }
      */
      
      // Hydrostatic suction considering steady rainfall rate. Add 06/02/14.
      if (domain->layer_top_depth >= water_table && domain->layer_top_depth == domain->groundwater_front[bin])
        {
          distance = 0.0;
        }
      else if (water_table <= domain->groundwater_front[bin])
        {
          // If groundwater falls below the water table there are numerical problems.  Move the groundwater to 10% of its capillary head above the water table.
          fprintf(stderr, "WARNING: Groundwater below water table moved up to 10 percent of capillary head above water table.\n");
          distance = (water_table - 0.1 * domain->parameters->bin_capillary_suction[bin]) - domain->groundwater_front[bin];
          printf("bin = %d water_table = %lf, domain->groundwater_front[bin] = %lf \n", bin, water_table, domain->groundwater_front[bin]); getchar();
          // FIXME, remove me, Feb, 09, 2015.
        }
      else
        {
          double suction = domain->parameters->bin_capillary_suction[bin];
          if (domain->parameters->cumulative_conductivity[1] < inflow_rate && inflow_rate < domain->parameters->cumulative_conductivity[bin]) 
            { // Calculate hydrostatic suction considering steady rainfall. 
              double suction_new = domain->parameters->bc_psib / (1.0 - inflow_rate / domain->parameters->cumulative_conductivity[domain->parameters->num_bins]) 
                                   + (domain->parameters->bin_capillary_suction[bin] - domain->parameters->bc_psib) / 
                                     (1.0 - 0.5 * inflow_rate / domain->parameters->cumulative_conductivity[bin]
                                          - 0.5 * inflow_rate / domain->parameters->cumulative_conductivity[domain->parameters->num_bins]);

              if (water_table - domain->layer_top_depth > suction_new || 
                  water_table - domain->layer_top_depth < domain->parameters->bin_capillary_suction[domain->parameters->num_bins])
                {
                  suction = suction_new;
                }
              else if (water_table > domain->layer_top_depth)
                {
                  suction = 0.99 * water_table;
                }
            }
          
          distance = ((domain->parameters->cumulative_conductivity[bin] - domain->parameters->cumulative_conductivity[first_bin - 1]) /
              (domain->parameters->bin_water_content[bin] - domain->parameters->bin_water_content[first_bin - 1])) *
                  ((-suction / (water_table - domain->groundwater_front[bin])) + 1) * dt;

          // Do not allow the groundwater to travel beyond hydrostatic.
          //double distance_to_hydrostatic = (water_table - domain->parameters->bin_capillary_suction[bin]) - domain->groundwater_front[bin];
          double distance_to_hydrostatic = (water_table - suction) - domain->groundwater_front[bin];

          if ((0.0 >= distance_to_hydrostatic && distance < distance_to_hydrostatic) || (0.0 <= distance_to_hydrostatic && distance > distance_to_hydrostatic))
            {
              distance = distance_to_hydrostatic;
            }
        }
    }

  return distance;
}

/* Process the groundwater step of the simulation.
 * Return TRUE if there is an error, FALSE otherwise.
 * Actually always returns FALSE.  No conditions generate an error.
 *
 * Parameters:
 *
 * domain               - A pointer to the t_o_domain struct.
 * dt                   - The duration of the timestep in seconds.
 * first_bin            - A scalar passed by reference containing the leftmost
 *                        bin that is not completely full of water.  Will be
 *                        updated to the new value after infiltration.
 * water_table          - The depth in meters of the water table.
 * ponded_water         - Whether the demand of the fully saturated bins was
 *                        completely satisfied from surface water.
 *                        Pass in the value set by satisfy_saturated_bins.
 *                        If this is TRUE then bins less than first_bin are
 *                        left pinned to the surface.  Otherwise, their
 *                        groundwater front is simulated the same as other
 *                        bins.
 * groundwater_recharge - A scalar passed by reference containing any
 *                        previously accumulated groundwater recharge in meters
 *                        of water.  Will be updated for the amount of water
 *                        that flowed between the Talbot-Ogden domain and
 *                        groundwater.  Positive means water flowed down out of
 *                        the Talbot-Ogden domain.  Negative means water flowed
 *                        up in to the Talbot-Ogden domain.
 * inflow_rate          - Flow rate through fully saturated bins in meters per second.
 */
int t_o_groundwater(t_o_domain* domain, double dt, int* first_bin, double water_table, int ponded_water, double* groundwater_recharge, double inflow_rate)
{
  int error = FALSE; // Error flag.
  int ii;            // Loop counter.

  assert(NULL != domain && 0.0 < dt && 0.0 <= water_table && NULL != groundwater_recharge);

  if (domain->yes_groundwater)
    {
      // If the demand of the fully saturated bins was completely satisfied from surface water then bins less than first_bin are pinned to the surface.
      // Otherwise, reduce first_bin until you reach a bin whoose hydrostatic equlibrium is at or above the surface.
      if (!ponded_water)
        {
          while (2 < *first_bin && domain->layer_top_depth < water_table - domain->parameters->bin_capillary_suction[*first_bin - 1])
            {
              (*first_bin)--;
            }
        }

      for (ii = *first_bin; ii <= domain->parameters->num_bins; ii++)
        {
          // FIXLATER if we want to update first_bin in this loop we will need to modify groundwater_distance to produce an array of distances for all bins
          // like infiltrate_distance does so that we can calculate them all before first_bin changes
          // FIXME, wencong 6/2/14, add inflow_rate to calculate groundwater distance, as inflow rate affects hydrostatic capillary height.
          double delta_z = groundwater_distance(domain, ii, *first_bin, dt, water_table, inflow_rate); // The distance groundwater wants to move this timestep.

          // Move the water.
          if (0.0 > delta_z)
            {
              // Groundwater moves upward.
              
              /* FIXME do we really want to do this?  Probably use this in real world case, but not in tests with arbitrary water table location.
              // Don't make a positive recharge negative, if recharge is already negative, do nothing.
              if (*groundwater_recharge <= 0.0)
                {
                  delta_z = 0.0;
                }
              else if (delta_z < -(*groundwater_recharge) / domain->parameters->delta_water_content)
                {
                  delta_z = -(*groundwater_recharge) / domain->parameters->delta_water_content;
                }
              */
            
              double final_depth = domain->groundwater_front[ii] + delta_z; // Depth in meters that groundwater wants to move to. 
                 
              // Groundwater cannot move above the surface.
              if (final_depth < domain->layer_top_depth)
                {
                  final_depth = domain->layer_top_depth;
                }

              // Accumulate the amount of water actually moved in meters of bin depth.
              delta_z = 0.0;

              // We must loop because groundwater might overtake multiple slugs or slugs and the surface front in one timestep.
              while (domain->groundwater_front[ii] > final_depth)
                {
                  if (NULL != domain->bot_slug[ii])
                    {
                      // A slug exists.
                      if (final_depth <= domain->bot_slug[ii]->bot)
                        {
                          // The groundwater hits the bottom slug.
                          delta_z += domain->bot_slug[ii]->bot - domain->groundwater_front[ii];
                          domain->groundwater_front[ii] = domain->bot_slug[ii]->top;
                          kill_slug(domain, ii, domain->bot_slug[ii]);
                        }
                      else
                        {
                          // The groundwater does not hit the bottom slug.
                          delta_z += final_depth - domain->groundwater_front[ii];
                          domain->groundwater_front[ii] = final_depth;
                        }
                    }
                  else
                    {
                      // No slug exists.
                      if (final_depth <= domain->surface_front[ii])
                        {
                          // The groundwater hits the surface front.
                          delta_z += domain->surface_front[ii] - domain->groundwater_front[ii];
                          domain->surface_front[ii]     = domain->layer_top_depth;
                          domain->groundwater_front[ii] = domain->layer_top_depth;
                        }
                      else
                        {
                          // The groundwater does not hit the surface front.
                          delta_z += final_depth - domain->groundwater_front[ii];
                          domain->groundwater_front[ii] = final_depth;
                        }
                    }
                }
            }
          else if (0.0 < delta_z)
            {
              // Groundwater moves downward.
              if (domain->groundwater_front[ii] + delta_z >= domain->layer_bottom_depth)
                {
                  // The groundwater hits the bottom of the domain.
                  delta_z = domain->layer_bottom_depth - domain->groundwater_front[ii];
                  domain->groundwater_front[ii] = domain->layer_bottom_depth;
                }
              else
                {
                  // The groundwater does not hit the bottom of the domain.
                  domain->groundwater_front[ii] += delta_z;
                }
            }

          *groundwater_recharge += delta_z * domain->parameters->delta_water_content;
        } // End loop over all bins

      // Force bin 1 to be completely full of water.
      if (domain->layer_top_depth < domain->groundwater_front[1])
        {
          fprintf(stderr, "WARNING: Groundwater in bin 1 wants to fall below the surface.  Groundwater in bin 1 is being pinned "
              "to the surface.  The simulation will be inaccurate unless you decrease residual_saturation.\n");
          *groundwater_recharge += -(domain->groundwater_front[1] - domain->layer_top_depth) * domain->parameters->delta_water_content;
          domain->groundwater_front[1] = domain->layer_top_depth;
        }
      
      // Find the new value of first_bin.  It could move to the right or left.
      *first_bin = find_first_bin(domain, 2);
    } // End if (domain->yes_groundwater)

  return error;
}

//Function passed to qsort() to sort the list from largest to smallest
int
compare_surface(const double *a, const double *b)
{
  double temp = *a - *b;
  if (temp > 0)
    return -1;
  else if (temp < 0)
    return 1;
  else
    return 0;
}

//Function passed to qsort() to sort the list from smallest to largest
int
compare_ground(const double *a, const double *b)
{
  double temp = *a - *b;
  if (temp > 0)
    return 1;
  else if (temp < 0)
    return -1;
  else
    return 0;
}

int
insert_into_list_top_sort(slug* (*list), slug* (*end), slug* (*element))
{
  if ((*list) == NULL )
    {
      (*element)->next = NULL;
      (*element)->prev = NULL;
      (*list) = (*element);
      (*end) = (*list);
    }
  else if ((*list)->top >= (*element)->top)
    {
      (*element)->prev = NULL;
      (*element)->next = (*list);
      (*list)->prev = (*element);
      (*list) = (*element);
    }
  else
    {
      slug* current = (*list);
      while (current->next != NULL && current->next->top < (*element)->top)
        {
          current = current->next;
        }
      //found correct insertion point
      if (current->next != NULL )
        {
          (*element)->next = current->next;
          current->next->prev = (*element);
          current->next = (*element);
          (*element)->prev = current;
        }
      else
        {
          current->next = (*element);
          (*element)->prev = current;
          (*element)->next = NULL;
          (*end) = (*element);
        }
    }
  return 0;
}

int
insert_into_list_bot_sort(slug* (*list), slug* (*end), slug* (*element))
{
  if ((*list) == NULL )
    {
      (*element)->next = NULL;
      (*element)->prev = NULL;
      (*list) = (*element);
      (*end) = (*list);
    }
  else if ((*list)->bot <= (*element)->bot)
    {
      (*element)->prev = NULL;
      (*element)->next = (*list);
      (*list)->prev = (*element);
      (*list) = (*element);
    }
  else
    {
      slug* current = (*list);
      while (current->next != NULL && current->next->bot > (*element)->bot)
        {
          current = current->next;
        }
      //found correct insertion point
      if (current->next != NULL )
        {
          (*element)->next = current->next;
          current->next->prev = (*element);
          current->next = (*element);
          (*element)->prev = current;
        }
      else
        {
          current->next = (*element);
          (*element)->prev = current;
          (*element)->next = NULL;
          (*end) = (*element);
        }
    }
  return 0;
}

int
cut_slugs(t_o_domain* domain, slug* (*all), slug* (*top_list),
    slug* (*top_list_end), slug* (*bot_list), slug* (*bot_list_end),
    slug* (*mid_list), slug* (*mid_list_end), int first_bin)
{
  int error = FALSE;
  slug* head = (*all);
  slug* next;
  while (head != NULL )
    {
      next = head->next;
      //TODO TRY TO COMBINE SOME OF THESE CASES???
      if (head->bot <= domain->surface_front[first_bin])
        { //this entire slug goes into the top list
          insert_into_list_top_sort(&(*top_list), &(*top_list_end), &head);
        }
      else if(domain->yes_groundwater && head->top >= domain->groundwater_front[first_bin])
        { //this entire slug goes into the bottom list
          insert_into_list_bot_sort(&(*bot_list), &(*bot_list_end), &head);
        }
      else if((!domain->yes_groundwater || head->bot <= domain->groundwater_front[first_bin]) && head->top >= domain->surface_front[first_bin])
        { //this entire slug goes into the middle list
          insert_into_list_top_sort(&(*mid_list), &(*mid_list_end), &head);
        }
      else
        {
          //we have to split the slug up, three cases: top/mid, mid/bot, top/mid/bot
          if(head->top < domain->surface_front[first_bin] && (!domain->yes_groundwater || head->bot <= domain->groundwater_front[first_bin]))
            {
              //hit the top/mid case, create one new slug
              slug* sl;
              slug_alloc(&sl, domain->surface_front[first_bin], head->bot);
              head->bot = sl->top;
              //put head into top list
              insert_into_list_top_sort(&(*top_list), &(*top_list_end), &head);
              //put new slug, sl, into middle list
              insert_into_list_top_sort(&(*mid_list), &(*mid_list_end), &sl);
            }
          else if((domain->yes_groundwater && head->bot > domain->groundwater_front[first_bin]) && head->top >= domain->surface_front[first_bin])
            {
              //hit the mid/bot case, create one new slug
              slug* sl;
              if(slug_alloc(&sl, head->top, domain->groundwater_front[first_bin]))
                {
                  error = TRUE;
                  break;
                }
              head->top = sl->bot;
              //put head into bot list
              insert_into_list_bot_sort(&(*bot_list), &(*bot_list_end), &head);
              //put new slug, sl, into middle list
              insert_into_list_top_sort(&(*mid_list), &(*mid_list_end), &sl);
            }
          else if(head->top < domain->surface_front[first_bin] && (domain->yes_groundwater && head->bot > domain->groundwater_front[first_bin]))
            {
              //hit the top/mid/bot case, create two new slugs
              slug* new_t;
              slug* new_b;
              slug_alloc(&new_t, head->top, domain->surface_front[first_bin]);
              slug_alloc(&new_b, domain->groundwater_front[first_bin], head->bot);
              head->top = new_t->bot;
              head->bot = new_b->top;

              //put new_t into top list
              insert_into_list_top_sort(&(*top_list), &(*top_list_end), &new_t);
              //put head into middle list
              insert_into_list_top_sort(&(*mid_list), &(*mid_list_end), &head);
              //put new_b into bottom list
              insert_into_list_bot_sort(&(*bot_list), &(*bot_list_end), &new_b);
            }
          else
            {
              assert(FALSE); //should never get to this case!
            }
        }
      head = next;
    }
  return error;
}

/* Helper function to deal with merging of ground and surface
 * fronts
 */
int
find_collisions(t_o_domain* domain, slug* (*top_list), slug* (*top_list_end),
    slug* (*bot_list), slug* (*bot_list_end), slug* (*mid_list),
    slug* (*mid_list_end), int *first_bin)
{
  int error = FALSE;
  int i;
  int new_first_bin = *first_bin;
  //Find the first non overlapping bin to pass to cut slugs
  for(i = *first_bin; i <= domain->parameters->num_bins; i++)
    {
      if(domain->surface_front[i] < domain->groundwater_front[i])
        {
          new_first_bin = i;
          break;
        }
    }

  for (i = *first_bin; i <= domain->parameters->num_bins; i++)
    {
      if (domain->surface_front[i] > domain->groundwater_front[i])
        {
          //We have overlapping water, create slug
          slug* sl;
          if(slug_alloc(&sl, domain->groundwater_front[i], domain->surface_front[i]))
            {
              error = TRUE;
              break;
            }
          cut_slugs(domain, &sl, &(*top_list),
              &(*top_list_end), &(*bot_list), &(*bot_list_end), &(*mid_list),
              &(*mid_list_end), new_first_bin);

          //SLUGS ARE CREATED, UPDATE FRONTS TO SHOW FULL BIN
          domain->surface_front[i]     = domain->layer_top_depth;
          domain->groundwater_front[i] = domain->layer_top_depth;
        }
      else if (domain->surface_front[i] == domain->groundwater_front[i])
        {
          //if the ground and surface are exactly equal, no slug is created,
          //but need to "fill" the bin
          domain->surface_front[i]     = domain->layer_top_depth;
          domain->groundwater_front[i] = domain->layer_top_depth;
        }
      else
        {
          break;
        }
    }
  *first_bin = new_first_bin;
  return error;
}

/*
 * helper function to add a slug to the bin, checks for equality
 * and adds to ground/surface if necessary
 *
 * This function assumes that bin_slug will fit without overlapping
 * between two fronts, or between a single slug and a front.  It does
 * handle equality.
 */
int
add_binned_slug(t_o_domain* domain, slug* (*bin_slug), int bin)
{
  slug* tmp_slug = domain->top_slug[bin];

  if ((*bin_slug)->top == domain->surface_front[bin]
                                                && (domain->yes_groundwater && (*bin_slug)->bot == domain->groundwater_front[bin]))
    {
      //not likely, but this could happen, if it does, we have filled bin
      domain->surface_front[bin]     = domain->layer_top_depth;
      domain->groundwater_front[bin] = domain->layer_top_depth;
      slug_dealloc(&(*bin_slug));
      return 0;
    }
  else if ((*bin_slug)->top == domain->surface_front[bin])
    {
      if (tmp_slug != NULL && tmp_slug->top == (*bin_slug)->bot)
        {
          domain->surface_front[bin] = tmp_slug->bot;
          kill_slug(domain, bin, tmp_slug);
        }
      else
        {
          domain->surface_front[bin] = (*bin_slug)->bot;
        }
      slug_dealloc(&(*bin_slug));
      return 0;
    }
  else if (domain->yes_groundwater && (*bin_slug)->bot == domain->groundwater_front[bin])
    {
      tmp_slug = domain->bot_slug[bin];
      if (tmp_slug != NULL && tmp_slug->bot == (*bin_slug)->top)
        {
          domain->groundwater_front[bin] = tmp_slug->top;
          kill_slug(domain, bin, tmp_slug);
        }
      else
        {
          domain->groundwater_front[bin] = (*bin_slug)->top;
        }
      slug_dealloc(&(*bin_slug));
      return 0;
    }

  //TODO MAY BE POSSIBLE TO REFACTOR THIS PORTION INTO A MORE CONSICE
  //BIT OF CODE
  if (tmp_slug == NULL )
    {
      //this is the first slug to be added to this bin
      //Make sure the slugs pointers are NULL
      (*bin_slug)->next = NULL;
      (*bin_slug)->prev = NULL;
      domain->top_slug[bin] = (*bin_slug);
      domain->bot_slug[bin] = (*bin_slug);
      return 0;
    }
  //Otherwise slugs exist and we need to insert this in
  //the correct position in the linked list
  while (tmp_slug != NULL )
    {
      if ((*bin_slug)->bot < tmp_slug->top)
        {
          //bin_slug goes here
          //rehook pointers
          (*bin_slug)->prev = tmp_slug->prev;
          tmp_slug->prev = (*bin_slug);
          (*bin_slug)->next = tmp_slug;
          if((*bin_slug)->prev == NULL)
            {
              domain->top_slug[bin] = (*bin_slug);
            }
          else
            {
              (*bin_slug)->prev->next = (*bin_slug);
            }
          //Slug is inserted now!
          return 0;
        }
      else if ((*bin_slug)->bot == tmp_slug->top)
        {
          //the two slugs merge into one slug
          tmp_slug->top = (*bin_slug)->top;
          slug_dealloc(&(*bin_slug));
          return 0;
        }
      else if ((*bin_slug)->top == tmp_slug->bot && NULL != tmp_slug->next && (*bin_slug)->bot == tmp_slug->next->top)
        {
          //the three slugs merge into one slug
          tmp_slug->bot = tmp_slug->next->bot;
          slug_dealloc(&(*bin_slug));
          kill_slug(domain, bin, tmp_slug->next);
          return 0;
        }
      else if ((*bin_slug)->top == tmp_slug->bot)
        {
          //the two slugs merge into one slug
          tmp_slug->bot = (*bin_slug)->bot;
          slug_dealloc(&(*bin_slug));
          return 0;
        }
      tmp_slug = tmp_slug->next;
    }
  //if we get to here, we know that the slug is at the end
  //of the binned list
  domain->bot_slug[bin]->next = (*bin_slug);
  (*bin_slug)->prev = domain->bot_slug[bin];
  (*bin_slug)->next = NULL;
  domain->bot_slug[bin] = (*bin_slug);
  return 0;
}

int
redistribute_top_slugs(t_o_domain* domain, slug* (*slugs_head), slug* (*slugs_end),
    int first_bin)
{
  double surface_max = domain->surface_front[first_bin];
  //Loop over all slugs that need to be re-arranged
  int i;
  slug* tmp;
  slug* collide;
  while ((*slugs_head) != NULL )
    {
      //find the first bin the bottom of the slug can contribute to
      for(i = first_bin; i <= domain->parameters->num_bins; i++)
        {
          if((*slugs_head)->bot > domain->surface_front[i])
            {
              break;
            }
        }

      while((*slugs_head) != NULL)
        {
          tmp = domain->bot_slug[i];
          collide = NULL;
          while(tmp != NULL)
            {
              //find the one, if any, slug that will potentially collide
              if(tmp->bot <= surface_max)
                {
                  //this slug is the last in the top of the bin
                  //this will be the only slug of interest;
                  collide = tmp;
                  break;
                }
              tmp = tmp->prev;
            }
          if(collide == NULL)
            {
              if((*slugs_head)->top >= domain->surface_front[i])
                {
                  //slug fits entirely in this bin
                  tmp = (*slugs_head)->next;
                  add_binned_slug(domain, &(*slugs_head), i);
                  (*slugs_head) = tmp;
                  break;
                }
              //otherwise it merges into surface front
              double new_bot = domain->surface_front[i];
              domain->surface_front[i] = (*slugs_head)->bot;
              (*slugs_head)->bot = new_bot;
              i++;
              continue;
            }
          else if((*slugs_head)->bot <= collide->bot)
            {
              i++;
              continue;
            }
          if((*slugs_head)->top >= collide->bot)
            {
              //slug fits entirely in this bin
              tmp = (*slugs_head)->next;
              add_binned_slug(domain, &(*slugs_head), i);
              (*slugs_head) = tmp;
              break;
            }//otherwise it may merge with the bottom slug
          else
            {
              //slug must collide
              //merge slugs and continue
              double new_bot = collide->bot;
              collide->bot = (*slugs_head)->bot;
              (*slugs_head)->bot = new_bot;
              i++;
              continue;
            }
        }
    }
  return 0;
}

int
redistribute_mid_slugs(t_o_domain* domain, slug* (*slugs_head), slug* (*slugs_end),
    int first_bin)
{
  //Loop over all slugs that need to be re-arranged
  int i;
  double ground_max = domain->yes_groundwater ? domain->groundwater_front[first_bin] : domain->layer_bottom_depth;
  slug* tmp;
  slug* collide;
  while ((*slugs_head) != NULL )
    {
      //find the first bin the bottom of the slug can contribute to
      for(i = first_bin; i <= domain->parameters->num_bins; i++)
        {
          if((!domain->yes_groundwater || (*slugs_head)->bot <= domain->groundwater_front[i])
              &&(*slugs_head)->bot > domain->surface_front[i])
            {
              //slug can put water in this bin under the surface front
              //and above the groundwater front, we check for > surface_front because
              //the bottom of the slug should contribute to next bin if there is =
              break;
            }
        }

      while((*slugs_head) != NULL)
        {
          //find the first slug in the middle cut, if one exists
          //this is the slug that can potentially collide
          collide = NULL;
          tmp = domain->bot_slug[i];
          while(tmp != NULL)
            {
              if(tmp->top < ground_max)
                {
                  collide = tmp;
                  break;
                }
              tmp = tmp->prev;
            }

          //now see if slugs_head can contribute
          if(collide == NULL)
            {
              if((*slugs_head)->top >= domain->surface_front[i])
                {
                  tmp = (*slugs_head)->next;
                  add_binned_slug(domain, &(*slugs_head), i);
                  (*slugs_head)= tmp;
                  break;
                }
              else
                {
                  //otherwise it merges into surface front
                  double new_bot = domain->surface_front[i];

                  // surface front might merge with groundwater.
                  if (domain->yes_groundwater && (*slugs_head)->bot == domain->groundwater_front[i])
                    {
                      domain->groundwater_front[i] = domain->layer_top_depth;
                      domain->surface_front[i]     = domain->layer_top_depth;
                    }
                  else if (NULL != domain->top_slug[i] && (*slugs_head)->bot == domain->top_slug[i]->top)
                    {
                      // surface front might merge with the slug below it.
                      domain->surface_front[i] = domain->top_slug[i]->bot;
                      kill_slug(domain, i, domain->top_slug[i]);
                    }
                  else
                    {
                      domain->surface_front[i] = (*slugs_head)->bot;
                    }
                  (*slugs_head)->bot = new_bot;
                  i++;
                  continue;
                }
            }
          if((*slugs_head)->bot <= collide->bot)
            {
              i++;
              continue;
            }
          else if((*slugs_head)->top >= collide->bot)
            {
              //entire slug fits under first_mid_slug
              tmp = (*slugs_head)->next;
              add_binned_slug(domain, &(*slugs_head), i);
              (*slugs_head)= tmp;
              break;
            }
          else
            {
              //merge and continue
              double new_bot = collide->bot;

              // collide might merge with groundwater.
              if (domain->yes_groundwater && (*slugs_head)->bot == domain->groundwater_front[i])
                {
                  domain->groundwater_front[i] = collide->top;
                  kill_slug(domain, i, collide);
                }
              else if (NULL != collide->next && (*slugs_head)->bot == collide->next->top)
                {
                  // collide might merge with the slug below it.
                  collide->bot = collide->next->bot;
                  kill_slug(domain, i, collide->next);
                }
              else
                {
                  collide->bot = (*slugs_head)->bot;
                }

              (*slugs_head)->bot = new_bot;
              i++;
              continue;
            }

        }
    }
  return 0;
}

int
redistribute_bot_slugs(t_o_domain* domain, slug* (*slugs_head), slug* (*slugs_end),
    int first_bin)
{
  double ground_max = domain->groundwater_front[first_bin];
  //Loop over all slugs that need to be re-arranged
  int i;
  slug* tmp;
  slug* collide;

  while ((*slugs_head) != NULL )
    {
      //find the first bin the top of the slug can contribute to
      for(i = first_bin; i <= domain->parameters->num_bins; i++)
        {
          if((*slugs_head)->top < domain->groundwater_front[i])
            {
              break;
            }
        }

      while((*slugs_head) != NULL)
        {
          collide = NULL;
          tmp = domain->top_slug[i];
          while(tmp != NULL)
            {
              //find the one, if any, slug that will potentially collide
              if(tmp->top >= ground_max)
                {
                  //this slug is the first in the bottom of the bin
                  //this will be the only slug of interest;
                  collide = tmp;
                  break;
                }
              tmp = tmp->next;
            }

          if(collide == NULL)
            {
              if((*slugs_head)->bot <= domain->groundwater_front[i])
                {
                  //slug fits entirely in this bin
                  tmp = (*slugs_head)->next;
                  add_binned_slug(domain, &(*slugs_head), i);
                  (*slugs_head) = tmp;
                  break;
                }
              //otherwise it merges into groundwater front
              double new_top = domain->groundwater_front[i];
              domain->groundwater_front[i] = (*slugs_head)->top;
              (*slugs_head)->top = new_top;
              i++;
              continue;
            }
          else if((*slugs_head)->top >= collide->top)
            {
              //water cannot go in this bin above the slug
              i++;
              continue;
            }
          if((*slugs_head)->bot <= collide->top)
            {
              //slug fits entirely in this bin
              tmp = (*slugs_head)->next;
              add_binned_slug(domain, &(*slugs_head), i);
              (*slugs_head) = tmp;
              break;
            }
          else
            {
              //slug must collide
              //merge slugs and continue
              double new_top = collide->top;
              collide->top = (*slugs_head)->top;
              (*slugs_head)->top = new_top;
              i++;
              continue;
            }
        }
    }
  return 0;
}

/* Redistribute water within the domain sideways-tetris-style
 * with no vertical movement of water so that at all depths
 * there is no wet bin to the right of a dry bin.
 * Return TRUE if there is an error, FALSE otherwise.
 * If there is an error some but not all of the redistribution
 * might have been done.
 *
 * Parameters:
 *
 * domain    - A pointer to the t_o_domain struct.
 * first_bin - The leftmost bin that is not completely full of water.
 *             t_o_redistribute can change first_bin, but we are not passing it
 *             by reference and updating it because this is the last step and
 *             we will call find_first_bin anew for the next timestep.
 */
int t_o_redistribute(t_o_domain* domain, int first_bin)
{
  //OPTIMIZATION, ONLY SORT FROM FIRST NON-ZERO BIN ONWARDS
  //ADDITIONALLY, ONLY NEED TO REDISTRIBUTE SLUGS FROM FIRST NON-ZERO BINS
  //TODO IF WE KNOW LAST BIN, WE ONLY HAVE TO SORT BETWEEN FIRST AND LAST BIN.
  int error = FALSE;
  int i;
  int old_first_bin = first_bin;

  //All bins are full, no redistribution necessary
  if (first_bin > domain->parameters->num_bins)
    {
      return 0;
    }
  //First sort the surface_front bins
  //TODO TRY MERGE SORT
  qsort((domain->surface_front) + first_bin,
      domain->parameters->num_bins - first_bin + 1,
      sizeof(*(domain->surface_front)), (void *) compare_surface);
  //Then sort the groundwater_front bins
  if(domain->yes_groundwater)
    {
      qsort((domain->groundwater_front) + first_bin,
          domain->parameters->num_bins - first_bin + 1,
          sizeof(*(domain->groundwater_front)), (void *) compare_ground);
    }
  //Redistribute slugs
  //This requires finding all slugs that exist in the domain
  //and combining them with slugs from the collision of
  //groundwater and surfacewater

  slug* slugs_head = NULL; // The head of a doubly linked list of slugs to insert into the domain.
  slug* slugs_end = NULL; // The tail of a doubly linked list of slugs to insert into the domain.

  slug* slugs_top = NULL;
  slug* slugs_top_end = NULL;
  slug* slugs_bot = NULL;
  slug* slugs_bot_end = NULL;
  slug* slugs_mid = NULL;
  slug* slugs_mid_end = NULL;

  if(domain->yes_groundwater && domain->surface_front[first_bin] != domain->layer_top_depth)
    {
      error = find_collisions(domain, &slugs_top, &slugs_top_end, &slugs_bot,
          &slugs_bot_end, &slugs_mid, &slugs_mid_end, &first_bin);
    }
  if(!error)
    {
      //need to use old first bin here in case there were slugs in a bin that was filled
      //by find_collisions
      for (i = domain->parameters->num_bins; i >= old_first_bin; i--)
        {
          if (domain->top_slug[i] != NULL )
            {
              slug* next = domain->top_slug[i];
              slug* tmp;
              while (next != NULL )
                {
                  tmp = next->next;
                  //TODO POSSIBLE OPTIMIZATION, CALL CUT_SLUGS ON EACH SLUG
                  //AS IT IS FOUND, AND INSERT INSIDE CUT_SLUGS
                  //FIXME DOING THIS OPTIMIZTION RESULTS IN A BUG??? AN INFINITE LOOP SOMEWHERE
                  insert_into_list_top_sort(&slugs_head, &slugs_end, &next);
                  next = tmp;
                }
              domain->top_slug[i] = NULL;
              domain->bot_slug[i] = NULL;
            }
        }

      //put all the slugs in their appropriate "sections"
      error = cut_slugs(domain, &slugs_head, &(slugs_top),
          &(slugs_top_end), &(slugs_bot), &(slugs_bot_end), &(slugs_mid),
          &(slugs_mid_end), first_bin);

      //We can now start to redistribute the slugs
      if(!error)
        {
          //we redistribute in this order so that we only need to check collisions
          //for middle slugs...
          redistribute_top_slugs(domain, &slugs_top, &slugs_top_end, first_bin);
          if(domain->yes_groundwater)
            {
              redistribute_bot_slugs(domain, &slugs_bot, &slugs_bot_end, first_bin);
            }
          else
            {
              assert(slugs_bot == NULL);
            }
          redistribute_mid_slugs(domain, &slugs_mid, &slugs_mid_end, first_bin);
        }
    }
  return error;
}

/* The operation of the code can shave slivers off of slugs creating extra slug
 * structs that take time to process, but are tiny and shouldn't really exist.
 * Move that water down to whatever is below it.
 *
 * Parameters:
 *
 * domain    - A pointer to the t_o_domain struct.
 */
void t_o_handle_sliver_slugs(t_o_domain* domain)
{
  int ii; // Loop counter.

  for (ii = 1; ii < domain->parameters->num_bins; ii++)
    {
      slug* temp_slug = domain->bot_slug[ii];

      while (NULL != temp_slug)
        {
          slug*  prev_slug = temp_slug->prev;
          double slug_size = temp_slug->bot - temp_slug->top;

          if (SLIVER_SLUG_SIZE >= slug_size)
            {
              if (NULL != temp_slug->next)
                {
                  // Put the water in the next lower slug.
                  temp_slug->next->top -= slug_size;
                  kill_slug(domain, ii, temp_slug);
                }
              else if (domain->yes_groundwater)
                {
                  // Put the water in the groundwater.
                  domain->groundwater_front[ii] -= slug_size;
                  kill_slug(domain, ii, temp_slug);
                }
              else
                {
                  // Put the water at the bottom of the domain.
                  temp_slug->top = domain->layer_bottom_depth - slug_size;
                  temp_slug->bot = domain->layer_bottom_depth;
                }
            }

          temp_slug = prev_slug;
        }
    }
}

/* Comment in .h file */
int t_o_timestep(t_o_domain* domain, double dt, double surfacewater_head, double* surfacewater_depth, double water_table, double* groundwater_recharge)
{
  int error = FALSE;    // Error flag.
  int ponded_water = 0; // Flag set by t_o_satisfy_saturated_bins that must be passed to t_o_groundwater.
  int first_bin;        // The leftmost bin that is not completely full of water.

  if (NULL == domain)
    {
      fprintf(stderr, "ERROR: domain must not be NULL\n");
      error = TRUE;
    }

  if (0.0 >= dt)
    {
      fprintf(stderr, "ERROR: dt must be greater than zero\n");
      error = TRUE;
    }

  if (0.0 > water_table)
    {
      fprintf(stderr, "ERROR: water_table must be greater than or equal to zero\n");
      error = TRUE;
    }

  if (NULL == surfacewater_depth)
    {
      fprintf(stderr, "ERROR: surfacewater_depth must not be NULL\n");
      error = TRUE;
    }
  else if (0.0 > *surfacewater_depth)
    {
      fprintf(stderr, "ERROR: surfacewater_depth must be greater than or equal to zero\n");
      error = TRUE;
    }

  if (NULL == groundwater_recharge)
    {
      fprintf(stderr, "ERROR: groundwater_recharge must not be NULL\n");
      error = TRUE;
    }

  // FIXME, wencong 6/2/14, add variable inflow_rate for t_o_groundwater.
  double inflow_rate = 0.0;  // Flow rate through fully saturated bins in unit of meters per second.
  int    no_flow     = FALSE;  // FIXME, add no flow lower boundary, Jan. 09, 2015. 
  if (!error)
    {
      first_bin = find_first_bin(domain, 2);
      if (!no_flow)
        {
           double recharge_old = *groundwater_recharge;
          // Add water_table as passing parameter in t_o_satisfy_saturated_bins() 06/17/14.
          error               = t_o_satisfy_saturated_bins(domain, dt, first_bin, surfacewater_depth, &ponded_water, groundwater_recharge, water_table);
          inflow_rate         = (*groundwater_recharge - recharge_old) / dt; 
        }
    }

  if (!error)
    { // FIXME, wencong, add ponded_water flag, infiltrate when ponded_water is TRUE even surfacewater_depth is zero.
      // error = t_o_infiltrate(domain, dt, &first_bin, surfacewater_head, surfacewater_depth, groundwater_recharge);
         error = t_o_infiltrate(domain, dt, &first_bin, surfacewater_head, surfacewater_depth, groundwater_recharge, ponded_water);
    }

  if (!error)
    {
      error = t_o_falling_slugs(domain, dt, first_bin, groundwater_recharge);
    }

  if (!error)
    {
      // FIXME, wencong 6/2/14, add inflow_rate.
      // error = t_o_groundwater(domain, dt, &first_bin, water_table, ponded_water, groundwater_recharge);
      error = t_o_groundwater(domain, dt, &first_bin, water_table, ponded_water, groundwater_recharge, inflow_rate);
    }

  if (!error)
    {
      t_o_handle_sliver_slugs(domain);
    }
  
  if (!error)
    {
      error = t_o_redistribute(domain, first_bin);
    }

  // FIXME Do we want to call this here?  It is also being called by adhydro_check_invariant.
#if (DEBUG_LEVEL & DEBUG_LEVEL_INTERNAL_ASSERTIONS)
  if (!error)
    {
      t_o_check_invariant(domain); 
    }
#endif // (DEBUG_LEVEL & DEBUG_LEVEL_INTERNAL_ASSERTIONS)

  return error;
}

/* Return a conservative estimate of the depth to fill to in order to add
 * groundwater_recharge to groundwater.  This estimate is achieved by assuming
 * that all of the space above groundwater is empty.  If it really is empty
 * this will produce the exact depth to fill groundwater to.  If there are some
 * slugs or surface front water then add_recharge will not add all of the water
 * and you will have to loop until all of the water is added or the domain is
 * full.
 * 
 * Parameters:
 * 
 * domain               - A pointer to the t_o_domain struct.
 * groundwater_recharge - The amount of water to add to the domain in meters of
 *                        water.
 */
double find_recharge_depth(t_o_domain* domain, double groundwater_recharge)
{
  int    ii              = domain->parameters->num_bins; // Loop counter.
  int    done            = FALSE;                        // Loop flag.
  double space_available = 0.0;                          // Space found so far in meters of water.
  double depth           = 0.0;                          // Meters.  Use zero if we never find enough space for the water.
  
  assert(NULL != domain && 0.0 <= groundwater_recharge);

  while(!done && 1 < ii && domain->layer_top_depth < domain->groundwater_front[ii])
    {
      double new_space= (domain->groundwater_front[ii] - domain->groundwater_front[ii - 1]) * (domain->parameters->num_bins + 1 - ii) *
                         domain->parameters->delta_water_content; // Meters of water.
      
      if (space_available + new_space < groundwater_recharge)
        {
          // Enough water to fill up to domain->groundwater_front[ii - 1].
          space_available += new_space;
        }
      else
        {
          depth = domain->groundwater_front[ii] - (((groundwater_recharge - space_available) / domain->parameters->delta_water_content) / (domain->parameters->num_bins + 1 - ii));
          done = TRUE;
        }
      
      ii--;
    }
  
  assert(domain->layer_top_depth <= depth && depth <= domain->layer_bottom_depth);
  
  return depth;
}

/* Move all groundwater up to at least depth while handling collisions with
 * slugs and surface front water.  Water required to fill the space is taken
 * from groundwater_recharge.
 *
 * Parameters:
 * 
 * domain               - A pointer to the t_o_domain struct.
 * depth                - The depth to fill the groundwater to.
 * groundwater_recharge - A scalar passed by reference containing an amount of
 *                        water in meters of water.  Water used to fill
 *                        groundwater to depth is subtracted from this amount
 *                        of water.  If more water is used than available
 *                        groundwater_recharge will go negative.  It is the
 *                        responsibility of the caller to ensure that filling
 *                        groundwater to depth will not use more water than
 *                        available.
 */
void add_recharge(t_o_domain* domain, double depth, double* groundwater_recharge)
{
  int ii = domain->parameters->num_bins; // Loop counter.
  
  assert(NULL != domain && 0.0 <= depth && NULL != groundwater_recharge && 0.0 <= *groundwater_recharge);
  
  // Loop from right to left until you hit a bin where the groundwater depth is already depth or higher.

  while (0 < ii && depth < domain->groundwater_front[ii])
    {
      // Check for collisions with slugs.
      while (NULL != domain->bot_slug[ii] && domain->bot_slug[ii]->bot >= depth)
        {
          *groundwater_recharge         -= (domain->groundwater_front[ii] - domain->bot_slug[ii]->bot) * domain->parameters->delta_water_content;
          domain->groundwater_front[ii]  = domain->bot_slug[ii]->top;
          kill_slug(domain, ii, domain->bot_slug[ii]);
        }
      
      // Check for collision with surface front water.
      if (domain->surface_front[ii] >= depth)
        {
          *groundwater_recharge         -= (domain->groundwater_front[ii] - domain->surface_front[ii]) * domain->parameters->delta_water_content;
          domain->surface_front[ii]      = domain->layer_top_depth;
          domain->groundwater_front[ii]  = domain->layer_top_depth;
        }
      else if (domain->groundwater_front[ii] > depth) // Must check in case a slug moved groundater.
        {
          *groundwater_recharge         -= (domain->groundwater_front[ii] - depth) * domain->parameters->delta_water_content;
          domain->groundwater_front[ii]  = depth;
        }
      
      ii--;
    }
}

/* Comment in .h file */
void t_o_add_groundwater(t_o_domain* domain, double* groundwater_recharge)
{
  int done = FALSE; // Loop flag.
  
  assert(NULL != domain && NULL != groundwater_recharge && 0.0 <= *groundwater_recharge);
  
  if(domain->yes_groundwater)
    {
      while (!done && epsilon_less(0.0, *groundwater_recharge))
        {
          double depth = find_recharge_depth(domain, *groundwater_recharge);
          
          add_recharge(domain, depth, groundwater_recharge);
          
          if (0.0 == depth)
            {
              done = TRUE;
            }
        }
    }
  else
    {
      fprintf(stderr, "WARNING: called t_o_add_groundwater on a t_o_domain with yes_groundwater FALSE.\n");
    }
}

/* Comment in .h file */
void t_o_take_groundwater(t_o_domain* domain, double water_table, double* groundwater_recharge)
{
  int ii; // Loop counter.

  if(domain->yes_groundwater)
    {
      for (ii = domain->parameters->num_bins; epsilon_greater(0.0, *groundwater_recharge) && ii >= 1; ii--)
        {
          double maximum_bin_depth = water_table - 0.1 * domain->parameters->bin_capillary_suction[ii]; // Meters.

          // If there is water to the right of the groundwater in this bin we can't orphan that water.
          if (ii < domain->parameters->num_bins)
            {
              if (domain->surface_front[ii + 1] > domain->groundwater_front[ii])
                {
                  maximum_bin_depth = domain->layer_top_depth;
                }

              // Find the highest slug, if any, that the groundwater can't drop below.
              slug* temp_slug = domain->bot_slug[ii + 1];

              while (NULL != temp_slug && NULL != temp_slug->prev && temp_slug->prev->bot > domain->groundwater_front[ii])
                {
                  temp_slug = temp_slug->prev;
                }

              if (NULL != temp_slug && temp_slug->bot > domain->groundwater_front[ii] && maximum_bin_depth > temp_slug->top)
                {
                  maximum_bin_depth = temp_slug->top;
                }
              
              if (maximum_bin_depth > domain->groundwater_front[ii + 1])
                {
                  maximum_bin_depth = domain->groundwater_front[ii + 1];
                }
            }
          
          // Don't let groundwater fall below the bottom of the layer.
          if (maximum_bin_depth > domain->layer_bottom_depth)
            {
              maximum_bin_depth = domain->layer_bottom_depth;
            }

          if (maximum_bin_depth > domain->groundwater_front[ii])
            {
              double water_available = (maximum_bin_depth - domain->groundwater_front[ii]) * domain->parameters->delta_water_content; // Meters of water.

              if (water_available <= -*groundwater_recharge)
                {
                  // There is not enough water.  Take it all.
                  domain->groundwater_front[ii]  = maximum_bin_depth;
                  *groundwater_recharge         += water_available;
                }
              else // if (water_available > -*groundwater_recharge)
                {
                  // There is enough water.  Take what you need.
                  domain->groundwater_front[ii] -= *groundwater_recharge / domain->parameters->delta_water_content;
                  *groundwater_recharge          = 0.0;
                }
            }
        }
    }
  else
    {
      fprintf(stderr, "WARNING: called t_o_take_groundwater on a t_o_domain with yes_groundwater FALSE.\n");
    }
}

/* Comment in .h file. */
double t_o_specific_yield(t_o_domain* domain, double water_table)
{
  double porosity               = domain->parameters->bin_water_content[domain->parameters->num_bins];
  double residual_saturation    = (domain->parameters->bin_water_content[1] - domain->parameters->delta_water_content);
  // FIXME store residual saturation in bin_water_content[0] or its own struct member?
  double specific_yield         = (porosity - residual_saturation) *
                                  (1.0 - pow(domain->parameters->bc_psib / (domain->parameters->bc_psib + water_table + 0.01), domain->parameters->bc_lambda));

  if (0.1 > specific_yield)
    {
      specific_yield = 0.1;
    }
  
  if (specific_yield < residual_saturation)
    {
      specific_yield = residual_saturation;
    }
  
  if (specific_yield > porosity)
    {
      specific_yield = porosity;
    }

  return specific_yield;
}

   /******************************************************************************/
  /* The code below is for an old version of t_o_redistribute.  It is only kept */
 /*  around to check the correctness of the new version of t_o_redistribute.   */
/******************************************************************************/

/* Return the next depth at which whether a bin has water will flip.
 * Will return domain->layer_depth if the bin never flips below depth.
 *
 * Parameters:
 *
 * domain - A pointer to the t_o_domain struct.
 * bin    - Which bin to check for flip depth.  One based indexing is used.
 * depth  - Find the next flip strictly below depth.
 */
double next_flip(t_o_domain* domain, int bin, double depth)
{
  double flip;
  // FIXME, wencong
  assert(NULL != domain && 0 < bin && bin <= domain->parameters->num_bins && 0.0 <= depth && depth <= domain->layer_bottom_depth);

  flip = domain->layer_bottom_depth;

  if (depth < domain->surface_front[bin])
    {
      // has_water will flip when it hits the bottom of the surface front water.
      flip = domain->surface_front[bin];
    }
  else
    {
      // Find the first slug deeper than depth, if any.
      slug* temp_slug = domain->top_slug[bin];

      while (NULL != temp_slug && depth >= temp_slug->bot)
        {
          temp_slug = temp_slug->next;
        }

      if (NULL != temp_slug && depth < temp_slug->top)
        {
          // has_water will flip when it hits the top of the slug.
          flip = temp_slug->top;
        }
      else if (NULL != temp_slug && depth < temp_slug->bot)
        {
          // has_water will flip when it hits the bottom of the slug.
          flip = temp_slug->bot;
        }
      else if (domain->yes_groundwater && depth < domain->groundwater_front[bin])
        {
          // has_water will flip when it hits the top of the groundwater.
          flip = domain->groundwater_front[bin];
        }
    }

  return flip;
}

/* Remove water between top and bot from bin.
 * This function assumes that bin has water from top to bot.
 * If yes_groundwater is FALSE, you cannot use this function to
 * remove water from bins less than or equal to initial_water_content.
 * Return TRUE if there is an error, FALSE otherwise.
 * If there is an error no water is removed.
 *
 * Parameters:
 *
 * domain - A pointer to the t_o_domain struct.
 * bin    - Which bin to remove water from.  One based indexing is used.
 * top    - The depth of the top of the region of water to remove.
 * bot    - The depth of the bottom of the region of water to remove.
 */
int remove_water(t_o_domain* domain, int bin, double top, double bot)
{
  int error = FALSE; // Error flag.
  //FIXME, wencong, layer_bottom_depth.
  assert(NULL != domain && 0 < bin && bin <= domain->parameters->num_bins && domain->layer_top_depth <= top && top < bot && bot <= domain->layer_bottom_depth);

  if (top < domain->surface_front[bin])
    {
      // The water is in the surface front water.
      assert(bot <= domain->surface_front[bin]);

      if (bot < domain->surface_front[bin])
        {
          // Need to create a slug from bot to surface_front.
          error = create_slug_after(domain, bin, NULL, bot, domain->surface_front[bin]);
        }

      if (!error)
        {
          // Surface front water now goes down to top.
          domain->surface_front[bin] = top;
        }
    }
  else
    {
      // The water is not in the surface front water.
      // Try to find it in a slug.
      slug* temp_slug = domain->top_slug[bin];

      while (NULL != temp_slug && top >= temp_slug->bot)
        {
          temp_slug = temp_slug->next;
        }

      if (NULL != temp_slug)
        {
          // The water is in temp_slug.
          assert(top >= temp_slug->top && bot <= temp_slug->bot);

          if (top == temp_slug->top && bot == temp_slug->bot)
            {
              // Get the whole slug.
              kill_slug(domain, bin, temp_slug);
            }
          else if (top == temp_slug->top)
            {
              // Get water from the top of the slug.
              temp_slug->top = bot;
            }
          else if (bot == temp_slug->bot)
            {
              // Get water from the bottom of the slug.
              temp_slug->bot = top;
            }
          else
            {
              // Get water from the middle of the slug.
              // need to create a slug from bot to temp_slug->bot.
              error = create_slug_after(domain, bin, temp_slug, bot, temp_slug->bot);

              if (!error)
                {
                  // The old slug now goes down to top.
                  temp_slug->bot = top;
                }
            }
        } // End the water is in temp_slug.
      else
        {
          // The water is in the groundwater.
          if (!domain->yes_groundwater)
            {
              fprintf(stderr, "ERROR: Function remove_water attempting to remove water from a bin less than or equal to initial_water_content with "
                  "yes_groundwater equal to FALSE.\n");
              error = TRUE;
            }

          if (!error)
            {
              assert(top >= domain->groundwater_front[bin]);

              if (top > domain->groundwater_front[bin])
                {
                  if (0.0 == domain->groundwater_front[bin])
                    {
                      // The water above the removed water is surface front water
                      domain->surface_front[bin] = top;
                    }
                  else
                    {
                      // need to create a slug from groundwater_front to top
                      error = create_slug_after(domain, bin, domain->bot_slug[bin], domain->groundwater_front[bin], top);
                    }
                }
            }

          if (!error)
            {
              // Groundwater now starts at bot.
              domain->groundwater_front[bin] = bot;
            }
        } // End the water is in the groundwater.
    } // End the water is not in the surface attched water.

  return error;
}

/* Add water between top and bot to bin.
 * This function assumes that the added water does not touch surface
 * front water or groundwater and that bin is empty from top to bot.
 * Return TRUE if there is an error, FALSE otherwise.
 * If there is an error no water is added.
 *
 * Parameters:
 *
 * domain - A pointer to the t_o_domain struct.
 * bin    - Which bin to put water in.  One based indexing is used.
 * top    - The depth of the top of the region of water to add.
 * bot    - The depth of the bottom of the region of water to add.
 */
int add_water_to_slug(t_o_domain* domain, int bin, double top, double bot)
{
  int error = FALSE; // Error flag.

  assert(NULL != domain && 0 < bin && bin <= domain->parameters->num_bins && domain->surface_front[bin] < top && top < bot);

  if (domain->yes_groundwater)
    {
      assert(bot < domain->groundwater_front[bin]);
    }
  else
    {
      // FIXME, wencong.
      assert(bot <= domain->layer_bottom_depth);
    }

  if (NULL == domain->top_slug[bin] || bot < domain->top_slug[bin]->top)
    {
      // Add the first slug in the bin or a slug above the top slug.
      error = create_slug_after(domain, bin, NULL, top, bot);
    }
  else if (bot == domain->top_slug[bin]->top)
    {
      // The water is attached to the top of the top slug.
      domain->top_slug[bin]->top = top;
    }
  else if (top > domain->bot_slug[bin]->bot)
    {
      // Add a slug below the bottom slug.
      error = create_slug_after(domain, bin, domain->bot_slug[bin], top, bot);
    }
  else if (top == domain->bot_slug[bin]->bot)
    {
      // The water is attached to the bottom of the bottom slug.
      domain->bot_slug[bin]->bot = bot;
    }
  else
    {
      // Loop over the gaps to find the gap that the slug is in.
      // Look at the gap between temp_slug and temp_slug->next.
      slug* temp_slug = domain->top_slug[bin];

      while (temp_slug->next->next != NULL && top >= temp_slug->next->bot)
        {
          temp_slug = temp_slug->next;
        }

      assert(top >= temp_slug->bot && bot <= temp_slug->next->top);

      if (top == temp_slug->bot && bot == temp_slug->next->top)
        {
          // Merge the two slugs.
          temp_slug->bot = temp_slug->next->bot;
          kill_slug(domain, bin, temp_slug->next);
        }
      else if (top == temp_slug->bot)
        {
          // Add the water to the bottom of temp_slug.
          temp_slug->bot = bot;
        }
      else if (bot == temp_slug->next->top)
        {
          // Add the water to the top of temp_slug->next.
          temp_slug->next->top = top;
        }
      else
        {
          // Add a new slug in between.
          error = create_slug_after(domain, bin, temp_slug, top, bot);
        }
    }

  return error;
}

/* Add water between top and bot to bin.
 * This function assumes that bin is empty from top to bot.
 * Return TRUE if there is an error, FALSE otherwise.
 * If there is an error no water is added.
 *
 * Parameters:
 *
 * domain - A pointer to the t_o_domain struct.
 * bin    - Which bin to put water in.  One based indexing is used.
 * top    - The depth of the top of the region of water to add.
 * bot    - The depth of the bottom of the region of water to add.
 */
int add_water(t_o_domain* domain, int bin, double top, double bot)
{
  int error = FALSE; // Error flag.
  //FIXME, wencong.
  assert(NULL != domain && 0 < bin && bin <= domain->parameters->num_bins && domain->surface_front[bin] <= top && top < bot && bot <= domain->layer_bottom_depth);

  if (top == domain->surface_front[bin])
    {
      // Add the water to the bottom of the surface front water.
      if (NULL != domain->top_slug[bin])
        {
          assert(bot <= domain->top_slug[bin]->top);

          if (bot == domain->top_slug[bin]->top)
            {
              // The surface front water has joined with the top slug.
              domain->surface_front[bin] = domain->top_slug[bin]->bot;
              kill_slug(domain, bin, domain->top_slug[bin]);
            }
          else
            {
              domain->surface_front[bin] = bot;
            }
        }
      else if (domain->yes_groundwater)
        {
          assert(bot <= domain->groundwater_front[bin]);

          if (bot == domain->groundwater_front[bin])
            {
              // The surface front water has joined with the groundwater.
              // FIXME, wencong, change two 0.0.
              domain->surface_front[bin]     = domain->layer_top_depth;
              domain->groundwater_front[bin] = domain->layer_top_depth;
            }
          else
            {
              domain->surface_front[bin] = bot;
            }
        }
      else
        {
          domain->surface_front[bin] = bot;
        }
    } // End add the water to the bottom of the surface front water.
  else
    {
      // Add the water to the groundwater or a slug.
      if (domain->yes_groundwater)
        {
          assert(bot <= domain->groundwater_front[bin]);

          if (bot == domain->groundwater_front[bin])
            {
              // Add the water to the top of the groundwater.
              if (NULL != domain->bot_slug[bin])
                {
                  assert(top >= domain->bot_slug[bin]->bot);

                  if (top == domain->bot_slug[bin]->bot)
                    {
                      // The groundwater has joined with the bottom slug.
                      domain->groundwater_front[bin] = domain->bot_slug[bin]->top;
                      kill_slug(domain, bin, domain->bot_slug[bin]);
                    }
                  else
                    {
                      domain->groundwater_front[bin] = top;
                    }
                }
              else
                {
                  domain->groundwater_front[bin] = top;
                }
            }
          else // The bottom of the water we are adding is not touching groundwater.
            {
              // Add the water to a slug.
              error = add_water_to_slug(domain, bin, top, bot);
            }
        }
      else // domain->yes_groundwater is FALSE
        {
          // Add the water to a slug.
          error = add_water_to_slug(domain, bin, top, bot);
        }
    } // End add the water to the groundwater or a slug.

  return error;
}

/* Move the water between top and bot from from_bin to to_bin.
 * This function assumes that from_bin has water from top to bot
 * and that to_bin is empty from top to bot.
 * Return TRUE if there is an error, FALSE otherwise.
 * If there is an error it might be that no water is moved, but
 * it is also possible that water is removed from from_bin and errors
 * prevent it from being put either into to_bin or back into from_bin.
 *
 * Parameters:
 *
 * domain   - A pointer to the t_o_domain struct.
 * from_bin - Which bin to take water from.  One based indexing is used.
 * to_bin   - Which bin to put water in.  One based indexing is used.
 * top      - The depth of the top of the region of water to move.
 * bot      - The depth of the bottom of the region of water to move.
 */
int move_water(t_o_domain* domain, int from_bin, int to_bin, double top, double bot)
{
  int error = FALSE; // Error flag.
 
  // FIXME, wencong.
  assert(NULL != domain && 0 < from_bin && from_bin <= domain->parameters->num_bins && 0 < to_bin && to_bin <= domain->parameters->num_bins &&
      from_bin != to_bin && domain->layer_top_depth <= top && top < bot && bot <= domain->layer_bottom_depth && 
      has_water_at_depth(domain, from_bin, top, bot) && !has_water_at_depth(domain, to_bin, top, bot));

  error = remove_water(domain, from_bin, top, bot);

  if (!error)
    {
      error = add_water(domain, to_bin, top, bot);

      if (error)
        {
          // Put the water back.
          if (add_water(domain, from_bin, top, bot))
            {
              // Can't put the water anywhere.
              fprintf(stderr, "ERROR: Lost water in move_water.\n");
            }
        }
    }

  return error;
}

/* This is an old, slower way of doing redistribution.  It is kept around
 * to check the correctness of the new redistribution algorithm.
 * It should be a drop-in replacement for redistribute().
 * See the comment of redistribute() for more details.
 *
 * Algorithm:
 *
 * top       - The top of the depth range currently being processed.
 * bot       - The bottom of the depth range currently being processed.
 *             top and bot will start at the surface and walk down the
 *             entire depth of the domain.
 * has_water - An array of flags, one per bin, representing whether that bin
 *             has water between top and bot.
 * flip      - An array of doubles, one per bin, representing the next depth
 *             at which the bin will flip the value of its has_water flag.
 *
 * Initialize top and bot to zero.
 * Initialize has_water to whether each bin has water at depth zero.
 * Initialize flip to the first depth where a bin flips whether it has water.
 *   This depth may be domain->full_depth if the bin never flips.
 *
 * Do
 *   Set top = bot.
 *   Set bot to the smallest value in flip.
 *   All bins are either entirely wet or entirely dry between top and bot.
 *   Do
 *     Move the water between top and bot from the rightmost bin
 *       that is wet to the leftmost bin that is dry.
 *   Until there is no wet bin to the right of a dry bin.
 *   For all bins where flip == bot:
 *     Flip the has_water flag.
 *     Find the next lower flip, may be domain->full_depth.
 * until bot == domain->full_depth.
 *
 * Parameters:
 *
 * domain - A pointer to the t_o_domain struct.
 */
int redistribute_slow(t_o_domain* domain)
{
  int    error = FALSE; // Error flag.
  int    ii;            // Loop counter.
  double top = 0.0;     // The top of the slice being processed.
  double bot = 0.0;     // The bottom of the slice being processed.

  if (NULL == domain)
    {
      fprintf(stderr, "ERROR: domain must not be NULL\n");
      error = TRUE;
    }

  // has_water[ii] is TRUE if bin ii has water between top and bot.
  int has_water[domain->parameters->num_bins + 1];

  // flip[ii] is the depth at which has_water[ii] will next change or
  // domain->layer_depth if has_water[ii] will never change again.
  double flip[domain->parameters->num_bins + 1];

  // Fill in has_water and flip.
  if (!error)
    {
      for (ii = 1; ii <= domain->parameters->num_bins; ii++)
        {
          // FIXME, wencong, change the following 0.0??
          has_water[ii] = has_water_at_depth(domain, ii, 0.0, 0.0);
          flip[ii]      = next_flip(domain, ii, 0.0);
        }
    }

  // Walk down the entire depth of the domain.
  // FIXME, wencong.
  while (!error && bot < domain->layer_bottom_depth)
    {
      // Set top to bot and bot to the closest flip.
      top = bot;
      bot = flip[1];

      for (ii = 2; ii <= domain->parameters->num_bins; ii++)
        {
          if (bot > flip[ii])
            {
              bot = flip[ii];
            }
        }

      assert(bot > top);

      // At this point all bins are either entirely wet
      // or entirely dry between top and bot.

      // Find the right most wet bin and the leftmost dry bin
      // and move the water from the wet bin to the dry bin.
      // Repeat until there are no wet bins to the right of a dry bin.
      int right_bin = domain->parameters->num_bins;
      int left_bin  = 1;

      while (!error && right_bin > left_bin)
        {
          if (!has_water[right_bin])
            {
              right_bin--;
            }
          else if (has_water[left_bin])
            {
              left_bin++;
            }
          else
            {
              error = move_water(domain, right_bin, left_bin, top, bot);
              right_bin--;
              left_bin++;
            }
        }

      // Flip has_water flag for bins that flip at bot.
      if (!error)
        {
          for (ii = 1; ii <= domain->parameters->num_bins; ii++)
            {
              if (bot == flip[ii])
                {
                  has_water[ii] = !has_water[ii];
                  flip[ii] = next_flip(domain, ii, bot);
                }
            }
        }
    }

  return error;
}

   /*********************************************************************************/
  /* The code below is for an old version of t_o_add_groundwater.  It is only kept */
 /*  around to check the correctness of the new version of t_o_add_groundwater.   */
/*********************************************************************************/

void t_o_add_groundwater_slow(t_o_domain* domain, double* groundwater_recharge)
{
  int ii;                  // Loop counter.
  int domain_full = FALSE; // Whether the domain is completely full of water.

  if (domain->yes_groundwater)
    {
      // FIXLATER This is brute force and could probably be optimized.
      while (!domain_full && 0.0 < *groundwater_recharge)
        {
          double top = 0.0; // The top of the layer to fill with water.
          double bot = 0.0; // The bottom of the layer to fill with water
          int    number_of_bins_at_bot = 0;

          // Find space in the bin(s) with the lowest groundwater.
          for (ii = domain->parameters->num_bins; 1 < ii; ii--)
            {
              if (bot < domain->groundwater_front[ii])
                {
                  bot = domain->groundwater_front[ii];
                  number_of_bins_at_bot = 1;
                }
              else if (bot == domain->groundwater_front[ii])
                {
                  number_of_bins_at_bot++;
                }
              else if (top < domain->groundwater_front[ii])
                {
                  top = domain->groundwater_front[ii];
                }

              if (NULL != domain->bot_slug[ii]
                                           && top < domain->bot_slug[ii]->bot)
                {
                  top = domain->bot_slug[ii]->bot;
                }

              if (top < domain->surface_front[ii])
                {
                  top = domain->surface_front[ii];
                }
            }

          // FIXME, wencong
          assert(domain->layer_top_depth <= top && top <= bot && bot <= domain->layer_bottom_depth);

          double space_available_in_bins = (bot - top) * number_of_bins_at_bot; // Meters of bin depth.

          double depth_to_move_to; // The depth that all bins at bot should move up to.

          // Figure out if that space can take all, some, or none of the available water.
          if (space_available_in_bins * domain->parameters->delta_water_content >= *groundwater_recharge)
            {
              depth_to_move_to = bot - (*groundwater_recharge / (number_of_bins_at_bot * domain->parameters->delta_water_content));
              *groundwater_recharge = 0.0;
            }
          else if (0.0 < space_available_in_bins)
            {
              depth_to_move_to = top;
              *groundwater_recharge -= space_available_in_bins * domain->parameters->delta_water_content;
            }
          else
            {
              domain_full = TRUE; // No more space to put water.
              depth_to_move_to = bot;
            }

          // Move the groundwater front up for the added water.
          for (ii = domain->parameters->num_bins; 1 < ii; ii--)
            {
              if (domain->groundwater_front[ii] == bot)
                {
                  if (domain->surface_front[ii] == depth_to_move_to)
                    {
                      // Groundwater hits the surface front.
                      // FIXME, wencong, change two 0.0.
                      domain->groundwater_front[ii] = domain->layer_top_depth;
                      domain->surface_front[ii]     = domain->layer_top_depth;
                    }
                  else if (NULL != domain->bot_slug[ii] && domain->bot_slug[ii]->bot == depth_to_move_to)
                    {
                      // Groundwater hits a slug.
                      domain->groundwater_front[ii] = domain->bot_slug[ii]->top;
                      kill_slug(domain, ii, domain->bot_slug[ii]);
                    }
                  else
                    {
                      domain->groundwater_front[ii] = depth_to_move_to;
                    }
                }
            }
        } // End while (!domain_full && 0.0 < *groundwater_recharge).
    } // End if (domain->yes_groundwater).
  else
    {
      fprintf(stderr, "WARNING: called t_o_add_groundwater on a t_o_domain with yes_groundwater FALSE.\n");
    }
}

   /*******************************************************************************/
  /* The code below is for debugging only.  It compares two Talbot-Ogden domains */
 /*  to see if they are equal to test that two implementations are equivalent.  */
/*******************************************************************************/

int t_o_domains_equal(t_o_domain* domain1, t_o_domain* domain2)
{
  int equal = TRUE; // Whether the two domains are equal.
  int ii;           // Loop counter.

  equal = NULL != domain1 && NULL != domain2;

  // Test if parameters are equal
  if (equal)
    {
      equal = domain1->parameters->num_bins                    == domain2->parameters->num_bins                    &&
          domain1->parameters->vg_alpha                    == domain2->parameters->vg_alpha                    &&
          domain1->parameters->vg_n                        == domain2->parameters->vg_n                        &&
          domain1->parameters->bc_lambda                   == domain2->parameters->bc_lambda                   &&
          domain1->parameters->bc_psib                     == domain2->parameters->bc_psib                     &&
          domain1->parameters->delta_water_content         == domain2->parameters->delta_water_content         &&
          domain1->parameters->effective_capillary_suction == domain2->parameters->effective_capillary_suction &&
          domain1->parameters->dry_depth_dt                == domain2->parameters->dry_depth_dt;
    }

  for (ii = 1; equal && ii <= domain1->parameters->num_bins; ii++)
    {
      equal = domain1->parameters->bin_water_content[ii]       == domain2->parameters->bin_water_content[ii]       &&
          domain1->parameters->cumulative_conductivity[ii] == domain2->parameters->cumulative_conductivity[ii] &&
          domain1->parameters->bin_capillary_suction[ii]   == domain2->parameters->bin_capillary_suction[ii]   &&
          domain1->parameters->bin_dry_depth[ii]           == domain2->parameters->bin_dry_depth[ii];
    }

  // Test if scalars are equal
  if (equal)
    {
      // FIXME,wencong, domain->layer_bottom_depth.
      equal = domain1->layer_bottom_depth            == domain2->layer_bottom_depth     &&
          domain1->yes_groundwater        == domain2->yes_groundwater &&
          (!domain1->yes_groundwater ||
              domain1->initial_water_content == domain2->initial_water_content);
    }

  // Test if surface fronts are equal
  for (ii = 1; equal && ii <= domain1->parameters->num_bins; ii++)
    {
      equal = domain1->surface_front[ii] == domain2->surface_front[ii];
    }

  // Test if slugs are equal
  for (ii = 1; equal && ii <= domain1->parameters->num_bins; ii++)
    {
      slug* slug1 = domain1->top_slug[ii];
      slug* slug2 = domain2->top_slug[ii];

      while (equal && NULL != slug1)
        {
          equal = NULL != slug2;

          if (equal)
            {
              equal = slug1->top == slug2->top && slug1->bot == slug2->bot;
              slug1 = slug1->next;
              slug2 = slug2->next;
            }
        }

      if (equal)
        {
          equal = NULL == slug2;
        }
    }

  // Test if groundwater fronts are equal
  for (ii = 1; domain1->yes_groundwater && equal && ii <= domain1->parameters->num_bins; ii++)
    {
      equal = domain1->groundwater_front[ii] == domain2->groundwater_front[ii];
    }

  return equal;
}

// FIXME, wencong, add ET, Dec. 10, 2014. A very simple ET function, need to separate evaporation and transpiration.
int t_o_ET(t_o_domain* domain, double dt, double root_depth, double PET, double field_capacity, double wilting_point, 
           int use_feddes, double field_capacity_suction, double wilting_point_suction, double* surfacewater_depth, double* evaporated_water)
{
  int error     = FALSE;
  int bare_soil = FALSE;
  int ii;
  
  if (root_depth < domain->layer_top_depth || PET <= 0.0)
    {
      return error;
    }
 
  if (FALSE == use_feddes)
    {
        assert(domain->parameters->bin_water_content[1] <= wilting_point && wilting_point < field_capacity && 
               field_capacity <= domain->parameters->bin_water_content[domain->parameters->num_bins]);
    }
  
  if (root_depth > domain->layer_bottom_depth)
    {
      root_depth = domain->layer_bottom_depth;
    }
  else if (root_depth < 0.1 + domain->layer_top_depth)
    {
      if (epsilon_equal(0.0, root_depth))
        {
          bare_soil = TRUE;
        }
      root_depth    = 0.1 + domain->layer_top_depth; // If root depth less than 0.1 m, set it as 0.1 m.
    }
  
  // Step 1, calculate actual ET form PET, based on water content of last bin, or water content of last slug with root depth.
  int first_bin        = find_first_bin(domain, 2);
  int last_bin         = find_last_bin(domain);
  double water_content = domain->parameters->bin_water_content[last_bin];
  double suction       = domain->parameters->bin_capillary_suction[last_bin];
  
  for (ii = domain->parameters->num_bins; ii > last_bin; ii--)
     {
       if (NULL != domain->top_slug[ii] && domain->top_slug[ii]->top < root_depth)
         {
           water_content = domain->parameters->bin_water_content[ii];
           suction       = domain->parameters->bin_capillary_suction[ii];
           break;
         }
     }  
  
  double actual_ET = PET;
  //use_feddes = TRUE;
  if (use_feddes && !bare_soil)
     {  // calculate actual ET using Feddes(1978) based on suction as in Hydrus-1D.
       //wilting_point_suction  = 150; //double suction_1 = 150.0 ; // 150 m, or 3 m.##############################
       //field_capacity_suction = 0.3; //double suction_2 = 0.3;
       if (suction > wilting_point_suction)
         {
           actual_ET = 0.0;
         }
       else if (suction > field_capacity_suction)
         {
           actual_ET = PET * (1.0 - (suction - field_capacity_suction) / (wilting_point_suction - field_capacity_suction));
         }
       else
         {
           // Do nothing, actual_ET = PET.
         }
     }
  else if (!bare_soil)
    { // calculate actual ET based on water content, instead of pressure as Feddes.
      if (water_content <= wilting_point)
        {
          actual_ET = 0.0;
        }
      else if (water_content < field_capacity)
        {
          actual_ET = PET * (water_content - wilting_point) / (field_capacity - wilting_point);
        }
      else if (water_content >= field_capacity)
        {
          // Do nothing, actual_ET = PET.
        }
    }
  else if (bare_soil)
    {
      // Bare soil, do nothing, actual_ET = PET.
    }

  // Step 2, remove water.
  double demand_ET    = actual_ET * dt;      // Demand ET water in meter of water.
  if (bare_soil)
    { // Take water from surfacewater_depth.
      if (demand_ET >= *surfacewater_depth)
        {
          demand_ET          -= *surfacewater_depth;
          *evaporated_water  += *surfacewater_depth;
          *surfacewater_depth = 0.0;
        }
      else
        {
          *surfacewater_depth -= demand_ET;
          *evaporated_water   += demand_ET;
          demand_ET            = 0.0;
        }
    }
  double demand_ET_dz = demand_ET / domain->parameters->delta_water_content;    // Demand ET water in meter of bin width water.
  ii = domain->parameters->num_bins;

  while (demand_ET_dz > 0.0 && ii > 1)
    { // Loop to satisfy ET demand.
      // Step 2.1, ET from surface front.
      if (domain->surface_front[ii] > domain->layer_top_depth)
        { 
          double bin_demand_ET_dz = demand_ET_dz;      // Water distnace in a bin to be remove is min(demand_ET_dz, root_depth).
          if (bin_demand_ET_dz > root_depth - domain->layer_top_depth)
            {
              bin_demand_ET_dz = root_depth - domain->layer_top_depth;
            }
          
          if (domain->surface_front[ii] - domain->layer_top_depth <= bin_demand_ET_dz)
            { // Water in a bin is less than demand, remove all.
              *evaporated_water        += (domain->surface_front[ii] - domain->layer_top_depth) * domain->parameters->delta_water_content;
              domain->surface_front[ii] = domain->layer_top_depth;
              demand_ET_dz             -= (domain->surface_front[ii] - domain->layer_top_depth);
            }
          else
            { // Remove water from top, and make it a slug.
              error = detach_slug(domain, ii);
              if (!error)
                {
                  if (domain->top_slug[ii]->top + bin_demand_ET_dz < domain->top_slug[ii]->bot)
                    {
                      *evaporated_water         += bin_demand_ET_dz * domain->parameters->delta_water_content;
                      domain->top_slug[ii]->top += bin_demand_ET_dz;
                      demand_ET_dz              -= bin_demand_ET_dz;
                    }
                }
            }
        }// End of surface front.
        
      // Step 2.2, ET from slugs.
      slug* temp_slug = domain->top_slug[ii]; 
      while (NULL != temp_slug && temp_slug->top < root_depth)
        {
          // Save a pointer to the slug below get_slug in case we need to kill get_slug.
          slug* next_slug = temp_slug->next;
                  
          double bin_demand_ET_dz = demand_ET_dz;      // Water distnace in a bin to be remove is min(demand_ET_dz, root_depth).
          if (bin_demand_ET_dz > root_depth - temp_slug->top)
            {
              bin_demand_ET_dz = root_depth - temp_slug->top;
            }
           
          if (bin_demand_ET_dz >= temp_slug->bot - temp_slug->top)
            {
              demand_ET_dz      -= (temp_slug->bot - temp_slug->top);
              *evaporated_water += (temp_slug->bot - temp_slug->top) * domain->parameters->delta_water_content;
              kill_slug(domain, ii, temp_slug);
            }
          else
            {
              *evaporated_water += bin_demand_ET_dz * domain->parameters->delta_water_content;
              temp_slug->top    += bin_demand_ET_dz;
              demand_ET_dz      -= bin_demand_ET_dz;
              break;
            }
          temp_slug = next_slug;
        } // End of slug.
      
      // Step 2.3, ET from groundwater front. 
      if (domain->yes_groundwater && demand_ET_dz > 0.0)
        {
          double bin_demand_ET_dz = demand_ET_dz;
          if (root_depth > domain->groundwater_front[ii] )
            {
              if (bin_demand_ET_dz > root_depth - domain->groundwater_front[ii])
                {
                  bin_demand_ET_dz = root_depth - domain->groundwater_front[ii];
                }
              
              // Modified Feb, 09, 2015. Originaly outside the loop and it was wrong.
              demand_ET_dz                  -= bin_demand_ET_dz;
              *evaporated_water             += bin_demand_ET_dz * domain->parameters->delta_water_content;
              domain->groundwater_front[ii] += bin_demand_ET_dz;
            }
        }
       
      ii--;
    } // End of while loop.
  
  
  // Step 3, call redistribution.
  if (!error)
    {
      error = t_o_redistribute(domain, first_bin);
    }
  
  return error;
}
