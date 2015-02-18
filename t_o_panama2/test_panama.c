#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <X11/Xlib.h>
#include "t_o.h"
#include "epsilon.h"
#include "all.h"
#include "memfunc.h"

extern int t_o_domains_equal(t_o_domain* domain1, t_o_domain* domain2);

#define MARGIN        (10)
#define FRAME_WIDTH   (1500)
#define FRAME_HEIGHT  (1000)
#define WINDOW_WIDTH  (FRAME_WIDTH  + 2 * MARGIN)
#define WINDOW_HEIGHT (FRAME_HEIGHT + 2 * MARGIN)
#define RED           (0xFF0000)
#define BLUE          (0x0000FF)
#define PURPLE        (RED | BLUE)
#define ONE_MINUTE    (60.0)
#define ONE_HOUR      (60.0 * ONE_MINUTE)
#define ONE_DAY       (24.0 * ONE_HOUR)
//#define YES_PLOT

void display_water(Display* the_display, Pixmap the_pixmap, GC the_gc, t_o_domain* domain, int bin, double top, double bot)
{
  int width_factor  = FRAME_WIDTH  / domain->parameters->num_bins; // Pixels per bin.
  int height_factor = FRAME_HEIGHT / domain->layer_bottom_depth;   // Pixels per meter.
  int x             = MARGIN + (bin - 1) * width_factor;           // Pixels.
  int width         = width_factor;                                // Pixels.
  int y_top         = MARGIN + top * height_factor;                // Pixels.
  int y_bot         = MARGIN + bot * height_factor;                // Pixels
  int height        = y_bot - y_top + 1;                           // Pixels.  Add one to height to force the rectangles to be at least 1 pixel tall.

  XFillRectangle(the_display, the_pixmap, the_gc, x, y_top, width, height);
}

void display_domain(Display* the_display, int the_screen, Window the_window, Pixmap the_pixmap, t_o_domain* domain)
{
  int ii; // Loop counter.
  GC the_gc = DefaultGC(the_display, the_screen);

  // Clear the pixmap.
  XSetForeground(the_display, the_gc, WhitePixel(the_display, the_screen));
  XFillRectangle(the_display, the_pixmap, the_gc, 0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);

  // Fill rectangles where water is.
  for (ii = 1; ii <= domain->parameters->num_bins; ii++)
    {
      if (!domain->yes_groundwater &&  domain->parameters->bin_water_content[ii] <= domain->initial_water_content)
        {
          // Display the completely saturated bin.
          XSetForeground(the_display, the_gc, BLUE);
          display_water(the_display, the_pixmap, the_gc, domain, ii, domain->layer_top_depth, domain->layer_bottom_depth);
        }
      else
        {
          // Display the surface attached water.
          if (0.0 < domain->surface_front[ii])
            {
              XSetForeground(the_display, the_gc, RED);
              display_water(the_display, the_pixmap, the_gc, domain, ii, domain->layer_top_depth, domain->surface_front[ii]);
            }

          // Display each slug.
          slug* temp_slug = domain->top_slug[ii];

          while (NULL != temp_slug)
            {
              XSetForeground(the_display, the_gc, PURPLE);
              display_water(the_display, the_pixmap, the_gc, domain, ii, temp_slug->top, temp_slug->bot);
              temp_slug = temp_slug->next;
            }

          // Display the groundwater.
          if (domain->yes_groundwater)
            {
              XSetForeground(the_display, the_gc, BLUE);
              display_water(the_display, the_pixmap, the_gc, domain, ii, domain->groundwater_front[ii], domain->layer_bottom_depth);
            }
        }
    }

  // Copy pixmap to window.
  XCopyArea(the_display, the_pixmap, the_window, the_gc, 0, 0, WINDOW_WIDTH, WINDOW_HEIGHT, 0, 0);
  XFlush(the_display);
}

/* A link list struct stores depth, theta with depth increasing order. */
struct t_o_profile_link
{
  double depth;
  double theta;
  double pressure;
  struct t_o_profile_link *next;
};
typedef struct t_o_profile_link t_o_profile_link;

/* Insert a node with (depth, theta) into link list (head_pt), return error = TRUE if error occurs. */
int insert_t_o_profile_link(t_o_profile_link** head_pt, double depth, double theta, double pressure, int yes_groundater)
{
  int error = FALSE;
  assert(NULL != *head_pt);
  
  t_o_profile_link* pt_insert;
  t_o_profile_link* pt_rear;
  t_o_profile_link* pt_front;
  
  if (NULL == (pt_insert = (t_o_profile_link*)malloc(sizeof(t_o_profile_link))))
    {
      printf("malloc failed in creating t_o_profile_link node in function insert_t_o_profile_link(), exit.\n");
      error = TRUE;
      return error;
    }
  pt_insert->depth    = depth;
  pt_insert->theta    = theta;
  pt_insert->pressure = pressure;
  pt_insert->next     = NULL;
  
  pt_rear  = *head_pt;
  pt_front = (*head_pt)->next;
  // Head pointer *head_pt is initialized with domain->top_depth, and pass in value (depth) is larger than domain->top_depth.
  while (NULL != pt_front)
    {
      if (pt_insert->depth >= pt_front->depth)
        {
          pt_rear  = pt_front;
          pt_front = pt_front->next;     
        }
      else
        {
          break;
        }
    }
  
  // Insert pt_insert only when depth is strictly larger than pt_rear->depth, so that depth in link is increasing without same value.
  if (epsilon_greater(pt_insert->depth, pt_rear->depth))
    {
      pt_rear->next   = pt_insert;
      pt_insert->next = pt_front;
    }
  else if (!yes_groundater && epsilon_equal(pt_insert->depth, pt_rear->depth) && pt_rear->theta < pt_insert->theta)
    { // Add 10/09/14, to handle no groundwate cases,
     //printf("theta = %lf \n", pt_insert->theta); getchar();
      pt_rear->theta    = pt_insert->theta;
      pt_rear->pressure = pt_insert->pressure;
    }
  return error;
}

/* Map t_o domain into water content in 1D space discretization, water_content[num_elements] and effective_porosity are output values.
 *
 * Parameters:
 * domain             - A pointer to the t_o_domain struct.
 * num_elements       - Number of element in 1d soil column.
 * soil_depth_z       - A pointer to 1d array of size num_elements contains depth of each element's lower bound, in unit of [meters], positive downward.
 *
 * Output:
 * water_content      - A pointer to 1d array of size num_elements contains water content of each element.
 * pressure_head      - A pointer to 1d array of size num_elements contains pressure head of each element.
 * effective_porosity - Scalar passed by reference, unitless.
 */
int get_t_o_domain_profile(t_o_domain* domain, int num_elements, double* soil_depth_z, 
                           double* water_content,  double* pressure_head, double* effective_porosity)
{
  int error = FALSE;
  int ii, jj, jj_start;
  t_o_profile_link *head_pt = NULL, *tmp_pt = NULL, tmp;
  
  // Initialize head pointer node, fake node.
  tmp.depth    = domain->layer_top_depth;
  tmp.theta    = domain->parameters->bin_water_content[domain->parameters->num_bins];
  tmp.pressure = 0.0;
  tmp.next     = NULL;
  head_pt      = &tmp;
  
  assert(NULL != domain);
  
  // Find effective porosity.
  *effective_porosity = domain->parameters->bin_water_content[domain->parameters->num_bins];
  
  // Step 1, map the to_domain water content into link list t_o_profile_link_link.
  // Loop over bins, start from bin 2, since bin 1 is always saturated.
  for (ii = 2; ii <= domain->parameters->num_bins; ii++)
     {
       // Surface front.
       if (domain->surface_front[ii] > domain->layer_top_depth && domain->surface_front[ii] <= domain->layer_bottom_depth) // Add "<=" 10/09/14.
         {
           error = insert_t_o_profile_link(&head_pt, domain->surface_front[ii], domain->parameters->bin_water_content[ii], 
                                                                               -domain->parameters->bin_capillary_suction[ii], domain->yes_groundwater);
         }
        
       // Slug.
       slug* temp_slug = domain->top_slug[ii];
       while (NULL != temp_slug)
         {
           // slug top.
           error = insert_t_o_profile_link(&head_pt, temp_slug->top, domain->parameters->bin_water_content[ii - 1], 
                                                                    -domain->parameters->bin_capillary_suction[ii - 1], domain->yes_groundwater);
           
           // slug bottom.
           error = insert_t_o_profile_link(&head_pt, temp_slug->bot, domain->parameters->bin_water_content[ii], 
                                                                   -domain->parameters->bin_capillary_suction[ii], domain->yes_groundwater);
           
           temp_slug = temp_slug->next;
         }
        
       // groundwater front. 
       if (domain->yes_groundwater)
         {
           if (domain->groundwater_front[ii] > domain->layer_top_depth)
             {
              error = insert_t_o_profile_link(&head_pt, domain->groundwater_front[ii], domain->parameters->bin_water_content[ii - 1], 
                                                                                 -domain->parameters->bin_capillary_suction[ii - 1], domain->yes_groundwater);
             }
         }
     } // End of bin loop.
  
  // Insert fully saturated bin;
  if (domain->yes_groundwater)
    {
      error = insert_t_o_profile_link(&head_pt, domain->layer_bottom_depth, domain->parameters->bin_water_content[domain->parameters->num_bins], 
                                                          0.0, domain->yes_groundwater);
    }
  else
    {
      int first_bin = 2; // The leftmost bin that is not completely full of water.
 
      while (first_bin <= domain->parameters->num_bins && domain->parameters->bin_water_content[first_bin] <= domain->initial_water_content)
       {
         first_bin++;
       }
      error = insert_t_o_profile_link(&head_pt, domain->layer_bottom_depth, domain->parameters->bin_water_content[first_bin - 1], 
                                                                         -domain->parameters->bin_capillary_suction[first_bin - 1], domain->yes_groundwater);
    }

  // Step 2, fill 1D array water content using t_o_profile_link_link.
  tmp_pt   = head_pt;
  jj_start = 1;
  while (NULL != tmp_pt)
    {
      for (jj = jj_start; jj <= num_elements; jj++)
         {
           if (soil_depth_z[jj] <= tmp_pt->depth)
             {
               water_content[jj] = tmp_pt->theta;
               pressure_head[jj] = tmp_pt->pressure;
             }
           else
             {
               jj_start = jj;
               break;
             }
         }
      tmp_pt = tmp_pt->next;
    } 
    
  // Set head_pt to head_pt->next, as the first node is not malloc. Then free link-list.
  head_pt = head_pt->next;
  while (NULL != head_pt)
    {
       tmp_pt  = head_pt;
       head_pt = head_pt->next;
       free(tmp_pt);
    }
  
  return error;
} // End of get_t_o_domain_profile().

// #####################################################################################################################################################
// #####################################################################################################################################################
int main(void)
{
    /*******************/
   /* Initialize X11. */
  /*******************/
#ifdef YES_PLOT
  Display* the_display = XOpenDisplay(NULL);

  if (NULL == the_display)
    {
      fprintf(stderr, "ERROR: Could not open X display.\n");
      exit(1);
    }
  int the_screen = DefaultScreen(the_display);

  Window the_window = XCreateSimpleWindow(the_display, RootWindow(the_display, the_screen), 0, 0, WINDOW_WIDTH, WINDOW_HEIGHT,
                                        1, BlackPixel(the_display, the_screen), WhitePixel(the_display, the_screen));

  XSelectInput(the_display, the_window, ExposureMask);
  XMapWindow(the_display, the_window);

  Pixmap the_pixmap = XCreatePixmap(the_display, the_window, WINDOW_WIDTH, WINDOW_HEIGHT, DefaultDepth(the_display, the_screen));

  // We have to wait for at least one expose event before it will display anything.
  XEvent the_event;

  XNextEvent(the_display, &the_event);
#endif

    /*******************************/
   /* Create Talbot-Ogden domain. */
  /*******************************/

  t_o_parameters* parameters;
  t_o_domain*     domain;
      int    ii, jj, error = FALSE;
      int    num_bins              = 300;             // Number of bins.     ################################3
      double conductivity          = 1.0 / 360000.0; // Meters per second.
      double porosity              = 0.4;     // Unitless fraction.
      double residual_saturation   = 0.07;    // Unitless fraction.
      int    van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      double vg_alpha              = 3.6;      // One over meters.
      double vg_n                  = 1.56;     // Unitless.
      double bc_lambda             = 5.5;    // Unitless.
      double bc_psib               = 0.37;     // Meters.
      double layer_top_depth       = 0.0;      // Meters.
      double layer_bottom_depth    = 1.0;      // Meters.     ##################################
      int    yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      int    yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      double initial_water_content = 0.08;      // Unitless fraction.
      double water_table           = layer_bottom_depth;      // Meters.
      
      double current_time          = 0.0;                                                     // Current time in seconds.
      double delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      double max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
      double frate,rech,runoff;
      
      double infiltration_rate     = 0.0;
      double groundwater_inf_rate = 0.0;
      double rainfall_input_time[22973];
      double rainfall_input_intensity[22973];
      double potential_ET[22973];
      char string[100];
      
      double evaporated_water      = 0.0;
      double surfacewater_depth    = 0.0;                                                     // Meters.
      double total_water           = 0.0;                                                     // The should be value for total water in meters of water.
      double groundwater_recharge  = 0.0;                                                     // The total amount supplied to groundwater in meters of water.
      double accum_infil           = 0.0;
      double surfacewater_depth_old   = 0.0;
      double groundwater_recharge_old = 0.0;
      double domain_initial_water = 0.0;
      double accu_PET             = 0.0;
      double accu_rain            = 0.0;
      double PET                  = 0.0;
      double rainfall_rate        = 0.0;
      double rainfall             = 0.0;
  
  
  int    test_id  = 104;  // ###################################################################################################################
  /*
  test_id =  1, origianl panama test.
  test_id =  2, sand infiltration.
  test_id = 101 - 112, tests using 12 USAD soil type with van Genucthen parameters.
  */
  if ( 1 == test_id)
    {
      num_bins              = 300;            // Number of bins.
      conductivity          = 1.0 / 360000.0; // Meters per second.
      porosity              = 0.4;     // Unitless fraction.
      residual_saturation   = 0.027;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 3.6;      // One over meters.
      vg_n                  = 1.56;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 2 == test_id)
    { // Sand.
      num_bins              = 300;            // Number of bins.
      conductivity          = 29.7 / 360000.0; // Meters per second.
      porosity              = 0.43;     // Unitless fraction.
      residual_saturation   = 0.045;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 14.5;      // One over meters.
      vg_n                  = 2.68;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 1.0;                                                     // The duration of the timestep in seconds.
      max_time              = 6.0 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 3 == test_id)
    { // Silt
      num_bins              = 400;            // Number of bins.
      conductivity          = 0.25 / 360000.0; // Meters per second.
      porosity              = 0.46;     // Unitless fraction.
      residual_saturation   = 0.034;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 1.6;      // One over meters.
      vg_n                  = 1.37;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 1.0;                                                     // The duration of the timestep in seconds.
      max_time              = 6.0 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 101 == test_id)
    { // Panama test using Sand. 
      conductivity          = 29.7 / 360000.0; // Meters per second.
      porosity              = 0.43;     // Unitless fraction.
      residual_saturation   = 0.045;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 14.5;      // One over meters.
      vg_n                  = 2.68;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 102 == test_id)
    {
      conductivity          = 14.5917 / 360000.0; // Meters per second.
      porosity              = 0.41;     // Unitless fraction.
      residual_saturation   = 0.057;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 12.4;      // One over meters.
      vg_n                  = 2.28;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 103 == test_id)
    {
      conductivity          = 4.42083 / 360000.0; // Meters per second.
      porosity              = 0.41;     // Unitless fraction.
      residual_saturation   = 0.065;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 7.5;      // One over meters.
      vg_n                  = 1.89;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 104 == test_id)
    {
      conductivity          = 1.04 / 360000.0; // Meters per second.
      porosity              = 0.43;     // Unitless fraction.
      residual_saturation   = 0.078;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 3.6;      // One over meters.
      vg_n                  = 1.56;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }  
  else if ( 105 == test_id)
    {
      conductivity          = 0.25 / 360000.0; // Meters per second.
      porosity              = 0.46;     // Unitless fraction.
      residual_saturation   = 0.034;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 1.6;      // One over meters.
      vg_n                  = 1.37;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 106 == test_id)
    {
      conductivity          = 0.45 / 360000.0; // Meters per second.
      porosity              = 0.45;     // Unitless fraction.
      residual_saturation   = 0.067;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 2.0;      // One over meters.
      vg_n                  = 1.41;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 107 == test_id)
    {
      conductivity          = 1.31 / 360000.0; // Meters per second.
      porosity              = 0.39;     // Unitless fraction.
      residual_saturation   = 0.1;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 5.9;      // One over meters.
      vg_n                  = 1.48;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 108 == test_id)
    {
      conductivity          = 0.26 / 360000.0; // Meters per second.
      porosity              = 0.41;     // Unitless fraction.
      residual_saturation   = 0.095;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 1.9;      // One over meters.
      vg_n                  = 1.31;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 109 == test_id)
    {
      conductivity          = 0.07 / 360000.0; // Meters per second.
      porosity              = 0.43;     // Unitless fraction.
      residual_saturation   = 0.089;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 1.0;      // One over meters.
      vg_n                  = 1.23;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 110 == test_id)
    {
      conductivity          = 0.12 / 360000.0; // Meters per second.
      porosity              = 0.38;     // Unitless fraction.
      residual_saturation   = 0.1;      // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 2.7;      // One over meters.
      vg_n                  = 1.23;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 111 == test_id)
    {
      conductivity          = 0.02 / 360000.0; // Meters per second.
      porosity              = 0.36;     // Unitless fraction.
      residual_saturation   = 0.07;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 0.5;      // One over meters.
      vg_n                  = 1.09;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 112 == test_id)
    {
      conductivity          = 0.2 / 360000.0; // Meters per second.
      porosity              = 0.38;     // Unitless fraction.
      residual_saturation   = 0.068;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 0.8;      // One over meters.
      vg_n                  = 1.09;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  
  if (t_o_parameters_alloc(&parameters, num_bins, conductivity, porosity, residual_saturation, van_genutchen,  vg_alpha, vg_n, bc_lambda, bc_psib))
    {
      fprintf(stderr, "ERROR: Could not allocate t_o_parameters.\n");
      exit(1);
    }

  if (t_o_domain_alloc(&domain, parameters, layer_top_depth, layer_bottom_depth, yes_groundwater, initial_water_content, yes_groundwater, water_table))
    {
      fprintf(stderr, "ERROR: Could not allocate t_o_domain.\n");
      exit(1);
    }

  t_o_check_invariant(domain);

#ifdef YES_PLOT
   display_domain(the_display, the_screen, the_window, the_pixmap, domain);
#endif

    /***********************/
   /* Run the simulation. */
  /***********************/
  if ( 1 == test_id || 101 <= test_id)
    {
      // Read input rainfall file.
      FILE*  rain_fptr      = NULL;
  
      if (NULL == (rain_fptr = fopen("example_rainfall_PET.txt", "r")))
        {
          printf("Error reading file example_rainfall_PET.txt \n");
          exit(1);
        }
        
      fgets(string, 100, rain_fptr);      // Ignore header.
      for (ii = 1; ii <= 22972; ii++)
         {
           fscanf(rain_fptr, "%*d %*d %*d %*d %*d %lf %lf", &rainfall_input_intensity[ii], &potential_ET[ii]);  // Original data in mm /15 min.
           rainfall_input_time[ii]       = (ii - 1) * 15.0 * 60.0;  // In s.
           rainfall_input_intensity[ii] /= 900000.0;     // In m/s. 
           rainfall_input_intensity[ii] *= 5.0;           // convert back to real values.
           potential_ET[ii]             /= 900000.0;     // In m/s.
         }
      fclose(rain_fptr);
    }
  
#define INFILTRATION_OUTPUT_FILE
#ifdef INFILTRATION_OUTPUT_FILE
  FILE*  f_fptr;                                                                          // File output of infiltration rate.
  FILE*  acc_depth_fptr;
  FILE* fptr_day_50;
  FILE* fptr_day_100;
  FILE* fptr_day_150;
  FILE* fptr_day_200;
  int get_1 = FALSE;
  int get_2 = FALSE;
  int get_3 = FALSE;
  int get_4 = FALSE;
  //double acc_depth = 0.0;
  
  if (1 == test_id  || 101 <= test_id)
    {
      if (NULL == (fptr_day_50 = fopen("panama_profile_day_50.txt", "w")) )
         {
           printf("ERROR: Could not open profile output file.\n");
           exit(1);
         }
       if (NULL == (fptr_day_100 = fopen("panama_profile_day_100.txt", "w")) )
         {
           printf("ERROR: Could not open profile output file.\n");
           exit(1);
         }
       if (NULL == (fptr_day_150 = fopen("panama_profile_day_150.txt", "w")) )
         {
           printf("ERROR: Could not open profile output file.\n");
           exit(1);
         }
       if (NULL == (fptr_day_200 = fopen("panama_profile_day_200.txt", "w")) )
         {
           printf("ERROR: Could not open profile output file.\n");
           exit(1);
         }
    }
         
  if ((f_fptr = fopen("f.out", "w")) == NULL)
    {
     fprintf(stderr, "ERROR: Could not open infiltration output file.\n");
     exit(1);
    }
  if ((acc_depth_fptr = fopen("accum_depth.out", "w")) == NULL)
    {
     fprintf(stderr, "ERROR: Could not open accumulate depth output file.\n");
     exit(1);
    }
#endif // INFILTRATION_OUTPUT_FILE
  
  int num_layers   = 1;
  int num_elements = 300;
  double** soil_depth_z;   
  double** water_content;   
  double** pressure_head;   
  double effective_porosity;
  
  error =  dtwo_alloc(&soil_depth_z,  num_layers, num_elements);
  error =  dtwo_alloc(&water_content, num_layers, num_elements);
  error =  dtwo_alloc(&pressure_head, num_layers, num_elements);
  
  for (ii = 1; ii <= num_layers ; ii++)
     {
       soil_depth_z[ii][0] = domain->layer_top_depth;
       for (jj = 1; jj <= num_elements; jj++)
          {
            soil_depth_z[ii][jj]  = soil_depth_z[ii][0] + jj * (domain->layer_bottom_depth  - domain->layer_top_depth) / num_elements;
          }
     }
  
  time_t time_start       = time(NULL); // Wall clock time.
  time_t time_end;                      // Wall clock time.
  domain_initial_water = t_o_total_water_in_domain(domain);
  accu_PET             = 0.0;
  accu_rain            = 0.0;
  PET                  = 0.0;
  rainfall_rate        = 0.0;
  rainfall             = 0.0;
  runoff               = 0.0;
  total_water          = surfacewater_depth + domain_initial_water + accu_PET + groundwater_recharge + accu_rain;
  while (current_time < max_time)
    {
      // Rainfall.  ###########################################################
      if (1 == test_id  || 101 <= test_id)
        {
          ii = 1;
          while (ii < 22972 && current_time >= rainfall_input_time[ii + 1])
            {
             ii++;
            }
          PET           = potential_ET[ii];                  // m/s.      ##################################### Set no PET.
          rainfall_rate = rainfall_input_intensity[ii];      // m/s.
          rainfall      = rainfall_rate * delta_time;        // Meters of water.
        }
      else if (2 == test_id)
        {
          PET           = 0.0;
          if ( current_time < 0.25 * ONE_HOUR || 
              (current_time > 3.0  * ONE_HOUR && current_time < 3.25  * ONE_HOUR))
            {
              rainfall_rate = 30 * 0.01 / ONE_HOUR;      // m/s.
            }
          else
            {
              rainfall_rate = 0.0;
            }
          
          rainfall      = rainfall_rate * delta_time; // Meters of water.
        }
     else if (3 == test_id)
        {
          PET           = 0.0;
          if ( current_time < 1.0 * ONE_HOUR || 
              (current_time > 3.0  * ONE_HOUR && current_time < 4.0  * ONE_HOUR))
            {
              rainfall_rate = 3.0 * 0.01 / ONE_HOUR;      // m/s.
            }
          else
            {
              rainfall_rate = 0.0;
            }
          
          rainfall      = rainfall_rate * delta_time; // Meters of water.
        }
      surfacewater_depth  += rainfall;
      total_water         += rainfall;
      accu_rain           += rainfall;
      accu_PET            += PET * delta_time;
      
      surfacewater_depth_old   = surfacewater_depth;
      groundwater_recharge_old = groundwater_recharge;
      // Moving water table
     // water_table = 1.0 - 0.7 * fabs(sin(2.0 * 3.14 * current_time / (24.0 * 3600.0)));  //
      
      if (t_o_timestep(domain, delta_time, surfacewater_depth, &surfacewater_depth, water_table, &groundwater_recharge))
        {
          fprintf(stderr, "ERROR: t_o_timestep returned error.\n");
          exit(1);
        }
      
      // FIXME, Add ET as seperate function, Dec. 10, 2014. Better put them inside t_o_timestep.
      if (1 == test_id)
        {
          // t_o_ET(domain, delta_time, root_depth, PET, field_capacity, wilting_ponint, 
                    // use_feddes, field_capacity_suction, wilting_ponint_suction, &surfacewater_depth, &evaporated_water);
          if (t_o_ET(domain, delta_time, 0.5,       PET, 0.32,           0.03,           TRUE, 0.27, 1527.68, &surfacewater_depth, &evaporated_water))
            {
              fprintf(stderr, "ERROR: t_o_timestep returned error.\n");
              exit(1);
            }
        }// End of ET for test 1.
      else if (101 <= test_id)
        {
          if (t_o_ET(domain, delta_time, 0.5,       PET, 0.32,           0.03,           TRUE, 0.3, 150.0, &surfacewater_depth, &evaporated_water))
            {
              fprintf(stderr, "ERROR: t_o_timestep returned error.\n");
              exit(1);
            }
        }
      
      accum_infil                += surfacewater_depth_old - surfacewater_depth;
#ifdef INFILTRATION_OUTPUT_FILE
      if (delta_time > (int)current_time % 60) // output every 10 min (600 s). ######################################################
      {
        infiltration_rate    = (surfacewater_depth_old - surfacewater_depth) / delta_time * 100.0 * ONE_HOUR; // cm/hr.
        groundwater_inf_rate = (groundwater_recharge - groundwater_recharge_old) / delta_time * 100 * ONE_HOUR; // cm/hr.
        fprintf(f_fptr, "%lf %lf %lf %lf %lf %lf %lf %lf \n", current_time, rainfall_rate * 360000.0, accu_rain * 100, infiltration_rate, accum_infil * 100.0,
                                                            groundwater_inf_rate, groundwater_recharge * 100.0, evaporated_water * 100.0);
      }
#else
      surfacewater_depth_old = surfacewater_depth_old; // Prevent unused variable warning.
#endif // INFILTRATION_OUTPUT_FILE

      frate = (surfacewater_depth_old - surfacewater_depth) / delta_time * 100.0 * ONE_HOUR; // cm/hr.
      rech  = (groundwater_recharge - groundwater_recharge_old) / delta_time * 100.0 * ONE_HOUR; // cm/hr.
      //printf("%lf %lf %lf %lf \n", current_time/86400, rainfall_rate * 360000.0, frate, rech);
  
      current_time += delta_time;
#ifdef YES_PLOT
      display_domain(the_display, the_screen, the_window, the_pixmap, domain);     
#endif
      //usleep(1000);
       //printf("simulation time = %lf hr.\n", current_time / ONE_HOUR);
      //printf("Total water should be %0.20lf, is %0.20lf\n", total_water, surfacewater_depth + groundwater_recharge + t_o_total_water_in_domain(domain));

      // FIXME,  one simple way to implement no flow lower boundary, if groundwater_recharge is positive, put them back into domain.
      // What about negative groundwater_recharge???
      /*if (1 == test_id)
        {
          if (groundwater_recharge > 0.0 ) //&& current_time < 28 * 24 * 3600.0)
           {
             //t_o_add_groundwater(domain, &groundwater_recharge);
           }
        }*/
      
      if(yes_runoff)  // new FLO 21 Dec. 2014 
        {
          runoff            += surfacewater_depth;
          surfacewater_depth = 0.0;
        }
     
     assert(epsilon_equal(total_water, evaporated_water + surfacewater_depth + groundwater_recharge + t_o_total_water_in_domain(domain) + runoff));
      
      fprintf(acc_depth_fptr, "%lf %lf\n", current_time, runoff * 100.0); // Time in s, runoff in cm.
      //fprintf(acc_depth_fptr,"%lf %lf\n",current_time/86400.0, t_o_total_water_in_domain(domain)); 
      
      if ( 1 == test_id || 101 <= test_id)
        {
          if (abs(current_time - 50 * ONE_DAY) < delta_time && !get_1)
             {
               get_t_o_domain_profile(domain, num_elements, soil_depth_z[1], water_content[1], pressure_head[1], &effective_porosity);
               for (jj = 1; jj <= num_elements; jj++)
                  {
                    fprintf(fptr_day_50, "%lf %lf %lf\n", soil_depth_z[1][jj], water_content[1][jj], pressure_head[1][jj]);
                  }
               get_1 = TRUE;
               //getchar();
             }
           else if (abs(current_time - 100 * ONE_DAY) < delta_time && !get_2)
             {
               get_t_o_domain_profile(domain, num_elements, soil_depth_z[1], water_content[1], pressure_head[1], &effective_porosity);
               for (jj = 1; jj <= num_elements; jj++)
                  {
                    fprintf(fptr_day_100, "%lf %lf %lf\n", soil_depth_z[1][jj], water_content[1][jj], pressure_head[1][jj]);
                  }
               get_2 = TRUE;
               //getchar(); // REMOVEM E
             }
           else if (abs(current_time - 150 * ONE_DAY) < delta_time && !get_3)
             {
               get_t_o_domain_profile(domain, num_elements, soil_depth_z[1], water_content[1], pressure_head[1], &effective_porosity);
               for (jj = 1; jj <= num_elements; jj++)
                  {
                    fprintf(fptr_day_150 ,"%lf %lf %lf\n", soil_depth_z[1][jj], water_content[1][jj], pressure_head[1][jj]);
                  }
               get_3 = TRUE;
               //getchar();
             }
           else if (abs(current_time - 200 * ONE_DAY) < delta_time && !get_4)
             {
               get_t_o_domain_profile(domain, num_elements, soil_depth_z[1], water_content[1], pressure_head[1], &effective_porosity);
               for (jj = 1; jj <= num_elements; jj++)
                  {
                    fprintf(fptr_day_200, "%lf %lf %lf\n", soil_depth_z[1][jj], water_content[1][jj], pressure_head[1][jj]);
                  }
               get_4 = TRUE;
               //getchar();
             } 
        } // End of ouput for test_id = 1;
       
    } // End of time loop.
  time_end = time(NULL);
  printf("\nTotal simulation time  = %lf hours\n", max_time / ONE_HOUR);
  printf("\nElapsed time = %lf seconds \n \n",(double)(time_end - time_start));
  printf("Mass balance info: \n");
  printf("Initial water in domain = %lf mm \n", domain_initial_water * 1000);
  printf("Accumulated rainfall    = %lf mm \n", accu_rain * 1000);
  printf("Accumulated PET         = %lf mm \n", accu_PET * 1000);
  printf("Actual ET               = %lf mm \n", evaporated_water * 1000);
  printf("Groundwater recharge    = %lf mm \n", groundwater_recharge * 1000);
  printf("Final water in domain   = %lf mm \n", t_o_total_water_in_domain(domain) * 1000);
  printf("Final surface water     = %lf mm \n", surfacewater_depth);
  printf("Total surface runoff    = %lf mm \n", runoff*1000.0);
  printf("Mass error              = %8.5e mm \n", (domain_initial_water + accu_rain - evaporated_water - groundwater_recharge - 
                                                 t_o_total_water_in_domain(domain) - surfacewater_depth-runoff) * 1000);
  
    /*************/
   /* Clean up. */
  /*************/
  if (1 == test_id)
    {
      fclose(fptr_day_50);
      fclose(fptr_day_100);
      fclose(fptr_day_150);
      fclose(fptr_day_200);
    }
  t_o_domain_dealloc(&domain);
  t_o_parameters_dealloc(&parameters);
#ifdef YES_PLOT
   XCloseDisplay(the_display);
#endif
#ifdef INFILTRATION_OUTPUT_FILE
  fclose(f_fptr);
  fclose(acc_depth_fptr);
#endif // INFILTRATION_OUTPUT_FILE

  return 0;
}
