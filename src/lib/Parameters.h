
/**
 * \file      Parameters.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Parameters class declaration
 */

/****************************************************************************
 * Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * Web: https://github.com/charlesrocabert/Evo2Sim
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#ifndef __EVOEVO__Parameters__
#define __EVOEVO__Parameters__

#include <iostream>
#include <cstring>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <zlib.h>
#include <assert.h>

#include "Macros.h"
#include "Enums.h"
#include "Structs.h"
#include "Prng.h"


class Parameters
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  Parameters( void );
  Parameters( size_t backup_time );
  Parameters( const Parameters& parameters );
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~Parameters( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  
  /*------------------------------------------------------------------ prng seed */
  
  inline unsigned long int get_seed( void ) const;
  
  /*------------------------------------------------------------------ parallel computing */
  
  inline bool get_parallel_computing( void ) const;
  
  /*------------------------------------------------------------------ simulation schemes */
  
  inline bool         get_energy_costs_scheme( void ) const;
  inline bool         get_membrane_permeability_scheme( void ) const;
  inline bool         get_metabolic_inheritance_scheme( void ) const;
  inline bool         get_enzymatic_inheritance_scheme( void ) const;
  inline bool         get_co_enzyme_activity_scheme( void ) const;
  inline score_scheme get_score_scheme( void ) const;
  inline double       get_selection_threshold( void ) const;
  
  /*------------------------------------------------------------------ space */
  
  inline size_t get_width( void ) const;
  inline size_t get_height( void ) const;
  
  /*------------------------------------------------------------------ output */
  
  inline size_t get_simulation_backup_step( void ) const;
  inline size_t get_figures_generation_step( void ) const;
  
  /*------------------------------------------------------------------ genome */
  
  inline bool get_load_genome_from_file( void ) const;
  
  inline variable_range* get_metabolite_tag_initial_range( void );
  inline variable_range* get_binding_site_tag_initial_range( void );
  inline variable_range* get_co_enzyme_tag_initial_range( void );
  inline variable_range* get_transcription_factor_tag_initial_range( void );
  
  inline size_t get_transcription_factor_binding_window( void ) const;
  
  inline size_t get_initial_number_of_NC_units( void ) const;
  inline size_t get_initial_number_of_E_units( void ) const;
  inline size_t get_initial_number_of_TF_units( void ) const;
  inline size_t get_initial_number_of_BS_units( void ) const;
  inline size_t get_initial_number_of_P_units( void ) const;
  
  inline double get_point_mutation_rate( void ) const;
  inline double get_duplication_rate( void ) const;
  inline double get_deletion_rate( void ) const;
  inline double get_translocation_rate( void ) const;
  inline double get_inversion_rate( void ) const;
  inline double get_transition_rate( void ) const;
  inline double get_breakpoint_rate( void ) const;
  
  inline double get_substrate_tag_mutation_size( void ) const;
  inline double get_product_tag_mutation_size( void ) const;
  inline double get_kcat_mutation_size( void ) const;
  inline double get_kcat_km_ratio_mutation_size( void ) const;
  inline double get_binding_site_tag_mutation_size( void ) const;
  inline double get_co_enzyme_tag_mutation_size( void ) const;
  inline double get_transcription_factor_tag_mutation_size( void ) const;
  inline double get_basal_expression_level_mutation_size( void ) const;
  
  inline double get_mutation_of_mutation_rates( void ) const;
  
  /*------------------------------------------------------------------ genetic regulation network */
  
  inline double get_genetic_regulation_network_timestep( void ) const;
  
  inline double get_hill_function_theta( void ) const;
  inline double get_hill_function_n( void ) const;
  inline double get_protein_degradation_rate( void ) const;
  
  /*------------------------------------------------------------------ metabolic network */
  
  inline double get_metabolism_timestep( void ) const;
  
  inline double get_essential_metabolites_toxicity_threshold( void );
  inline double get_non_essential_metabolites_toxicity_threshold( void );
  
  inline double get_initial_metabolites_amount_in_cells( void ) const;
  
  inline size_t get_maximum_reaction_size( void ) const;
  
  /*------------------------------------------------------------------ energy */
  
  inline double get_energy_transcription_cost( void ) const;
  inline double get_energy_degradation_cost( void ) const;
  inline double get_energy_enzymatic_cost( void ) const;
  inline double get_energy_pumping_cost( void ) const;
  
  inline double get_energy_dissipation_rate( void ) const;
  
  inline double get_energy_toxicity_threshold( void ) const;
  
  inline double get_initial_energy_amount_in_cells( void ) const;
  
  /*------------------------------------------------------------------ cell */
  
  inline double get_membrane_permeability( void ) const;
  
  /*------------------------------------------------------------------ population */
  
  inline double get_death_probability( void ) const;
  inline double get_migration_rate( void ) const;
  inline double get_hgt_rate( void ) const;
  
  /*------------------------------------------------------------------ environment */
  
  inline environment_properties* get_environment_properties( void );
  
  /*------------------------------------------------------------------ prime numbers */
  
  inline int* get_prime_numbers( void );
  
  /*------------------------------------------------------------------ prngs */
  
  inline Prng* get_simulation_prng( void );
  inline Prng* get_population_prng( size_t pos );
  inline Prng* get_environment_prng( void );
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  
  /*------------------------------------------------------------------ prng seed */
  
  inline void set_seed( unsigned long int seed );
  
  /*------------------------------------------------------------------ parallel computing */
  
  inline void set_parallel_computing( bool parallel_computing );
  
  /*------------------------------------------------------------------ simulation schemes */
  
  inline void set_energy_costs_scheme( bool energy_costs_scheme );
  inline void set_membrane_permeability_scheme( bool membrane_permeability_scheme );
  inline void set_metabolic_inheritance_scheme( bool metabolic_inheritance_scheme );
  inline void set_enzymatic_inheritance_scheme( bool enzymatic_inheritance_scheme );
  inline void set_co_enzyme_activity_scheme( bool co_enzyme_activity_scheme );
  inline void set_score_scheme( score_scheme scheme );
  inline void set_selection_threshold( double selection_threshold );
  
  /*------------------------------------------------------------------ space */
  
  inline void set_width( size_t width );
  inline void set_height( size_t height );
  
  /*------------------------------------------------------------------ output */
  
  inline void set_simulation_backup_step( size_t step );
  inline void set_figures_generation_step( size_t step );
  
  /*------------------------------------------------------------------ genome */
  
  inline void set_load_genome_from_file( bool choice );
  
  inline void set_metabolite_tag_initial_range( const variable_range* range );
  inline void set_binding_site_tag_initial_range( const variable_range* range );
  inline void set_co_enzyme_tag_initial_range( const variable_range* range );
  inline void set_transcription_factor_tag_initial_range( const variable_range* range );
  
  inline void set_transcription_factor_binding_window( size_t window );
  
  inline void set_initial_number_of_NC_units( size_t number_of_NC_units );
  inline void set_initial_number_of_E_units( size_t number_of_E_units );
  inline void set_initial_number_of_TF_units( size_t number_of_TF_units );
  inline void set_initial_number_of_BS_units( size_t number_of_BS_units );
  inline void set_initial_number_of_P_units( size_t number_of_P_units );
  
  inline void set_point_mutation_rate( double rate );
  inline void set_duplication_rate( double rate );
  inline void set_deletion_rate( double rate );
  inline void set_translocation_rate( double rate );
  inline void set_inversion_rate( double rate );
  inline void set_transition_rate( double rate );
  inline void set_breakpoint_rate( double rate );
  
  inline void set_substrate_tag_mutation_size( double size );
  inline void set_product_tag_mutation_size( double size );
  inline void set_kcat_mutation_size( double size );
  inline void set_kcat_km_ratio_mutation_size( double size );
  inline void set_binding_site_tag_mutation_size( double size );
  inline void set_co_enzyme_tag_mutation_size( double size );
  inline void set_transcription_factor_tag_mutation_size( double size );
  inline void set_basal_expression_level_mutation_size( double size );
  
  inline void set_mutation_of_mutation_rates( double rate );
  
  /*------------------------------------------------------------------ genetic regulation network */
  
  inline void set_genetic_regulation_network_timestep( double genetic_regulation_network_timestep );
  
  inline void set_hill_function_theta( double hill_function_theta );
  inline void set_hill_function_n( double hill_function_n );
  inline void set_protein_degradation_rate( double protein_degradation_rate );
  
  /*------------------------------------------------------------------ metabolic network */
  
  inline void set_metabolism_timestep( double metabolism_timestep );
  
  inline void set_essential_metabolites_toxicity_threshold( double toxicity_threshold );
  inline void set_non_essential_metabolites_toxicity_threshold( double toxicity_threshold );
  
  inline void set_initial_metabolites_amount_in_cells( double metabolites_amount_in_cells );
  
  inline void set_maximum_reaction_size( size_t size );
  
  /*------------------------------------------------------------------ energy */
  
  inline void set_energy_transcription_cost( double cost );
  inline void set_energy_degradation_cost( double cost );
  inline void set_energy_enzymatic_cost( double cost );
  inline void set_energy_pumping_cost( double cost );
  
  inline void set_energy_dissipation_rate( double rate );
  
  inline void set_energy_toxicity_threshold( double threshold );
  
  inline void set_initial_energy_amount_in_cells( double amount );
  
  /*------------------------------------------------------------------ cell */
  
  inline void set_membrane_permeability( double permeability );
  
  /*------------------------------------------------------------------ population */
  
  inline void set_death_probability( double death_probability );
  inline void set_migration_rate( double migration_rate );
  inline void set_hgt_rate( double hgt_rate );
  
  /*------------------------------------------------------------------ environment */
  
  inline void set_environment_properties( const environment_properties* properties );
  
  /*------------------------------------------------------------------ prngs */
  
  inline void set_simulation_prng( Prng* prng );
  inline void set_population_prng( Prng* prng, size_t pos );
  inline void set_environment_prng( Prng* prng );
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void load_parameters_from_file( std::string filename );
  void save( size_t backup_time );
  void write( std::string filename );
  
  int    draw_initial_metabolite_tag( void );
  int    draw_initial_binding_site_tag( void );
  int    draw_initial_co_enzyme_tag( void );
  int    draw_initial_transcription_factor_tag( void );
  int    draw_environment_number_of_species( void );
  int    draw_environment_species_tag( void );
  double draw_environment_concentration( void );
  
  void draw_random_genetic_unit( genetic_unit& unit, genetic_unit_type type );
  void load_genome_from_file( std::vector<genetic_unit>& sequence );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  bool parse_line( std::vector<std::string>* words, std::string line );
  void build_prime_numbers_list( size_t maximum );
  void initialize_prngs( void );
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  
  /*------------------------------------------------------------------ prng seed */
  
  unsigned long int _seed; /*!< Seed of the prng */
  
  /*------------------------------------------------------------------ parallel computing */
  
  bool _parallel_computing; /*!< Indicates if the parallel computing is activated */
  
  /*------------------------------------------------------------------ simulation schemes */
  
  bool         _energy_costs_scheme;          /*!< Indicates if energy costs are activated              */
  bool         _membrane_permeability_scheme; /*!< Indicates if membrane permeability is activated      */
  bool         _metabolic_inheritance_scheme; /*!< Indicates if metabolic inheritance is activated      */
  bool         _enzymatic_inheritance_scheme; /*!< Indicates if enzymatic inheritance is activated      */
  bool         _co_enzyme_activity_scheme;    /*!< Indicates if co-enzyme activity is activated         */
  score_scheme _score_scheme;                 /*!< Score scheme used to compute the score               */
  double       _selection_threshold;          /*!< Selection threshold when selection scheme is fitprop */
  
  /*------------------------------------------------------------------ space */
  
  size_t _width;  /*!< Spatial grid width  */
  size_t _height; /*!< Spatial grid height */
  
  /*------------------------------------------------------------------ output */
  
  size_t _simulation_backup_step;  /*!< Step of simulation backups */
  size_t _figures_generation_step; /*!< Step of figures generation */
  
  /*------------------------------------------------------------------ genome */
  
  bool _load_genome_from_file; /*!< Load the initial genome from a text file */
  
  variable_range _metabolite_tag_initial_range;           /*!< Metabolite tag initial range           */
  variable_range _binding_site_tag_initial_range;         /*!< Binding site tag initial range         */
  variable_range _co_enzyme_tag_initial_range;            /*!< Co-enzyme tag initial range            */
  variable_range _transcription_factor_tag_initial_range; /*!< Transcription factor tag initial range */
  
  size_t         _transcription_factor_binding_window;    /*!< Binding window of a TF on a BS         */
  
  size_t         _initial_number_of_NC_units;             /*!< Initial number of NC units             */
  size_t         _initial_number_of_E_units;              /*!< Initial number of E units              */
  size_t         _initial_number_of_TF_units;             /*!< Initial number of TF units             */
  size_t         _initial_number_of_BS_units;             /*!< Initial number of BS units             */
  size_t         _initial_number_of_P_units;              /*!< Initial number of P units              */
  
  double         _point_mutation_rate;                    /*!< Point mutation rate                    */
  double         _duplication_rate;                       /*!< Duplication rate                       */
  double         _deletion_rate;                          /*!< Deletion rate                          */
  double         _translocation_rate;                     /*!< Translocation rate                     */
  double         _inversion_rate;                         /*!< Inversion rate                         */
  double         _transition_rate;                        /*!< Transition rate                        */
  double         _breakpoint_rate;                        /*!< Breakpoint rate                        */
  
  double         _substrate_tag_mutation_size;            /*!< Substrate tag mutation size            */
  double         _product_tag_mutation_size;              /*!< Product tag mutation size              */
  double         _kcat_mutation_size;                     /*!< Kcat constant mutation size            */
  double         _kcat_km_ratio_mutation_size;            /*!< Kcat/Km ratio constant mutation size   */
  double         _binding_site_tag_mutation_size;         /*!< Binding site tag mutation size         */
  double         _co_enzyme_tag_mutation_size;            /*!< Co-enzyme tag mutation size            */
  double         _transcription_factor_tag_mutation_size; /*!< Transcription factor tag mutation size */
  double         _basal_expression_level_mutation_size;   /*!< Basal expression level mutation size   */
  
  double         _mutation_of_mutation_rates;             /*!< Mutation of mutation rates             */
  
  /*------------------------------------------------------------------ genetic regulation network */
  
  double _genetic_regulation_network_timestep; /*!< Genetic regulation network ODE timestep */
  
  double _hill_function_theta;                 /*!< Parameter theta of the Hill function    */
  double _hill_function_n;                     /*!< Parameter n of the Hill function        */
  double _protein_degradation_rate;            /*!< Protein degradation rate                */
  
  /*------------------------------------------------------------------ metabolic network */
  
  double _metabolism_timestep;                          /*!< Metabolism ODE timestep                      */
  
  double _essential_metabolites_toxicity_threshold;     /*!< Essential metabolites toxicity threshold     */
  double _non_essential_metabolites_toxicity_threshold; /*!< Non essential metabolites toxicity threshold */
  
  double _initial_metabolites_amount_in_cells;          /*!< Initial metabolites amount in cells          */
  
  size_t _maximum_reaction_size;                        /*!< Maximum reaction size possible               */
  
  /*------------------------------------------------------------------ energy */
  
  double _energy_transcription_cost;      /*!< Energetical cost of transcription       */
  double _energy_degradation_cost;        /*!< Energetical cost of protein degradation */
  double _energy_enzymatic_cost;          /*!< Energetical cost of enzymatic activity  */
  double _energy_pumping_cost;            /*!< Energetical cost of pumping activity    */
  
  double _energy_dissipation_rate;        /*!< Energy dissipation rate                 */
  
  double _energy_toxicity_threshold;      /*!< Energy toxicity threshold               */
  
  double _initial_energy_amount_in_cells; /*!< Initial energy amount in cells          */
                                           
  /*------------------------------------------------------------------ cell */
  
  double _membrane_permeability; /*!< Membrane permeability (diffusion coefficient) */
  
  /*------------------------------------------------------------------ population */
  
  double _death_probability; /*!< Death probability */
  double _migration_rate;    /*!< Migration rate    */
  double _hgt_rate;          /*!< HGT rate          */
  
  /*------------------------------------------------------------------ environment */
  
  environment_properties _environment_properties; /*!< Environment properties */
  
  /*------------------------------------------------------------------ prime numbers */
  
  int* _prime_numbers; /*!< List of prime numbers */
  
  /*------------------------------------------------------------------ prngs */
  
  Prng*  _simulation_prng;  /*!< Simulation prng          */
  Prng** _population_prng;  /*!< List of population prngs */
  Prng*  _environment_prng; /*!< Environment prng         */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/*------------------------------------------------------------------ prng seed */

/**
 * \brief    Get PRNG seed
 * \details  Return the seed of the prng
 * \param    void
 * \return   \e unsigned long int
 */
inline unsigned long int Parameters::get_seed( void ) const
{
  return _seed;
}

/*------------------------------------------------------------------ parallel computing */

/**
 * \brief    Get parallel computing
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Parameters::get_parallel_computing( void ) const
{
  return _parallel_computing;
}

/*------------------------------------------------------------------ simulation schemes */

/**
 * \brief    Get the energy costs scheme
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Parameters::get_energy_costs_scheme( void ) const
{
  return _energy_costs_scheme;
}

/**
 * \brief    Get the membrane permeability scheme
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Parameters::get_membrane_permeability_scheme( void ) const
{
  return _membrane_permeability_scheme;
}

/**
 * \brief    Get the metabolic inheritance scheme
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Parameters::get_metabolic_inheritance_scheme( void ) const
{
  return _metabolic_inheritance_scheme;
}

/**
 * \brief    Get the enzymatic inheritance scheme
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Parameters::get_enzymatic_inheritance_scheme( void ) const
{
  return _enzymatic_inheritance_scheme;
}

/**
 * \brief    Get the co-enzyme activity scheme
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Parameters::get_co_enzyme_activity_scheme( void ) const
{
  return _co_enzyme_activity_scheme;
}

/**
 * \brief    Get the score scheme
 * \details  --
 * \param    void
 * \return   \e score_scheme
 */
inline score_scheme Parameters::get_score_scheme( void ) const
{
  return _score_scheme;
}

/**
 * \brief    Get the selection threshold
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_selection_threshold( void ) const
{
  return _selection_threshold;
}

/*------------------------------------------------------------------ space */

/**
 * \brief    Get population's grid width
 * \details  Return the population's grid width
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_width( void ) const
{
  return _width;
}

/**
 * \brief    Get population's grid height
 * \details  Return the population's grid height
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_height( void ) const
{
  return _height;
}

/*------------------------------------------------------------------ output */

/**
 * \brief    Get simulation backup step
 * \details  Return the step of each simulation backup
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_simulation_backup_step( void ) const
{
  return _simulation_backup_step;
}

/**
 * \brief    Get figures generation step
 * \details  Return the step of each figures generation
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_figures_generation_step( void ) const
{
  return _figures_generation_step;
}

/*------------------------------------------------------------------ genome */

/**
 * \brief    Get the load genome from file boolean
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Parameters::get_load_genome_from_file( void ) const
{
  return _load_genome_from_file;
}

/**
 * \brief    Get metabolite tag initial range
 * \details  --
 * \param    void
 * \return   \e variable_range*
 */
inline variable_range* Parameters::get_metabolite_tag_initial_range( void )
{
  return &_metabolite_tag_initial_range;
}

/**
 * \brief    Binding site tag initial range
 * \details  --
 * \param    void
 * \return   \e variable_range*
 */
inline variable_range* Parameters::get_binding_site_tag_initial_range( void )
{
  return &_binding_site_tag_initial_range;
}

/**
 * \brief    Co-enzyme tag initial range
 * \details  --
 * \param    void
 * \return   \e variable_range*
 */
inline variable_range* Parameters::get_co_enzyme_tag_initial_range( void )
{
  return &_co_enzyme_tag_initial_range;
}

/**
 * \brief    Get transcription factor tag initial range
 * \details  --
 * \param    void
 * \return   \e variable_range*
 */
inline variable_range* Parameters::get_transcription_factor_tag_initial_range( void )
{
  return &_transcription_factor_tag_initial_range;
}

/**
 * \brief    Get binding window of a TF on a BS
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_transcription_factor_binding_window( void ) const
{
  return _transcription_factor_binding_window;
}

/**
 * \brief    Get initial number of non coding units (NC)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_initial_number_of_NC_units( void ) const
{
  return _initial_number_of_NC_units;
}

/**
 * \brief    Get initial number of enzyme units (E)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_initial_number_of_E_units( void ) const
{
  return _initial_number_of_E_units;
}

/**
 * \brief    Get initial number of transcription factor units (TF)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_initial_number_of_TF_units( void ) const
{
  return _initial_number_of_TF_units;
}

/**
 * \brief    Get initial number of binding site units (BS)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_initial_number_of_BS_units( void ) const
{
  return _initial_number_of_BS_units;
}

/**
 * \brief    Get initial number of promoter units (P)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_initial_number_of_P_units( void ) const
{
  return _initial_number_of_P_units;
}

/**
 * \brief    Get point mutation rate
 * \details  Return point mutation rate used to initialize individuals
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_point_mutation_rate( void ) const
{
  return _point_mutation_rate;
}

/**
 * \brief    Get duplication rate
 * \details  Return duplication rate used to initialize individuals
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_duplication_rate( void ) const
{
  return _duplication_rate;
}

/**
 * \brief    Get deletion rate
 * \details  Return deletion rate used to initialize individuals
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_deletion_rate( void ) const
{
  return _deletion_rate;
}

/**
 * \brief    Get translocation rate
 * \details  Return translocation rate used to initialize individuals
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_translocation_rate( void ) const
{
  return _translocation_rate;
}

/**
 * \brief    Get inversion rate
 * \details  Return inversion rate used to initialize individuals
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_inversion_rate( void ) const
{
  return _inversion_rate;
}

/**
 * \brief    Get initial transition rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_transition_rate( void ) const
{
  return _transition_rate;
}

/**
 * \brief    Get initial breakpoint rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_breakpoint_rate( void ) const
{
  return _breakpoint_rate;
}

/**
 * \brief    Get substrate tag mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_substrate_tag_mutation_size( void ) const
{
  return _substrate_tag_mutation_size;
}

/**
 * \brief    Get product tag mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_product_tag_mutation_size( void ) const
{
  return _product_tag_mutation_size;
}

/**
 * \brief    Get Kcat constant mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_kcat_mutation_size( void ) const
{
  return _kcat_mutation_size;
}

/**
 * \brief    Get Kcat/Km ratio constant mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_kcat_km_ratio_mutation_size( void ) const
{
  return _kcat_km_ratio_mutation_size;
}

/**
 * \brief    Get binding site tag mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_binding_site_tag_mutation_size( void ) const
{
  return _binding_site_tag_mutation_size;
}

/**
 * \brief    Get co-enzyme tag mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_co_enzyme_tag_mutation_size( void ) const
{
  return _co_enzyme_tag_mutation_size;
}

/**
 * \brief    Get transcription factor tag mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_transcription_factor_tag_mutation_size( void ) const
{
  return _transcription_factor_tag_mutation_size;
}

/**
 * \brief    Get basal expression level mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_basal_expression_level_mutation_size( void ) const
{
  return _basal_expression_level_mutation_size;
}

/**
 * \brief    Get mutation of mutation rates
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_mutation_of_mutation_rates( void ) const
{
  return _mutation_of_mutation_rates;
}

/*------------------------------------------------------------------ genetic regulation network */

/**
 * \brief    Get the genetic regulation network ODE timestep
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_genetic_regulation_network_timestep( void ) const
{
  return _genetic_regulation_network_timestep;
}

/**
 * \brief    Get the parameter theta of the Hill function
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_hill_function_theta( void ) const
{
  return _hill_function_theta;
}

/**
 * \brief    Get the parameter n of the Hill function
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_hill_function_n( void ) const
{
  return _hill_function_n;
}

/**
 * \brief    Get the protein degradation rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_protein_degradation_rate( void ) const
{
  return _protein_degradation_rate;
}

/*------------------------------------------------------------------ metabolic network */

/**
 * \brief    Get metabolism ODE timestep
 * \details  Return the timestep of ODE solver per simulation timestep
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_metabolism_timestep( void ) const
{
  return _metabolism_timestep;
}

/**
 * \brief    Get essential metabolites toxicity threshold
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_essential_metabolites_toxicity_threshold( void )
{
  return _essential_metabolites_toxicity_threshold;
}

/**
 * \brief    Get non essential metabolites toxicity threshold
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_non_essential_metabolites_toxicity_threshold( void )
{
  return _non_essential_metabolites_toxicity_threshold;
}

/**
 * \brief    Get initial metabolites amount in cells
 * \details  Return the initial metabolite concentration in cells
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_initial_metabolites_amount_in_cells( void ) const
{
  return _initial_metabolites_amount_in_cells;
}

/**
 * \brief    Set the maximum reaction size
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_maximum_reaction_size( void ) const
{
  return _maximum_reaction_size;
}

/*------------------------------------------------------------------ energy */

/**
 * \brief    Get energetical cost of transcription
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_energy_transcription_cost( void ) const
{
  return _energy_transcription_cost;
}

/**
 * \brief    Get energetical cost of protein degradation
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_energy_degradation_cost( void ) const
{
  return _energy_degradation_cost;
}

/**
 * \brief    Get energetical cost of enzymatic activity
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_energy_enzymatic_cost( void ) const
{
  return _energy_enzymatic_cost;
}

/**
 * \brief    Get energetical cost of pumping activity
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_energy_pumping_cost( void ) const
{
  return _energy_pumping_cost;
}

/**
 * \brief    Get energy dissipation rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_energy_dissipation_rate( void ) const
{
  return _energy_dissipation_rate;
}

/**
 * \brief    Get energy toxicity threshold
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_energy_toxicity_threshold( void ) const
{
  return _energy_toxicity_threshold;
}

/**
 * \brief    Get initial energy amount in cells
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_initial_energy_amount_in_cells( void ) const
{
  return _initial_energy_amount_in_cells;
}

/*------------------------------------------------------------------ cell */

/**
 * \brief    Get membrane permeability
 * \details  Defines the initial permeability of the cell membrane
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_membrane_permeability( void ) const
{
  return _membrane_permeability;
}

/*------------------------------------------------------------------ population */

/**
 * \brief    Get death probability
 * \details  Return individual death probability
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_death_probability( void ) const
{
  return _death_probability;
}

/**
 * \brief    Get the migration rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_migration_rate( void ) const
{
  return _migration_rate;
}

/**
 * \brief    Get the HGT rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_hgt_rate( void ) const
{
  return _hgt_rate;
}

/*------------------------------------------------------------------ environment */

/**
 * \brief    Get the environment properties
 * \details  --
 * \param    void
 * \return   \e environment_properties*
 */
inline environment_properties* Parameters::get_environment_properties( void )
{
  return &_environment_properties;
}

/*------------------------------------------------------------------ prime numbers */

/**
 * \brief    Get list of prime numbers
 * \details  Return the first 25 prime numbers. This list is used to evaluate individuals if fitness evaluation is 'PRIME_NUMBERS'
 * \param    void
 * \return   \e int*
 */
inline int* Parameters::get_prime_numbers( void )
{
  return _prime_numbers;
}

/*------------------------------------------------------------------ prngs */

/**
 * \brief    Get the simulation prng
 * \details  --
 * \param    void
 * \return   \e Prng*
 */
inline Prng* Parameters::get_simulation_prng( void )
{
  return _simulation_prng;
}

/**
 * \brief    Get the population prng at position pos
 * \details  --
 * \param    void
 * \return   \e Prng*
 */
inline Prng* Parameters::get_population_prng( size_t pos )
{
  assert(pos < _width*_height);
  return _population_prng[pos];
}

/**
 * \brief    Get the environment prng
 * \details  --
 * \param    void
 * \return   \e Prng*
 */
inline Prng* Parameters::get_environment_prng( void )
{
  return _environment_prng;
}

/*----------------------------
 * SETTERS
 *----------------------------*/

/*------------------------------------------------------------------ prng seed */

/**
 * \brief    Set PRNG seed
 * \details  --
 * \param    unsigned long int seed
 * \return   \e void
 */
inline void Parameters::set_seed( unsigned long int seed )
{
  assert(seed > 0);
  assert(seed <= MAXIMUM_SEED);
  _seed = seed;
  initialize_prngs();
}

/*------------------------------------------------------------------ parallel computing */

/**
 * \brief    Set parallel computing
 * \details  --
 * \param    bool parallel_computing
 * \return   \e void
 */
inline void Parameters::set_parallel_computing( bool parallel_computing )
{
  _parallel_computing = parallel_computing;
}

/*------------------------------------------------------------------ simulation schemes */

/**
 * \brief    Set the energy costs scheme
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Parameters::set_energy_costs_scheme( bool energy_costs_scheme )
{
  _energy_costs_scheme = energy_costs_scheme;
}

/**
 * \brief    Set the membrane permeability scheme
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Parameters::set_membrane_permeability_scheme( bool membrane_permeability_scheme )
{
  _membrane_permeability_scheme = membrane_permeability_scheme;
}

/**
 * \brief    Set the metabolic inheritance scheme
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Parameters::set_metabolic_inheritance_scheme( bool metabolic_inheritance_scheme )
{
  _metabolic_inheritance_scheme = metabolic_inheritance_scheme;
}

/**
 * \brief    Set the enzymatic inheritance scheme
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Parameters::set_enzymatic_inheritance_scheme( bool enzymatic_inheritance_scheme )
{
  _enzymatic_inheritance_scheme = enzymatic_inheritance_scheme;
}

/**
 * \brief    Set the co-enzyme activity scheme
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Parameters::set_co_enzyme_activity_scheme( bool co_enzyme_activity_scheme )
{
  _co_enzyme_activity_scheme = co_enzyme_activity_scheme;
}

/**
 * \brief    Set the score scheme
 * \details  --
 * \param    score_scheme scheme
 * \return   \e void
 */
inline void Parameters::set_score_scheme( score_scheme scheme )
{
  _score_scheme = scheme;
}

/**
 * \brief    Set the selection threshold
 * \details  --
 * \param    double selection_threshold
 * \return   \e void
 */
inline void Parameters::set_selection_threshold( double selection_threshold )
{
  assert(selection_threshold >= 0.0);
  assert(selection_threshold <= 1.0);
  _selection_threshold = selection_threshold;
}

/*------------------------------------------------------------------ space */

/**
 * \brief    Set population's grid width
 * \details  --
 * \param    size_t width
 * \return   \e void
 */
inline void Parameters::set_width( size_t width )
{
  _width = width;
}

/**
 * \brief    Set population's grid height
 * \details  --
 * \param    size_t height
 * \return   \e void
 */
inline void Parameters::set_height( size_t height )
{
  _height = height;
}

/*------------------------------------------------------------------ output */

/**
 * \brief    Set simulation backup step
 * \details  --
 * \param    size_t step
 * \return   \e void
 */
inline void Parameters::set_simulation_backup_step( size_t step )
{
  _simulation_backup_step = step;
}

/**
 * \brief    Set figures generation step
 * \details  --
 * \param    size_t step
 * \return   \e void
 */
inline void Parameters::set_figures_generation_step( size_t step )
{
  _figures_generation_step = step;
}

/*------------------------------------------------------------------ genome */

/**
 * \brief    Set the load genome from file boolean
 * \details  --
 * \param    bool choice
 * \return   \e void
 */
inline void Parameters::set_load_genome_from_file( bool choice )
{
  _load_genome_from_file = choice;
}

/**
 * \brief    Set metabolite tag initial range
 * \details  --
 * \param    variable_range* range
 * \return   \e void
 */
inline void Parameters::set_metabolite_tag_initial_range( const variable_range* range )
{
  _metabolite_tag_initial_range.min = range->min;
  _metabolite_tag_initial_range.max = range->max;
}

/**
 * \brief    Set binding site tag initial range
 * \details  --
 * \param    variable_range* range
 * \return   \e void
 */
inline void Parameters::set_binding_site_tag_initial_range( const variable_range* range )
{
  _binding_site_tag_initial_range.min = range->min;
  _binding_site_tag_initial_range.max = range->max;
}

/**
 * \brief    Set co-enzyme tag initial range
 * \details  --
 * \param    variable_range* range
 * \return   \e void
 */
inline void Parameters::set_co_enzyme_tag_initial_range( const variable_range* range )
{
  _co_enzyme_tag_initial_range.min = range->min;
  _co_enzyme_tag_initial_range.max = range->max;
}

/**
 * \brief    Set transcription factor tag initial range
 * \details  --
 * \param    variable_range* range
 * \return   \e void
 */
inline void Parameters::set_transcription_factor_tag_initial_range( const variable_range* range )
{
  _transcription_factor_tag_initial_range.min = range->min;
  _transcription_factor_tag_initial_range.max = range->max;
}

/**
 * \brief    Set binding window of a TF on a BS
 * \details  --
 * \param    size_t window
 * \return   \e void
 */
inline void Parameters::set_transcription_factor_binding_window( size_t window )
{
  _transcription_factor_binding_window = window;
}

/**
 * \brief    Set initial number of non coding units (NC)
 * \details  --
 * \param    size_t number_of_NC_units
 * \return   \e void
 */
inline void Parameters::set_initial_number_of_NC_units( size_t number_of_NC_units )
{
  _initial_number_of_NC_units = number_of_NC_units;
}

/**
 * \brief    Set initial number of enzyme units (E)
 * \details  --
 * \param    size_t number_of_E_units
 * \return   \e void
 */
inline void Parameters::set_initial_number_of_E_units( size_t number_of_E_units )
{
  _initial_number_of_E_units = number_of_E_units;
}

/**
 * \brief    Set initial number of transcription factor units (TF)
 * \details  --
 * \param    size_t number_of_TF_units
 * \return   \e void
 */
inline void Parameters::set_initial_number_of_TF_units( size_t number_of_TF_units )
{
  _initial_number_of_TF_units = number_of_TF_units;
}

/**
 * \brief    Set initial number of binding site units (BS)
 * \details  --
 * \param    size_t number_of_BS_units
 * \return   \e void
 */
inline void Parameters::set_initial_number_of_BS_units( size_t number_of_BS_units )
{
  _initial_number_of_BS_units = number_of_BS_units;
}

/**
 * \brief    Set initial number of promoter units (P)
 * \details  --
 * \param    size_t number_of_P_units
 * \return   \e void
 */
inline void Parameters::set_initial_number_of_P_units( size_t number_of_P_units )
{
  _initial_number_of_P_units = number_of_P_units;
}

/**
 * \brief    Set point mutation rates
 * \details  --
 * \param    double rate
 * \return   \e void
 */
inline void Parameters::set_point_mutation_rate( double rate )
{
  assert(rate >= 0.0);
  assert(rate <= 1.0);
  _point_mutation_rate = rate;
}

/**
 * \brief    Set duplication rate
 * \details  --
 * \param    double rate
 * \return   \e void
 */
inline void Parameters::set_duplication_rate( double rate )
{
  assert(rate >= 0.0);
  assert(rate <= 1.0);
  _duplication_rate = rate;
}

/**
 * \brief    Set deletion rate
 * \details  --
 * \param    double rate
 * \return   \e void
 */
inline void Parameters::set_deletion_rate( double rate )
{
  assert(rate >= 0.0);
  assert(rate <= 1.0);
  _deletion_rate = rate;
}

/**
 * \brief    Set translocation rate
 * \details  --
 * \param    double rate
 * \return   \e void
 */
inline void Parameters::set_translocation_rate( double rate )
{
  assert(rate >= 0.0);
  assert(rate <= 1.0);
  _translocation_rate = rate;
}

/**
 * \brief    Set inversion rate
 * \details  --
 * \param    double rate
 * \return   \e void
 */
inline void Parameters::set_inversion_rate( double rate )
{
  assert(rate >= 0.0);
  assert(rate <= 1.0);
  _inversion_rate = rate;
}

/**
 * \brief    Set transition rate
 * \details  --
 * \param    double rate
 * \return   \e void
 */
inline void Parameters::set_transition_rate( double rate )
{
  assert(rate >= 0.0);
  assert(rate <= 1.0);
  _transition_rate = rate;
}


/**
 * \brief    Set breakpoint rate
 * \details  --
 * \param    double rate
 * \return   \e void
 */
inline void Parameters::set_breakpoint_rate( double rate )
{
  assert(rate >= 0.0);
  assert(rate <= 1.0);
  _breakpoint_rate = rate;
}

/**
 * \brief    Set substrate tag mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_substrate_tag_mutation_size( double size )
{
  assert(size >= 0.0);
  _substrate_tag_mutation_size = size;
}

/**
 * \brief    Set product tag mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_product_tag_mutation_size( double size )
{
  assert(size >= 0.0);
  _product_tag_mutation_size = size;
}

/**
 * \brief    Set Kcat constant mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_kcat_mutation_size( double size )
{
  assert(size >= 0.0);
  _kcat_mutation_size = size;
}

/**
 * \brief    Set Kcat/Km ratio constant mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_kcat_km_ratio_mutation_size( double size )
{
  assert(size >= 0.0);
  _kcat_km_ratio_mutation_size = size;
}

/**
 * \brief    Set binding site tag mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_binding_site_tag_mutation_size( double size )
{
  assert(size >= 0.0);
  _binding_site_tag_mutation_size = size;
}

/**
 * \brief    Set co-enzyme tag mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_co_enzyme_tag_mutation_size( double size )
{
  assert(size >= 0.0);
  _co_enzyme_tag_mutation_size = size;
}

/**
 * \brief    Set transcription factor tag mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_transcription_factor_tag_mutation_size( double size )
{
  assert(size >= 0.0);
  _transcription_factor_tag_mutation_size = size;
}

/**
 * \brief    Set basal expression level mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_basal_expression_level_mutation_size( double size )
{
  assert(size >= 0.0);
  _basal_expression_level_mutation_size = size;
}

/**
 * \brief    Set mutation of mutation rates
 * \details  --
 * \param    double rate
 * \return   \e void
 */
inline void Parameters::set_mutation_of_mutation_rates( double rate )
{
  assert(rate <= 1.0);
  assert(rate >= 0.0);
  _mutation_of_mutation_rates = rate;
}

/*------------------------------------------------------------------ genetic regulation network */

/**
 * \brief    Set the genetic regulation network ODE timestep
 * \details  --
 * \param    double genetic_regulation_network_timestep
 * \return   \e void
 */
inline void Parameters::set_genetic_regulation_network_timestep( double genetic_regulation_network_timestep )
{
  assert(genetic_regulation_network_timestep > 0.0);
  _genetic_regulation_network_timestep = genetic_regulation_network_timestep;
}

/**
 * \brief    Set the parameter theta of the Hill function
 * \details  --
 * \param    double hill_function_theta
 * \return   \e void
 */
inline void Parameters::set_hill_function_theta( double hill_function_theta )
{
  assert(hill_function_theta >= 0.0);
  //assert(hill_function_theta <= 1.0);
  _hill_function_theta = hill_function_theta;
}

/**
 * \brief    Set the parameter n of the Hill function
 * \details  --
 * \param    double hill_function_n
 * \return   \e void
 */
inline void Parameters::set_hill_function_n( double hill_function_n )
{
  assert(hill_function_n >= 0.0);
  _hill_function_n = hill_function_n;
}

/**
 * \brief    Set the protein degradation rate
 * \details  --
 * \param    double protein_degradation_rate
 * \return   \e void
 */
inline void Parameters::set_protein_degradation_rate( double protein_degradation_rate )
{
  assert(protein_degradation_rate >= 0.0);
  assert(protein_degradation_rate <= 1.0);
  _protein_degradation_rate = protein_degradation_rate;
}

/*------------------------------------------------------------------ metabolic network */

/**
 * \brief    Set metabolism ODE timestep
 * \details  --
 * \param    double metabolism_timestep
 * \return   \e void
 */
inline void Parameters::set_metabolism_timestep( double metabolism_timestep )
{
  assert(metabolism_timestep > 0.0);
  _metabolism_timestep = metabolism_timestep;
}

/**
 * \brief    Set essential metabolites toxicity threshold
 * \details  --
 * \param    double toxicity_threshold
 * \return   \e void
 */
inline void Parameters::set_essential_metabolites_toxicity_threshold( double toxicity_threshold )
{
  assert(toxicity_threshold > 0.0);
  _essential_metabolites_toxicity_threshold = toxicity_threshold;
}

/**
 * \brief    Set non essential metabolites toxicity threshold
 * \details  --
 * \param    double toxic_threshold
 * \return   \e void
 */
inline void Parameters::set_non_essential_metabolites_toxicity_threshold( double toxicity_threshold )
{
  assert(toxicity_threshold > 0.0);
  _non_essential_metabolites_toxicity_threshold = toxicity_threshold;
}

/**
 * \brief    Set initial metabolites amount in cells
 * \details  --
 * \param    double amount
 * \return   \e void
 */
inline void Parameters::set_initial_metabolites_amount_in_cells( double amount )
{
  assert(amount >= 0.0);
  _initial_metabolites_amount_in_cells = amount;
}

/**
 * \brief    Set the maximum reaction size
 * \details  --
 * \param    size_t size
 * \return   \e void
 */
inline void Parameters::set_maximum_reaction_size( size_t size )
{
  _maximum_reaction_size = size;
}

/*------------------------------------------------------------------ energy */

/**
 * \brief    Set energetical cost of transcription
 * \details  --
 * \param    double cost
 * \return   \e void
 */
inline void Parameters::set_energy_transcription_cost( double cost )
{
  assert(cost >= 0.0);
  _energy_transcription_cost = cost;
}

/**
 * \brief    Set energetical cost of protein degradation
 * \details  --
 * \param    double cost
 * \return   \e void
 */
inline void Parameters::set_energy_degradation_cost( double cost )
{
  assert(cost >= 0.0);
  _energy_degradation_cost = cost;
}

/**
 * \brief    Set energetical cost of enzymatic activity
 * \details  --
 * \param    double cost
 * \return   \e void
 */
inline void Parameters::set_energy_enzymatic_cost( double cost )
{
  assert(cost >= 0.0);
  _energy_enzymatic_cost = cost;
}

/**
 * \brief    Set energetical cost of pumping activity
 * \details  --
 * \param    double cost
 * \return   \e void
 */
inline void Parameters::set_energy_pumping_cost( double cost )
{
  assert(cost >= 0.0);
  _energy_pumping_cost = cost;
}

/**
 * \brief    Set energy dissipation rate
 * \details  --
 * \param    double rate
 * \return   \e void
 */
inline void Parameters::set_energy_dissipation_rate( double rate )
{
  assert(rate >= 0.0);
  assert(rate <= 1.0);
  _energy_dissipation_rate = rate;
}

/**
 * \brief    Set energy toxicity threshold
 * \details  --
 * \param    double threshold
 * \return   \e void
 */
inline void Parameters::set_energy_toxicity_threshold( double threshold )
{
  assert(threshold >= 0.0);
  _energy_toxicity_threshold = threshold;
}

/**
 * \brief    Set initial energy amount in cells
 * \details  --
 * \param    double amount
 * \return   \e void
 */
inline void Parameters::set_initial_energy_amount_in_cells( double amount )
{
  assert(amount >= 0.0);
  _initial_energy_amount_in_cells = amount;
}

/*------------------------------------------------------------------ cell */

/**
 * \brief    Set membrane permeability
 * \details  Defines the permeability of the cell membrane
 * \param    double permeability
 * \return   \e void
 */
inline void Parameters::set_membrane_permeability( double permeability )
{
  assert(permeability >= 0.0);
  assert(permeability <= 1.0);
  _membrane_permeability = permeability;
}

/*------------------------------------------------------------------ population */

/**
 * \brief    Set death probability
 * \details  --
 * \param    double death_probability
 * \return   \e void
 */
inline void Parameters::set_death_probability( double death_probability )
{
  assert(death_probability >= 0.0);
  assert(death_probability <= 1.0);
  _death_probability = death_probability;
}

/**
 * \brief    Set the migration rate
 * \details  --
 * \param    double migration_rate
 * \return   \e void
 */
inline void Parameters::set_migration_rate( double migration_rate )
{
  assert(migration_rate >= 0.0);
  assert(migration_rate <= 1.0);
  _migration_rate = migration_rate;
}

/**
 * \brief    Set the hgt rate
 * \details  --
 * \param    double hgt_rate
 * \return   \e void
 */
inline void Parameters::set_hgt_rate( double hgt_rate )
{
  assert(hgt_rate >= 0.0);
  assert(hgt_rate <= 1.0);
  _hgt_rate = hgt_rate;
}

/*------------------------------------------------------------------ environment */

/**
 * \brief    Set the environment properties
 * \details  --
 * \param    const environment_properties* properties
 * \return   \e void
 */
inline void Parameters::set_environment_properties( const environment_properties* properties )
{
  _environment_properties.number_of_init_cycles = properties->number_of_init_cycles;
  
  _environment_properties.species_tag_range.min = properties->species_tag_range.min;
  _environment_properties.species_tag_range.max = properties->species_tag_range.max;
  
  _environment_properties.concentration_range.min = properties->concentration_range.min;
  _environment_properties.concentration_range.max = properties->concentration_range.max;
  
  _environment_properties.number_of_species_range.min = properties->number_of_species_range.min;
  _environment_properties.number_of_species_range.max = properties->number_of_species_range.max;
  
  _environment_properties.interaction_scheme  = properties->interaction_scheme;
  _environment_properties.renewal_scheme      = properties->renewal_scheme;
  _environment_properties.variation_scheme    = properties->variation_scheme;
  _environment_properties.localization_scheme = properties->localization_scheme;
  _environment_properties.metabolic_scheme    = properties->metabolic_scheme;
  
  _environment_properties.introduction_rate     = properties->introduction_rate;
  _environment_properties.diffusion_coefficient = properties->diffusion_coefficient;
  _environment_properties.degradation_rate      = properties->degradation_rate;
}

/*------------------------------------------------------------------ prngs */

/**
 * \brief    Set the simulation prng
 * \details  --
 * \param    Prng* prng
 * \return   \e void
 */
inline void Parameters::set_simulation_prng( Prng* prng )
{
  delete _simulation_prng;
  _simulation_prng = NULL;
  _simulation_prng = new Prng(*prng);
}

/**
 * \brief    Set the population prng at position pos
 * \details  --
 * \param    Prng* prng
 * \param    size_t pos
 * \return   \e void
 */
inline void Parameters::set_population_prng( Prng* prng, size_t pos )
{
  assert(pos < _width*_height);
  delete _population_prng[pos];
  _population_prng[pos] = NULL;
  _population_prng[pos] = new Prng(*prng);
}

/**
 * \brief    Set the environment prng
 * \details  --
 * \param    Prng* prng
 * \return   \e void
 */
inline void Parameters::set_environment_prng( Prng* prng )
{
  delete _environment_prng;
  _environment_prng = NULL;
  _environment_prng = new Prng(*prng);
}


#endif /* defined(__EVOEVO__Parameters__) */
