
/**
 * \file      Parameters.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Parameters class definition
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

#include "Parameters.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Default constructor
 * \details  --
 * \param    void
 * \return   \e void
 */
Parameters::Parameters( void )
{
  /*------------------------------------------------------------------ prng seed */
  
  _seed = 0;
  
  /*------------------------------------------------------------------ parallel computing */
  
  _parallel_computing = false;
  
  /*------------------------------------------------------------------ simulation schemes */
  
  _energy_costs_scheme          = false;
  _membrane_permeability_scheme = false;
  _metabolic_inheritance_scheme = false;
  _enzymatic_inheritance_scheme = false;
  _co_enzyme_activity_scheme    = false;
  _score_scheme                 = ESSENTIAL_METABOLITES_SUM;
  _selection_threshold          = 0.0;
  
  /*------------------------------------------------------------------ space */
  
  _width  = 0;
  _height = 0;
  
  /*------------------------------------------------------------------ output */
  
  _simulation_backup_step  = 0;
  _figures_generation_step = 0;
  
  /*------------------------------------------------------------------ genome */
  
  _load_genome_from_file = false;
  
  _metabolite_tag_initial_range.min = 0;
  _metabolite_tag_initial_range.max = 0;
  
  _binding_site_tag_initial_range.min = 0;
  _binding_site_tag_initial_range.max = 0;
  
  _co_enzyme_tag_initial_range.min = 0;
  _co_enzyme_tag_initial_range.max = 0;
  
  _transcription_factor_tag_initial_range.min = 0;
  _transcription_factor_tag_initial_range.max = 0;
  
  _transcription_factor_binding_window = 0;
  
  _initial_number_of_NC_units = 0;
  _initial_number_of_E_units  = 0;
  _initial_number_of_TF_units = 0;
  _initial_number_of_BS_units = 0;
  _initial_number_of_P_units  = 0;
  
  _point_mutation_rate = 0.0;
  _duplication_rate    = 0.0;
  _deletion_rate       = 0.0;
  _translocation_rate  = 0.0;
  _inversion_rate      = 0.0;
  _transition_rate     = 0.0;
  _breakpoint_rate     = 0.0;
  
  _substrate_tag_mutation_size            = 0;
  _product_tag_mutation_size              = 0;
  _kcat_mutation_size                     = 0.0;
  _kcat_km_ratio_mutation_size            = 0.0;
  _binding_site_tag_mutation_size         = 0;
  _co_enzyme_tag_mutation_size            = 0;
  _transcription_factor_tag_mutation_size = 0;
  _basal_expression_level_mutation_size   = 0.0;
  
  _mutation_of_mutation_rates = 0.0;
  
  /*------------------------------------------------------------------ genetic regulation network */
  
  _genetic_regulation_network_timestep = 0.0;
  _hill_function_theta                 = 0.0;
  _hill_function_n                     = 0.0;
  _protein_degradation_rate            = 0.0;
  
  /*------------------------------------------------------------------ metabolic network */
  
  _metabolism_timestep                          = 0.0;
  _essential_metabolites_toxicity_threshold     = 0.0;
  _non_essential_metabolites_toxicity_threshold = 0.0;
  _initial_metabolites_amount_in_cells          = 0.0;
  _maximum_reaction_size                        = 0;
  
  /*------------------------------------------------------------------ energy */
  
  _energy_transcription_cost = 0.0;
  _energy_degradation_cost   = 0.0;
  _energy_enzymatic_cost     = 0.0;
  _energy_pumping_cost       = 0.0;
  
  _energy_dissipation_rate = 0.0;
  
  _energy_toxicity_threshold = 0.0;
  
  _initial_energy_amount_in_cells = 0.0;
  
  /*------------------------------------------------------------------ cell */
  
  _membrane_permeability = 0.0;
  
  /*------------------------------------------------------------------ population */
  
  _death_probability = 0.0;
  _migration_rate    = 0.0;
  _hgt_rate          = 0.0;
  
  /*------------------------------------------------------------------ environment */
  
  _environment_properties.number_of_init_cycles = 0;
  
  _environment_properties.species_tag_range.min = 0;
  _environment_properties.species_tag_range.max = 0;
  
  _environment_properties.concentration_range.min = 0;
  _environment_properties.concentration_range.max = 0;
  
  _environment_properties.number_of_species_range.min = 0;
  _environment_properties.number_of_species_range.max = 0;
  
  _environment_properties.interaction_scheme  = INTERACTION;
  _environment_properties.renewal_scheme      = KEEP_MATTER;
  _environment_properties.variation_scheme    = PERIODIC_SCHEME;
  _environment_properties.localization_scheme = GLOBAL_LOCALIZATION;
  _environment_properties.metabolic_scheme    = MULTIPLE_METABOLITES;
  
  _environment_properties.introduction_rate     = 0.0;
  _environment_properties.diffusion_coefficient = 0.0;
  _environment_properties.degradation_rate      = 0.0;
  
  /*------------------------------------------------------------------ prime numbers */
  
  build_prime_numbers_list(PRIME_NUMBERS_LIST_LENGTH);
  
  /*------------------------------------------------------------------ prngs */
  
  _simulation_prng  = NULL;
  _population_prng  = NULL;
  _environment_prng = NULL;
  
}

/**
 * \brief    Contructor from backup file
 * \details  Load Parameter class from backup file
 * \param    size_t backup_time
 * \return   \e void
 */
Parameters::Parameters( size_t backup_time )
{
  std::stringstream backup_file_name;
  backup_file_name << "./parameters/parameters_" << backup_time;
  gzFile backup_file = gzopen(backup_file_name.str().c_str(), "r");
  
  /*------------------------------------------------------------------ prng seed */
  
  gzread( backup_file, &_seed, sizeof(_seed) );
  
  /*------------------------------------------------------------------ parallel computing */
  
  gzread( backup_file, &_parallel_computing, sizeof(_parallel_computing) );
  
  /*------------------------------------------------------------------ simulation schemes */
  
  gzread( backup_file, &_energy_costs_scheme,          sizeof(_energy_costs_scheme) );
  gzread( backup_file, &_membrane_permeability_scheme, sizeof(_membrane_permeability_scheme) );
  gzread( backup_file, &_metabolic_inheritance_scheme, sizeof(_metabolic_inheritance_scheme) );
  gzread( backup_file, &_enzymatic_inheritance_scheme, sizeof(_enzymatic_inheritance_scheme) );
  gzread( backup_file, &_co_enzyme_activity_scheme,    sizeof(_co_enzyme_activity_scheme) );
  gzread( backup_file, &_score_scheme,                 sizeof(_score_scheme) );
  gzread( backup_file, &_selection_threshold,          sizeof(_selection_threshold) );
  
  /*------------------------------------------------------------------ space */
  
  gzread( backup_file, &_width,  sizeof(_width) );
  gzread( backup_file, &_height, sizeof(_height) );
  
  /*------------------------------------------------------------------ output */
  
  gzread( backup_file, &_simulation_backup_step,  sizeof(_simulation_backup_step) );
  gzread( backup_file, &_figures_generation_step, sizeof(_figures_generation_step) );
  
  /*------------------------------------------------------------------ genome */
  
  gzread( backup_file, &_load_genome_from_file, sizeof(_load_genome_from_file) );
  
  gzread( backup_file, &_metabolite_tag_initial_range.min, sizeof(_metabolite_tag_initial_range.min) );
  gzread( backup_file, &_metabolite_tag_initial_range.max, sizeof(_metabolite_tag_initial_range.max) );
  
  gzread( backup_file, &_binding_site_tag_initial_range.min, sizeof(_binding_site_tag_initial_range.min) );
  gzread( backup_file, &_binding_site_tag_initial_range.max, sizeof(_binding_site_tag_initial_range.max) );
  
  gzread( backup_file, &_co_enzyme_tag_initial_range.min, sizeof(_co_enzyme_tag_initial_range.min) );
  gzread( backup_file, &_co_enzyme_tag_initial_range.max, sizeof(_co_enzyme_tag_initial_range.max) );
  
  gzread( backup_file, &_transcription_factor_tag_initial_range.min, sizeof(_transcription_factor_tag_initial_range.min) );
  gzread( backup_file, &_transcription_factor_tag_initial_range.max, sizeof(_transcription_factor_tag_initial_range.max) );
  
  gzread( backup_file, &_transcription_factor_binding_window, sizeof(_transcription_factor_binding_window) );
  
  gzread( backup_file, &_initial_number_of_NC_units, sizeof(_initial_number_of_NC_units) );
  gzread( backup_file, &_initial_number_of_E_units,  sizeof(_initial_number_of_E_units) );
  gzread( backup_file, &_initial_number_of_TF_units, sizeof(_initial_number_of_TF_units) );
  gzread( backup_file, &_initial_number_of_BS_units, sizeof(_initial_number_of_BS_units) );
  gzread( backup_file, &_initial_number_of_P_units,  sizeof(_initial_number_of_P_units) );
  
  gzread( backup_file, &_point_mutation_rate, sizeof(_point_mutation_rate) );
  gzread( backup_file, &_duplication_rate,    sizeof(_duplication_rate) );
  gzread( backup_file, &_deletion_rate,       sizeof(_deletion_rate) );
  gzread( backup_file, &_translocation_rate,  sizeof(_translocation_rate) );
  gzread( backup_file, &_inversion_rate,      sizeof(_inversion_rate) );
  gzread( backup_file, &_transition_rate,     sizeof(_transition_rate) );
  gzread( backup_file, &_breakpoint_rate,     sizeof(_breakpoint_rate) );
  
  gzread( backup_file, &_substrate_tag_mutation_size,            sizeof(_substrate_tag_mutation_size) );
  gzread( backup_file, &_product_tag_mutation_size,              sizeof(_product_tag_mutation_size) );
  gzread( backup_file, &_kcat_mutation_size,                     sizeof(_kcat_mutation_size) );
  gzread( backup_file, &_kcat_km_ratio_mutation_size,            sizeof(_kcat_km_ratio_mutation_size) );
  gzread( backup_file, &_binding_site_tag_mutation_size,         sizeof(_binding_site_tag_mutation_size) );
  gzread( backup_file, &_co_enzyme_tag_mutation_size,            sizeof(_co_enzyme_tag_mutation_size) );
  gzread( backup_file, &_transcription_factor_tag_mutation_size, sizeof(_transcription_factor_tag_mutation_size) );
  gzread( backup_file, &_basal_expression_level_mutation_size,   sizeof(_basal_expression_level_mutation_size) );
  
  gzread( backup_file, &_mutation_of_mutation_rates, sizeof(_mutation_of_mutation_rates) );
  
  /*------------------------------------------------------------------ genetic regulation network */
  
  gzread( backup_file, &_genetic_regulation_network_timestep, sizeof(_genetic_regulation_network_timestep) );
  gzread( backup_file, &_hill_function_theta,                 sizeof(_hill_function_theta) );
  gzread( backup_file, &_hill_function_n,                     sizeof(_hill_function_n) );
  gzread( backup_file, &_protein_degradation_rate,            sizeof(_protein_degradation_rate) );
  
  /*------------------------------------------------------------------ metabolic network */
  
  gzread( backup_file, &_metabolism_timestep,                          sizeof(_metabolism_timestep) );
  gzread( backup_file, &_essential_metabolites_toxicity_threshold,     sizeof(_essential_metabolites_toxicity_threshold) );
  gzread( backup_file, &_non_essential_metabolites_toxicity_threshold, sizeof(_non_essential_metabolites_toxicity_threshold) );
  gzread( backup_file, &_initial_metabolites_amount_in_cells,          sizeof(_initial_metabolites_amount_in_cells) );
  
  /*------------------------------------------------------------------ energy */
  
  gzread( backup_file, &_energy_transcription_cost, sizeof(_energy_transcription_cost) );
  gzread( backup_file, &_energy_degradation_cost,   sizeof(_energy_degradation_cost) );
  gzread( backup_file, &_energy_enzymatic_cost,     sizeof(_energy_enzymatic_cost) );
  gzread( backup_file, &_energy_pumping_cost,       sizeof(_energy_pumping_cost) );
  
  gzread( backup_file, &_energy_dissipation_rate, sizeof(_energy_dissipation_rate) );
  
  gzread( backup_file, &_energy_toxicity_threshold, sizeof(_energy_toxicity_threshold) );
  
  gzread( backup_file, &_initial_energy_amount_in_cells, sizeof(_initial_energy_amount_in_cells) );
  
  gzread( backup_file, &_maximum_reaction_size, sizeof(_maximum_reaction_size) );
  
  /*------------------------------------------------------------------ cell */
  
  gzread( backup_file, &_membrane_permeability, sizeof(_membrane_permeability) );
  
  /*------------------------------------------------------------------ population */
  
  gzread( backup_file, &_death_probability, sizeof(_death_probability) );
  gzread( backup_file, &_migration_rate,    sizeof(_migration_rate) );
  gzread( backup_file, &_hgt_rate,          sizeof(_hgt_rate) );
  
  /*------------------------------------------------------------------ environment */
  
  gzread( backup_file, &_environment_properties.number_of_init_cycles, sizeof(_environment_properties.number_of_init_cycles) );
  
  gzread( backup_file, &_environment_properties.species_tag_range.min, sizeof(_environment_properties.species_tag_range.min) );
  gzread( backup_file, &_environment_properties.species_tag_range.max, sizeof(_environment_properties.species_tag_range.max) );
  
  gzread( backup_file, &_environment_properties.concentration_range.min, sizeof(_environment_properties.concentration_range.min) );
  gzread( backup_file, &_environment_properties.concentration_range.max, sizeof(_environment_properties.concentration_range.max) );
  
  gzread( backup_file, &_environment_properties.number_of_species_range.min, sizeof(_environment_properties.number_of_species_range.min) );
  gzread( backup_file, &_environment_properties.number_of_species_range.max, sizeof(_environment_properties.number_of_species_range.max) );
  
  gzread( backup_file, &_environment_properties.interaction_scheme,  sizeof(_environment_properties.interaction_scheme) );
  gzread( backup_file, &_environment_properties.renewal_scheme,      sizeof(_environment_properties.renewal_scheme) );
  gzread( backup_file, &_environment_properties.variation_scheme,    sizeof(_environment_properties.variation_scheme) );
  gzread( backup_file, &_environment_properties.localization_scheme, sizeof(_environment_properties.localization_scheme) );
  gzread( backup_file, &_environment_properties.metabolic_scheme,    sizeof(_environment_properties.metabolic_scheme) );
  
  gzread( backup_file, &_environment_properties.introduction_rate,     sizeof(_environment_properties.introduction_rate) );
  gzread( backup_file, &_environment_properties.diffusion_coefficient, sizeof(_environment_properties.diffusion_coefficient) );
  gzread( backup_file, &_environment_properties.degradation_rate,      sizeof(_environment_properties.degradation_rate) );
  
  gzclose(backup_file);
  
  /*------------------------------------------------------------------ prime numbers */
  
  build_prime_numbers_list(PRIME_NUMBERS_LIST_LENGTH);
  
  /*------------------------------------------------------------------ prngs */
  
  std::stringstream prng_file_name;
  prng_file_name << "./prng/prng_" << backup_time;
  FILE * prng_file = fopen(prng_file_name.str().c_str(), "r");
  _simulation_prng = new Prng(prng_file);
  _population_prng = new Prng*[_width*_height];
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    _population_prng[pos] = new Prng(prng_file);
  }
  _environment_prng = new Prng(prng_file);
  fclose(prng_file);
}

/**
 * \brief    Copy constructor
 * \details  --
 * \param    const Parameters& parameters
 * \return   \e void
 */
Parameters::Parameters( const Parameters& parameters )
{
  /*------------------------------------------------------------------ prng seed */
  
  _seed = parameters._seed;
  
  /*------------------------------------------------------------------ parallel computing */
  
  _parallel_computing = parameters._parallel_computing;
  
  /*------------------------------------------------------------------ simulation schemes */
  
  _energy_costs_scheme          = parameters._energy_costs_scheme;
  _membrane_permeability_scheme = parameters._membrane_permeability_scheme;
  _metabolic_inheritance_scheme = parameters._metabolic_inheritance_scheme;
  _enzymatic_inheritance_scheme = parameters._enzymatic_inheritance_scheme;
  _co_enzyme_activity_scheme    = parameters._co_enzyme_activity_scheme;
  _score_scheme                 = parameters._score_scheme;
  _selection_threshold          = parameters._selection_threshold;
  
  /*------------------------------------------------------------------ space */
  
  _width  = parameters._width;
  _height = parameters._height;
  
  /*------------------------------------------------------------------ output */
  
  _simulation_backup_step  = parameters._simulation_backup_step;
  _figures_generation_step = parameters._figures_generation_step;
  
  /*------------------------------------------------------------------ genome */
  
  _load_genome_from_file = parameters._load_genome_from_file;
  
  _metabolite_tag_initial_range.min = parameters._metabolite_tag_initial_range.min;
  _metabolite_tag_initial_range.max = parameters._metabolite_tag_initial_range.max;
  
  _binding_site_tag_initial_range.min = parameters._binding_site_tag_initial_range.min;
  _binding_site_tag_initial_range.max = parameters._binding_site_tag_initial_range.max;
  
  _co_enzyme_tag_initial_range.min = parameters._co_enzyme_tag_initial_range.min;
  _co_enzyme_tag_initial_range.max = parameters._co_enzyme_tag_initial_range.max;
  
  _transcription_factor_tag_initial_range.min = parameters._transcription_factor_tag_initial_range.min;
  _transcription_factor_tag_initial_range.max = parameters._transcription_factor_tag_initial_range.max;
  
  _transcription_factor_binding_window = parameters._transcription_factor_binding_window;
  
  _initial_number_of_NC_units = parameters._initial_number_of_NC_units;
  _initial_number_of_E_units  = parameters._initial_number_of_E_units;
  _initial_number_of_TF_units = parameters._initial_number_of_TF_units;
  _initial_number_of_BS_units = parameters._initial_number_of_BS_units;
  _initial_number_of_P_units  = parameters._initial_number_of_P_units;
  
  _point_mutation_rate = parameters._point_mutation_rate;
  _duplication_rate    = parameters._duplication_rate;
  _deletion_rate       = parameters._deletion_rate;
  _translocation_rate  = parameters._translocation_rate;
  _inversion_rate      = parameters._inversion_rate;
  _transition_rate     = parameters._transition_rate;
  _breakpoint_rate     = parameters._breakpoint_rate;
  
  _substrate_tag_mutation_size            = parameters._substrate_tag_mutation_size;
  _product_tag_mutation_size              = parameters._product_tag_mutation_size;
  _kcat_mutation_size                     = parameters._kcat_mutation_size;
  _kcat_km_ratio_mutation_size            = parameters._kcat_km_ratio_mutation_size;
  _binding_site_tag_mutation_size         = parameters._binding_site_tag_mutation_size;
  _co_enzyme_tag_mutation_size            = parameters._co_enzyme_tag_mutation_size;
  _transcription_factor_tag_mutation_size = parameters._transcription_factor_tag_mutation_size;
  _basal_expression_level_mutation_size   = parameters._basal_expression_level_mutation_size;
  
  _mutation_of_mutation_rates = parameters._mutation_of_mutation_rates;
  
  /*------------------------------------------------------------------ genetic regulation network */
  
  _genetic_regulation_network_timestep = parameters._genetic_regulation_network_timestep;
  _hill_function_theta                 = parameters._hill_function_theta;
  _hill_function_n                     = parameters._hill_function_n;
  _protein_degradation_rate            = parameters._protein_degradation_rate;
  
  /*------------------------------------------------------------------ metabolic network */
  
  _metabolism_timestep                          = parameters._metabolism_timestep;
  _essential_metabolites_toxicity_threshold     = parameters._essential_metabolites_toxicity_threshold;
  _non_essential_metabolites_toxicity_threshold = parameters._non_essential_metabolites_toxicity_threshold;
  _initial_metabolites_amount_in_cells          = parameters._initial_metabolites_amount_in_cells;
  _maximum_reaction_size                        = parameters._maximum_reaction_size;
  
  /*------------------------------------------------------------------ energy */
  
  _energy_transcription_cost = parameters._energy_transcription_cost;
  _energy_degradation_cost   = parameters._energy_degradation_cost;
  _energy_enzymatic_cost     = parameters._energy_enzymatic_cost;
  _energy_pumping_cost       = parameters._energy_pumping_cost;
  
  _energy_dissipation_rate = parameters._energy_dissipation_rate;
  
  _energy_toxicity_threshold = parameters._energy_toxicity_threshold;
  
  _initial_energy_amount_in_cells = parameters._initial_energy_amount_in_cells;
  
  /*------------------------------------------------------------------ cell */
  
  _membrane_permeability = parameters._membrane_permeability;
  
  /*------------------------------------------------------------------ population */
  
  _death_probability = parameters._death_probability;
  _migration_rate    = parameters._migration_rate;
  _hgt_rate          = parameters._hgt_rate;
  
  /*------------------------------------------------------------------ environment */
  
  _environment_properties.number_of_init_cycles = parameters._environment_properties.number_of_init_cycles;
  
  _environment_properties.species_tag_range.min = parameters._environment_properties.species_tag_range.min;
  _environment_properties.species_tag_range.max = parameters._environment_properties.species_tag_range.max;
  
  _environment_properties.concentration_range.min = parameters._environment_properties.concentration_range.min;
  _environment_properties.concentration_range.max = parameters._environment_properties.concentration_range.max;
  
  _environment_properties.number_of_species_range.min = parameters._environment_properties.number_of_species_range.min;
  _environment_properties.number_of_species_range.max = parameters._environment_properties.number_of_species_range.max;
  
  _environment_properties.interaction_scheme  = parameters._environment_properties.interaction_scheme;
  _environment_properties.renewal_scheme      = parameters._environment_properties.renewal_scheme;
  _environment_properties.variation_scheme    = parameters._environment_properties.variation_scheme;
  _environment_properties.localization_scheme = parameters._environment_properties.localization_scheme;
  _environment_properties.metabolic_scheme    = parameters._environment_properties.metabolic_scheme;
  
  _environment_properties.introduction_rate     = parameters._environment_properties.introduction_rate;
  _environment_properties.diffusion_coefficient = parameters._environment_properties.diffusion_coefficient;
  _environment_properties.degradation_rate      = parameters._environment_properties.degradation_rate;
  
  /*------------------------------------------------------------------ prime numbers */
  
  build_prime_numbers_list(PRIME_NUMBERS_LIST_LENGTH);
  
  /*------------------------------------------------------------------ prngs */
  
  _simulation_prng = new Prng(*parameters._simulation_prng);
  _population_prng = new Prng*[_width*_height];
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    _population_prng[pos] = new Prng(*parameters._population_prng[pos]);
  }
  _environment_prng = new Prng(*parameters._environment_prng);
  
}

/*----------------------------
 * DESTRUCTORS
 *----------------------------*/

/**
 * \brief    Destructor
 * \details  --
 * \param    void
 * \return   \e void
 */
Parameters::~Parameters( void )
{
  delete[] _prime_numbers;
  _prime_numbers = NULL;
  delete _simulation_prng;
  _simulation_prng = NULL;
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    delete _population_prng[pos];
    _population_prng[pos] = NULL;
  }
  delete[] _population_prng;
  _population_prng = NULL;
  delete _environment_prng;
  _environment_prng = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Load parameters from a parameters file
 * \details  Parameters are loaded from a text file (default name: parameters.txt)
 * \param    std::string filename
 * \return   \e void
 */
void Parameters::load_parameters_from_file( std::string filename )
{
  std::ifstream file(filename.c_str(), std::ios::in);
  if(!file)
  {
    std::cout << "Error: " << filename << " file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::string line;
    std::string param_name;
    while(getline(file, line))
    {
      std::vector<std::string> words;
      if(parse_line(&words, line))
      {
        /*------------------------------------------------------------------ prng seed */
        
        if ( strcmp(words[0].c_str(), "SEED") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _seed;
          assert(_seed > 0);
          assert(_seed <= MAXIMUM_SEED);
        }
        
        /*------------------------------------------------------------------ parallel computing */
        
        else if ( strcmp(words[0].c_str(), "PARALLEL_COMPUTING") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> state;
          if (strcmp(state.c_str(), "NO") == 0 || strcmp(state.c_str(), "no") == 0)
          {
            _parallel_computing = false;
          }
          else if (strcmp(state.c_str(), "YES") == 0 || strcmp(state.c_str(), "yes") == 0)
          {
            _parallel_computing = true;
          }
          else
          {
            std::cout << "Error : PARALLEL_COMPUTING wrong value at line:\n " << line.c_str() << ".\n\n";
            exit(EXIT_FAILURE);
          }
        }
        
        /*------------------------------------------------------------------ simulation schemes */
        
        else if ( strcmp(words[0].c_str(), "ENERGY_COSTS_SCHEME") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> state;
          if (strcmp(state.c_str(), "NO") == 0 || strcmp(state.c_str(), "no") == 0)
          {
            _energy_costs_scheme = false;
          }
          else if (strcmp(state.c_str(), "YES") == 0 || strcmp(state.c_str(), "yes") == 0)
          {
            _energy_costs_scheme = true;
          }
          else
          {
            std::cout << "Error : ENERGY_COSTS_SCHEME wrong value at line:\n " << line.c_str() << ".\n\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "MEMBRANE_PERMEABILITY_SCHEME") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> state;
          if (strcmp(state.c_str(), "NO") == 0 || strcmp(state.c_str(), "no") == 0)
          {
            _membrane_permeability_scheme = false;
          }
          else if (strcmp(state.c_str(), "YES") == 0 || strcmp(state.c_str(), "yes") == 0)
          {
            _membrane_permeability_scheme = true;
          }
          else
          {
            std::cout << "Error : MEMBRANE_PERMEABILITY_SCHEME wrong value at line:\n " << line.c_str() << ".\n\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "METABOLIC_INHERITANCE_SCHEME") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> state;
          if (strcmp(state.c_str(), "NO") == 0 || strcmp(state.c_str(), "no") == 0)
          {
            _metabolic_inheritance_scheme = false;
          }
          else if (strcmp(state.c_str(), "YES") == 0 || strcmp(state.c_str(), "yes") == 0)
          {
            _metabolic_inheritance_scheme = true;
          }
          else
          {
            std::cout << "Error : METABOLIC_INHERITANCE_SCHEME wrong value at line:\n " << line.c_str() << ".\n\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "ENZYMATIC_INHERITANCE_SCHEME") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> state;
          if (strcmp(state.c_str(), "NO") == 0 || strcmp(state.c_str(), "no") == 0)
          {
            _enzymatic_inheritance_scheme = false;
          }
          else if (strcmp(state.c_str(), "YES") == 0 || strcmp(state.c_str(), "yes") == 0)
          {
            _enzymatic_inheritance_scheme = true;
          }
          else
          {
            std::cout << "Error : ENZYMATIC_INHERITANCE_SCHEME wrong value at line:\n " << line.c_str() << ".\n\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "CO_ENZYME_ACTIVITY_SCHEME") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> state;
          if (strcmp(state.c_str(), "NO") == 0 || strcmp(state.c_str(), "no") == 0)
          {
            _co_enzyme_activity_scheme = false;
          }
          else if (strcmp(state.c_str(), "YES") == 0 || strcmp(state.c_str(), "yes") == 0)
          {
            _co_enzyme_activity_scheme = true;
          }
          else
          {
            std::cout << "Error : CO_ENZYME_ACTIVITY_SCHEME wrong value at line:\n " << line.c_str() << ".\n\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "SCORE_SCHEME") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string law;
          flux >> param_name >> law;
          if (strcmp(law.c_str(), "SUM") == 0)
          {
            _score_scheme = ESSENTIAL_METABOLITES_SUM;
          }
          else if (strcmp(law.c_str(), "SUM_MINUS_DEV") == 0)
          {
            _score_scheme = ESSENTIAL_METABOLITES_SUM_MINUS_DEVIATION;
          }
          else if (strcmp(law.c_str(), "COMBINATORIAL") == 0)
          {
            _score_scheme = ESSENTIAL_METABOLITES_COMBINATORIAL_CONTRIBUTION;
          }
          else
          {
            std::cout << "Error : SCORE_SCHEME wrong value at line:\n " << line.c_str() << ".\n\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "SELECTION_THRESHOLD") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _selection_threshold;
          assert(_selection_threshold >= 0.0);
          assert(_selection_threshold <= 1.0);
        }
        
        /*------------------------------------------------------------------ space */
        
        else if ( strcmp(words[0].c_str(), "WIDTH") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _width;
          assert(_width > 0);
        }
        else if ( strcmp(words[0].c_str(), "HEIGHT") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _height;
          assert(_height > 0);
        }
        
        /*------------------------------------------------------------------ output */
        
        else if ( strcmp(words[0].c_str(), "SIMULATION_BACKUP_STEP") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _simulation_backup_step;
        }
        else if ( strcmp(words[0].c_str(), "FIGURES_GENERATION_STEP") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _figures_generation_step;
        }
        
        /*------------------------------------------------------------------ genome */
        
        else if ( strcmp(words[0].c_str(), "LOAD_GENOME_FROM_FILE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> state;
          if (strcmp(state.c_str(), "NO") == 0 || strcmp(state.c_str(), "no") == 0)
          {
            _load_genome_from_file = false;
          }
          else if (strcmp(state.c_str(), "YES") == 0 || strcmp(state.c_str(), "yes") == 0)
          {
            _load_genome_from_file = true;
          }
          else
          {
            std::cout << "Error : LOAD_GENOME_FROM_FILE wrong value at line:\n " << line.c_str() << ".\n\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "METABOLITE_TAG_INITIAL_RANGE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _metabolite_tag_initial_range.min >> _metabolite_tag_initial_range.max;
        }
        else if ( strcmp(words[0].c_str(), "BINDING_SITE_TAG_INITIAL_RANGE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _binding_site_tag_initial_range.min >> _binding_site_tag_initial_range.max;
        }
        else if ( strcmp(words[0].c_str(), "CO_ENZYME_TAG_INITIAL_RANGE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _co_enzyme_tag_initial_range.min >> _co_enzyme_tag_initial_range.max;
        }
        else if ( strcmp(words[0].c_str(), "TRANSCRIPTION_FACTOR_TAG_INITIAL_RANGE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _transcription_factor_tag_initial_range.min >> _transcription_factor_tag_initial_range.max;
        }
        else if ( strcmp(words[0].c_str(), "TRANSCRIPTION_FACTOR_BINDING_WINDOW") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _transcription_factor_binding_window;
        }
        else if ( strcmp(words[0].c_str(), "INITIAL_NUMBER_OF_NON_CODING_UNITS") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _initial_number_of_NC_units;
        }
        else if ( strcmp(words[0].c_str(), "INITIAL_NUMBER_OF_ENZYME_UNITS") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _initial_number_of_E_units;
        }
        else if ( strcmp(words[0].c_str(), "INITIAL_NUMBER_OF_TRANSCRIPTION_FACTOR_UNITS") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _initial_number_of_TF_units;
        }
        else if ( strcmp(words[0].c_str(), "INITIAL_NUMBER_OF_BINDING_SITE_UNITS") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _initial_number_of_BS_units;
        }
        else if ( strcmp(words[0].c_str(), "INITIAL_NUMBER_OF_PROMOTER_UNITS") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _initial_number_of_P_units;
        }
        else if ( strcmp(words[0].c_str(), "POINT_MUTATION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _point_mutation_rate;
          assert(_point_mutation_rate >= 0.0);
          assert(_point_mutation_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "DUPLICATION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _duplication_rate;
          assert(_duplication_rate >= 0.0);
          assert(_duplication_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "DELETION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _deletion_rate;
          assert(_deletion_rate >= 0.0);
          assert(_deletion_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "TRANSLOCATION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _translocation_rate;
          assert(_translocation_rate >= 0.0);
          assert(_translocation_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "INVERSION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _inversion_rate;
          assert(_inversion_rate >= 0.0);
          assert(_inversion_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "TRANSITION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _transition_rate;
          assert(_transition_rate >= 0.0);
          assert(_transition_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "BREAKPOINT_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _breakpoint_rate;
          assert(_breakpoint_rate >= 0.0);
          assert(_breakpoint_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "SUBSTRATE_TAG_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _substrate_tag_mutation_size;
          assert(_substrate_tag_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "PRODUCT_TAG_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _product_tag_mutation_size;
          assert(_product_tag_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "KCAT_KM_RATIO_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _kcat_km_ratio_mutation_size;
          assert(_kcat_km_ratio_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "KCAT_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _kcat_mutation_size;
          assert(_kcat_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "BINDING_SITE_TAG_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _binding_site_tag_mutation_size;
          assert(_binding_site_tag_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "CO_ENZYME_TAG_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _co_enzyme_tag_mutation_size;
          assert(_co_enzyme_tag_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _transcription_factor_tag_mutation_size;
          assert(_transcription_factor_tag_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "BASAL_EXPRESSION_LEVEL_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _basal_expression_level_mutation_size;
          assert(_basal_expression_level_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "MUTATION_OF_MUTATION_RATES") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _mutation_of_mutation_rates;
          assert(_mutation_of_mutation_rates >= 0.0);
          assert(_mutation_of_mutation_rates <= 1.0);
        }
        
        /*------------------------------------------------------------------ genetic regulation network */
        
        else if ( strcmp(words[0].c_str(), "GENETIC_REGULATION_NETWORK_TIMESTEP") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _genetic_regulation_network_timestep;
          assert(_genetic_regulation_network_timestep > 0.0);
        }
        else if ( strcmp(words[0].c_str(), "HILL_FUNCTION_THETA") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _hill_function_theta;
          assert(_hill_function_theta >= 0.0);
          //assert(_hill_function_theta <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "HILL_FUNCTION_N") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _hill_function_n;
          assert(_hill_function_n >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "PROTEIN_DEGRADATION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _protein_degradation_rate;
          assert(_protein_degradation_rate >= 0.0);
          assert(_protein_degradation_rate <= 1.0);
        }
        
        /*------------------------------------------------------------------ metabolic network */
        
        else if ( strcmp(words[0].c_str(), "METABOLISM_TIMESTEP") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _metabolism_timestep;
          assert(_metabolism_timestep > 0.0);
        }
        else if ( strcmp(words[0].c_str(), "ESSENTIAL_METABOLITES_TOXICITY_THRESHOLD") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _essential_metabolites_toxicity_threshold;
          assert(_essential_metabolites_toxicity_threshold > 0.0);
        }
        else if ( strcmp(words[0].c_str(), "NON_ESSENTIAL_METABOLITES_TOXICITY_THRESHOLD") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _non_essential_metabolites_toxicity_threshold;
          assert(_non_essential_metabolites_toxicity_threshold > 0.0);
        }
        else if ( strcmp(words[0].c_str(), "INITIAL_METABOLITES_AMOUNT_IN_CELLS") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _initial_metabolites_amount_in_cells;
          assert(_initial_metabolites_amount_in_cells >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "MAXIMUM_REACTION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _maximum_reaction_size;
        }
        
        /*------------------------------------------------------------------ energy */
        
        else if ( strcmp(words[0].c_str(), "ENERGY_TRANSCRIPTION_COST") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _energy_transcription_cost;
          assert(_energy_transcription_cost >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "ENERGY_DEGRADATION_COST") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _energy_degradation_cost;
          assert(_energy_degradation_cost >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "ENERGY_ENZYMATIC_COST") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _energy_enzymatic_cost;
          assert(_energy_enzymatic_cost >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "ENERGY_PUMPING_COST") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _energy_pumping_cost;
          assert(_energy_pumping_cost >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "ENERGY_DISSIPATION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _energy_dissipation_rate;
          assert(_energy_dissipation_rate >= 0.0);
          assert(_energy_dissipation_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "ENERGY_TOXICITY_THRESHOLD") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _energy_toxicity_threshold;
          assert(_energy_toxicity_threshold >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "INITIAL_ENERGY_AMOUNT_IN_CELLS") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _initial_energy_amount_in_cells;
          assert(_initial_energy_amount_in_cells >= 0.0);
        }
        
        /*------------------------------------------------------------------ cell */
        
        else if ( strcmp(words[0].c_str(), "MEMBRANE_PERMEABILITY") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _membrane_permeability;
          assert(_membrane_permeability >= 0.0);
          assert(_membrane_permeability <= 1.0);
        }
        
        /*------------------------------------------------------------------ population */
        
        else if ( strcmp(words[0].c_str(), "DEATH_PROBABILITY") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _death_probability;
          assert(_death_probability >= 0.0);
          assert(_death_probability <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "MIGRATION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _migration_rate;
          assert(_migration_rate >= 0.0);
          assert(_migration_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "HGT_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _hgt_rate;
          assert(_hgt_rate >= 0.0);
          assert(_hgt_rate <= 1.0);
        }
        
        /*------------------------------------------------------------------ environment */
        
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_INITIALIZATION_CYCLES") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _environment_properties.number_of_init_cycles;
          assert(_environment_properties.number_of_init_cycles >= 1);
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_SPECIES_TAG_RANGE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _environment_properties.species_tag_range.min >> _environment_properties.species_tag_range.max;
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_CONCENTRATION_RANGE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _environment_properties.concentration_range.min >> _environment_properties.concentration_range.max;
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_NUMBER_OF_SPECIES_RANGE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _environment_properties.number_of_species_range.min >> _environment_properties.number_of_species_range.max;
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_INTERACTION_SCHEME") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string scheme;
          flux >> param_name >> scheme;
          if (strcmp(scheme.c_str(), "NO_INTERACTION") == 0)
          {
            _environment_properties.interaction_scheme = NO_INTERACTION;
          }
          else if (strcmp(scheme.c_str(), "INTERACTION") == 0)
          {
            _environment_properties.interaction_scheme = INTERACTION;
          }
          else
          {
            std::cout << "Error : ENVIRONMENT_INTERACTION_SCHEME wrong value at line:\n " << line.c_str() << ".\n\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_RENEWAL_SCHEME") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string scheme;
          flux >> param_name >> scheme;
          if (strcmp(scheme.c_str(), "KEEP_MATTER") == 0)
          {
            _environment_properties.renewal_scheme = KEEP_MATTER;
          }
          else if (strcmp(scheme.c_str(), "CLEAR_MATTER") == 0)
          {
            _environment_properties.renewal_scheme = CLEAR_MATTER;
          }
          else
          {
            std::cout << "Error : ENVIRONMENT_RENEWAL_SCHEME wrong value at line:\n " << line.c_str() << ".\n\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_VARIATION_SCHEME") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string scheme;
          flux >> param_name >> scheme;
          if (strcmp(scheme.c_str(), "RANDOM") == 0)
          {
            _environment_properties.variation_scheme = RANDOM_SCHEME;
          }
          else if (strcmp(scheme.c_str(), "PERIODIC") == 0)
          {
            _environment_properties.variation_scheme = PERIODIC_SCHEME;
          }
          else if (strcmp(scheme.c_str(), "CYCLIC") == 0)
          {
            _environment_properties.variation_scheme = CYCLIC_SCHEME;
          }
          else
          {
            std::cout << "Error : ENVIRONMENT_VARIATION_SCHEME wrong value at line:\n " << line.c_str() << ".\n\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_LOCALIZATION_SCHEME") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string scheme;
          flux >> param_name >> scheme;
          if (strcmp(scheme.c_str(), "GLOBAL") == 0)
          {
            _environment_properties.localization_scheme = GLOBAL_LOCALIZATION;
          }
          else if (strcmp(scheme.c_str(), "RANDOM") == 0)
          {
            _environment_properties.localization_scheme = RANDOM_LOCALIZATION;
          }
          else if (strcmp(scheme.c_str(), "SPOT") == 0)
          {
            _environment_properties.localization_scheme = SPOT_LOCALIZATION;
          }
          else if (strcmp(scheme.c_str(), "CENTER") == 0)
          {
            _environment_properties.localization_scheme = CENTER_LOCALIZATION;
          }
          else
          {
            std::cout << "Error : ENVIRONMENT_LOCALIZATION_SCHEME wrong value at line:\n " << line.c_str() << ".\n\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_METABOLIC_SCHEME") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string scheme;
          flux >> param_name >> scheme;
          if (strcmp(scheme.c_str(), "UNIQUE") == 0)
          {
            _environment_properties.metabolic_scheme = UNIQUE_METABOLITE;
          }
          else if (strcmp(scheme.c_str(), "MULTIPLE") == 0)
          {
            _environment_properties.metabolic_scheme = MULTIPLE_METABOLITES;
          }
          else if (strcmp(scheme.c_str(), "BOUNDARIES") == 0)
          {
            _environment_properties.metabolic_scheme = BOUNDARIES;
          }
          else
          {
            std::cout << "Error : ENVIRONMENT_METABOLIC_SCHEME wrong value at line:\n " << line.c_str() << ".\n\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_INTRODUCTION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _environment_properties.introduction_rate;
          assert(_environment_properties.introduction_rate >= 0.0);
          assert(_environment_properties.introduction_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_DIFFUSION_COEFFICIENT") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _environment_properties.diffusion_coefficient;
          assert(_environment_properties.diffusion_coefficient >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_DEGRADATION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _environment_properties.degradation_rate;
          assert(_environment_properties.degradation_rate >= 0.0);
        }
        else
        {
          std::cout << "Unknown parameter '" << words[0] << "'. Exit.\n\n";
          exit(EXIT_FAILURE);
        }
      }
    }
    file.close();
    initialize_prngs();
  }
}

/**
 * \brief    Save in backup file
 * \details  --
 * \param    size_t backup_time
 * \return   \e void
 */
void Parameters::save( size_t backup_time )
{
  std::stringstream backup_file_name;
  backup_file_name << "./parameters/parameters_" << backup_time;
  gzFile backup_file = gzopen(backup_file_name.str().c_str(), "w");
  
  /*------------------------------------------------------------------ prng seed */
  
  gzwrite( backup_file, &_seed, sizeof(_seed) );
  
  /*------------------------------------------------------------------ parallel computing */
  
  gzwrite( backup_file, &_parallel_computing, sizeof(_parallel_computing) );
  
  /*------------------------------------------------------------------ simulation schemes */
  
  gzwrite( backup_file, &_energy_costs_scheme,          sizeof(_energy_costs_scheme) );
  gzwrite( backup_file, &_membrane_permeability_scheme, sizeof(_membrane_permeability_scheme) );
  gzwrite( backup_file, &_metabolic_inheritance_scheme, sizeof(_metabolic_inheritance_scheme) );
  gzwrite( backup_file, &_enzymatic_inheritance_scheme, sizeof(_enzymatic_inheritance_scheme) );
  gzwrite( backup_file, &_co_enzyme_activity_scheme,    sizeof(_co_enzyme_activity_scheme) );
  gzwrite( backup_file, &_score_scheme,                 sizeof(_score_scheme) );
  gzwrite( backup_file, &_selection_threshold,          sizeof(_selection_threshold) );
  
  /*------------------------------------------------------------------ space */
  
  gzwrite( backup_file, &_width,  sizeof(_width) );
  gzwrite( backup_file, &_height, sizeof(_height) );
  
  /*------------------------------------------------------------------ output */
  
  gzwrite( backup_file, &_simulation_backup_step,  sizeof(_simulation_backup_step) );
  gzwrite( backup_file, &_figures_generation_step, sizeof(_figures_generation_step) );
  
  /*------------------------------------------------------------------ genome */
  
  gzwrite( backup_file, &_load_genome_from_file, sizeof(_load_genome_from_file) );
  
  gzwrite( backup_file, &_metabolite_tag_initial_range.min, sizeof(_metabolite_tag_initial_range.min) );
  gzwrite( backup_file, &_metabolite_tag_initial_range.max, sizeof(_metabolite_tag_initial_range.max) );
  
  gzwrite( backup_file, &_binding_site_tag_initial_range.min, sizeof(_binding_site_tag_initial_range.min) );
  gzwrite( backup_file, &_binding_site_tag_initial_range.max, sizeof(_binding_site_tag_initial_range.max) );
  
  gzwrite( backup_file, &_co_enzyme_tag_initial_range.min, sizeof(_co_enzyme_tag_initial_range.min) );
  gzwrite( backup_file, &_co_enzyme_tag_initial_range.max, sizeof(_co_enzyme_tag_initial_range.max) );
  
  gzwrite( backup_file, &_transcription_factor_tag_initial_range.min, sizeof(_transcription_factor_tag_initial_range.min) );
  gzwrite( backup_file, &_transcription_factor_tag_initial_range.max, sizeof(_transcription_factor_tag_initial_range.max) );
  
  gzwrite( backup_file, &_transcription_factor_binding_window, sizeof(_transcription_factor_binding_window) );
  
  gzwrite( backup_file, &_initial_number_of_NC_units, sizeof(_initial_number_of_NC_units) );
  gzwrite( backup_file, &_initial_number_of_E_units,  sizeof(_initial_number_of_E_units) );
  gzwrite( backup_file, &_initial_number_of_TF_units, sizeof(_initial_number_of_TF_units) );
  gzwrite( backup_file, &_initial_number_of_BS_units, sizeof(_initial_number_of_BS_units) );
  gzwrite( backup_file, &_initial_number_of_P_units,  sizeof(_initial_number_of_P_units) );
  
  gzwrite( backup_file, &_point_mutation_rate, sizeof(_point_mutation_rate) );
  gzwrite( backup_file, &_duplication_rate,    sizeof(_duplication_rate) );
  gzwrite( backup_file, &_deletion_rate,       sizeof(_deletion_rate) );
  gzwrite( backup_file, &_translocation_rate,  sizeof(_translocation_rate) );
  gzwrite( backup_file, &_inversion_rate,      sizeof(_inversion_rate) );
  gzwrite( backup_file, &_transition_rate,     sizeof(_transition_rate) );
  gzwrite( backup_file, &_breakpoint_rate,     sizeof(_breakpoint_rate) );
  
  gzwrite( backup_file, &_substrate_tag_mutation_size,            sizeof(_substrate_tag_mutation_size) );
  gzwrite( backup_file, &_product_tag_mutation_size,              sizeof(_product_tag_mutation_size) );
  gzwrite( backup_file, &_kcat_mutation_size,                     sizeof(_kcat_mutation_size) );
  gzwrite( backup_file, &_kcat_km_ratio_mutation_size,            sizeof(_kcat_km_ratio_mutation_size) );
  gzwrite( backup_file, &_binding_site_tag_mutation_size,         sizeof(_binding_site_tag_mutation_size) );
  gzwrite( backup_file, &_co_enzyme_tag_mutation_size,            sizeof(_co_enzyme_tag_mutation_size) );
  gzwrite( backup_file, &_transcription_factor_tag_mutation_size, sizeof(_transcription_factor_tag_mutation_size) );
  gzwrite( backup_file, &_basal_expression_level_mutation_size,   sizeof(_basal_expression_level_mutation_size) );
  
  gzwrite( backup_file, &_mutation_of_mutation_rates, sizeof(_mutation_of_mutation_rates) );
  
  /*------------------------------------------------------------------ genetic regulation network */
  
  gzwrite( backup_file, &_genetic_regulation_network_timestep, sizeof(_genetic_regulation_network_timestep) );
  gzwrite( backup_file, &_hill_function_theta,                 sizeof(_hill_function_theta) );
  gzwrite( backup_file, &_hill_function_n,                     sizeof(_hill_function_n) );
  gzwrite( backup_file, &_protein_degradation_rate,            sizeof(_protein_degradation_rate) );
  
  /*------------------------------------------------------------------ metabolic network */
  
  gzwrite( backup_file, &_metabolism_timestep,                          sizeof(_metabolism_timestep) );
  gzwrite( backup_file, &_essential_metabolites_toxicity_threshold,     sizeof(_essential_metabolites_toxicity_threshold) );
  gzwrite( backup_file, &_non_essential_metabolites_toxicity_threshold, sizeof(_non_essential_metabolites_toxicity_threshold) );
  
  gzwrite( backup_file, &_initial_metabolites_amount_in_cells,          sizeof(_initial_metabolites_amount_in_cells) );
  
  /*------------------------------------------------------------------ energy */
  
  gzwrite( backup_file, &_energy_transcription_cost, sizeof(_energy_transcription_cost) );
  gzwrite( backup_file, &_energy_degradation_cost,   sizeof(_energy_degradation_cost) );
  gzwrite( backup_file, &_energy_enzymatic_cost,     sizeof(_energy_enzymatic_cost) );
  gzwrite( backup_file, &_energy_pumping_cost,       sizeof(_energy_pumping_cost) );
  
  gzwrite( backup_file, &_energy_dissipation_rate, sizeof(_energy_dissipation_rate) );
  
  gzwrite( backup_file, &_energy_toxicity_threshold, sizeof(_energy_toxicity_threshold) );
  
  gzwrite( backup_file, &_initial_energy_amount_in_cells, sizeof(_initial_energy_amount_in_cells) );
  
  gzwrite( backup_file, &_maximum_reaction_size, sizeof(_maximum_reaction_size) );
  
  /*------------------------------------------------------------------ cell */
  
  gzwrite( backup_file, &_membrane_permeability, sizeof(_membrane_permeability) );
  
  /*------------------------------------------------------------------ population */
  
  gzwrite( backup_file, &_death_probability, sizeof(_death_probability) );
  gzwrite( backup_file, &_migration_rate,    sizeof(_migration_rate) );
  gzwrite( backup_file, &_hgt_rate,          sizeof(_hgt_rate) );
  
  /*------------------------------------------------------------------ environment */
  
  gzwrite( backup_file, &_environment_properties.number_of_init_cycles, sizeof(_environment_properties.number_of_init_cycles) );
  
  gzwrite( backup_file, &_environment_properties.species_tag_range.min, sizeof(_environment_properties.species_tag_range.min) );
  gzwrite( backup_file, &_environment_properties.species_tag_range.max, sizeof(_environment_properties.species_tag_range.max) );
  
  gzwrite( backup_file, &_environment_properties.concentration_range.min, sizeof(_environment_properties.concentration_range.min) );
  gzwrite( backup_file, &_environment_properties.concentration_range.max, sizeof(_environment_properties.concentration_range.max) );
  
  gzwrite( backup_file, &_environment_properties.number_of_species_range.min, sizeof(_environment_properties.number_of_species_range.min) );
  gzwrite( backup_file, &_environment_properties.number_of_species_range.max, sizeof(_environment_properties.number_of_species_range.max) );
  
  gzwrite( backup_file, &_environment_properties.interaction_scheme,  sizeof(_environment_properties.interaction_scheme) );
  gzwrite( backup_file, &_environment_properties.renewal_scheme,      sizeof(_environment_properties.renewal_scheme) );
  gzwrite( backup_file, &_environment_properties.variation_scheme,    sizeof(_environment_properties.variation_scheme) );
  gzwrite( backup_file, &_environment_properties.localization_scheme, sizeof(_environment_properties.localization_scheme) );
  gzwrite( backup_file, &_environment_properties.metabolic_scheme,    sizeof(_environment_properties.metabolic_scheme) );
  
  gzwrite( backup_file, &_environment_properties.introduction_rate,     sizeof(_environment_properties.introduction_rate) );
  gzwrite( backup_file, &_environment_properties.diffusion_coefficient, sizeof(_environment_properties.diffusion_coefficient) );
  gzwrite( backup_file, &_environment_properties.degradation_rate,      sizeof(_environment_properties.degradation_rate) );
  
  gzclose(backup_file);
  
  /*------------------------------------------------------------------ prngs */
  
  std::stringstream prng_file_name;
  prng_file_name << "./prng/prng_" << backup_time;
  FILE * prng_file = fopen(prng_file_name.str().c_str(), "w");
  _simulation_prng->save(prng_file);
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    _population_prng[pos]->save(prng_file);
  }
  _environment_prng->save(prng_file);
  fclose(prng_file);
}

/**
 * \brief    Write parameters in a parameter file
 * \details  --
 * \param    std::string filename
 * \return   \e void
 */
void Parameters::write( std::string filename )
{
  std::ofstream file(filename.c_str(), std::ios::out | std::ios::trunc);
  file << "\n";
  file << "#***************************************************************************\n";
  file << "# Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon\n";
  file << "# Web: https://github.com/charlesrocabert/Evo2Sim\n";
  file << "#\n";
  file << "# This program is free software: you can redistribute it and/or modify\n";
  file << "# it under the terms of the GNU General Public License as published by\n";
  file << "# the Free Software Foundation, either version 3 of the License, or\n";
  file << "# (at your option) any later version.\n";
  file << "#\n";
  file << "# This program is distributed in the hope that it will be useful,\n";
  file << "# but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
  file << "# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n";
  file << "# GNU General Public License for more details.\n";
  file << "#\n";
  file << "# You should have received a copy of the GNU General Public License\n";
  file << "# along with this program.  If not, see <http://www.gnu.org/licenses/>.\n";
  file << "#***************************************************************************\n";
  file << "\n";
  file << "########################################################\n";
  file << "# PSEUDORANDOM NUMBERS GENERATOR\n";
  file << "########################################################\n";
  file << "\n";
  file << "SEED  " << _seed << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# PARALLEL COMPUTING\n";
  file << "########################################################\n";
  file << "\n";
  file << "PARALLEL_COMPUTING  ";
  if (!_parallel_computing)
  {
    file << "NO ";
  }
  else
  {
    file << "YES";
  }
  file << " (YES, NO)\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# SIMULATION SCHEMES\n";
  file << "########################################################\n";
  file << "\n";
  
  file << "ENERGY_COSTS_SCHEME           ";
  if (!_energy_costs_scheme)
  {
    file << "NO ";
  }
  else
  {
    file << "YES";
  }
  file << " (YES, NO)\n";
  
  file << "MEMBRANE_PERMEABILITY_SCHEME  ";
  if (!_membrane_permeability_scheme)
  {
    file << "NO ";
  }
  else
  {
    file << "YES";
  }
  file << " (YES, NO)\n";
  
  file << "METABOLIC_INHERITANCE_SCHEME  ";
  if (!_metabolic_inheritance_scheme)
  {
    file << "NO ";
  }
  else
  {
    file << "YES";
  }
  file << " (YES, NO)\n";
  
  file << "ENZYMATIC_INHERITANCE_SCHEME  ";
  if (!_enzymatic_inheritance_scheme)
  {
    file << "NO ";
  }
  else
  {
    file << "YES";
  }
  file << " (YES, NO)\n";
  
  file << "CO_ENZYME_ACTIVITY_SCHEME     ";
  if (!_co_enzyme_activity_scheme)
  {
    file << "NO ";
  }
  else
  {
    file << "YES";
  }
  file << " (YES, NO)\n";
  
  file << "\n";
  file << "SCORE_SCHEME         ";
  switch (_score_scheme)
  {
    case ESSENTIAL_METABOLITES_SUM:
      file << "SUM (SUM, SUM_MINUS_DEV, COMBINATORIAL)\n";
      break;
    case ESSENTIAL_METABOLITES_SUM_MINUS_DEVIATION:
      file << "SUM_MINUS_DEV (SUM, SUM_MINUS_DEV, COMBINATORIAL)\n";
      break;
    case ESSENTIAL_METABOLITES_COMBINATORIAL_CONTRIBUTION:
      file << "COMBINATORIAL (SUM, SUM_MINUS_DEV, COMBINATORIAL)\n";
      break;
  }
  file << "SELECTION_THRESHOLD  " << _selection_threshold << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# SPACE\n";
  file << "########################################################\n";
  file << "\n";
  file << "WIDTH   " << _width << "\n";
  file << "HEIGHT  " << _height << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# OUTPUT\n";
  file << "########################################################\n";
  file << "\n";
  file << "SIMULATION_BACKUP_STEP   " << _simulation_backup_step << "\n";
  file << "FIGURES_GENERATION_STEP  " << _figures_generation_step << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# GENOME\n";
  file << "########################################################\n";
  file << "\n";
  file << "LOAD_GENOME_FROM_FILE  ";
  if (!_load_genome_from_file)
  {
    file << "NO ";
  }
  else
  {
    file << "YES";
  }
  file << " (YES, NO)\n";
  file << "\n";
  file << "METABOLITE_TAG_INITIAL_RANGE            " << _metabolite_tag_initial_range.min << " " << _metabolite_tag_initial_range.max << "\n";
  file << "BINDING_SITE_TAG_INITIAL_RANGE          " << _binding_site_tag_initial_range.min << " " << _binding_site_tag_initial_range.max << "\n";
  
  file << "CO_ENZYME_TAG_INITIAL_RANGE             " << _co_enzyme_tag_initial_range.min << " " << _co_enzyme_tag_initial_range.max << "\n";
  file << "TRANSCRIPTION_FACTOR_TAG_INITIAL_RANGE  " << _transcription_factor_tag_initial_range.min << " " << _transcription_factor_tag_initial_range.max << "\n";
  file << "\n";
  file << "TRANSCRIPTION_FACTOR_BINDING_WINDOW  " << _transcription_factor_binding_window << "\n";
  file << "\n";
  file << "INITIAL_NUMBER_OF_NON_CODING_UNITS            " << _initial_number_of_NC_units << "\n";
  file << "INITIAL_NUMBER_OF_ENZYME_UNITS                " << _initial_number_of_E_units << "\n";
  file << "INITIAL_NUMBER_OF_TRANSCRIPTION_FACTOR_UNITS  " << _initial_number_of_TF_units << "\n";
  file << "INITIAL_NUMBER_OF_BINDING_SITE_UNITS          " << _initial_number_of_BS_units << "\n";
  file << "INITIAL_NUMBER_OF_PROMOTER_UNITS              " << _initial_number_of_P_units << "\n";
  file << "\n";
  file << "POINT_MUTATION_RATE  " << _point_mutation_rate << "\n";
  file << "DUPLICATION_RATE     " << _duplication_rate << "\n";
  file << "DELETION_RATE        " << _deletion_rate << "\n";
  file << "TRANSLOCATION_RATE   " << _translocation_rate << "\n";
  file << "INVERSION_RATE       " << _inversion_rate << "\n";
  file << "TRANSITION_RATE      " << _transition_rate << "\n";
  file << "BREAKPOINT_RATE      " << _breakpoint_rate << "\n";
  file << "\n";
  file << "SUBSTRATE_TAG_MUTATION_SIZE             " << _substrate_tag_mutation_size << "\n";
  file << "PRODUCT_TAG_MUTATION_SIZE               " << _product_tag_mutation_size << "\n";
  file << "KCAT_MUTATION_SIZE                      " << _kcat_mutation_size << "\n";
  file << "KCAT_KM_RATIO_MUTATION_SIZE             " << _kcat_km_ratio_mutation_size << "\n";
  file << "BINDING_SITE_TAG_MUTATION_SIZE          " << _binding_site_tag_mutation_size << "\n";
  file << "CO_ENZYME_TAG_MUTATION_SIZE             " << _co_enzyme_tag_mutation_size << "\n";
  file << "TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE  " << _transcription_factor_tag_mutation_size << "\n";
  file << "BASAL_EXPRESSION_LEVEL_MUTATION_SIZE    " << _basal_expression_level_mutation_size << "\n";
  file << "\n";
  file << "MUTATION_OF_MUTATION_RATES  " << _mutation_of_mutation_rates << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# GENETIC REGULATION NETWORK\n";
  file << "########################################################\n";
  file << "\n";
  file << "GENETIC_REGULATION_NETWORK_TIMESTEP  " << _genetic_regulation_network_timestep << "\n";
  file << "\n";
  file << "HILL_FUNCTION_THETA       " << _hill_function_theta << "\n";
  file << "HILL_FUNCTION_N           " << _hill_function_n << "\n";
  file << "PROTEIN_DEGRADATION_RATE  " << _protein_degradation_rate << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# METABOLIC NETWORK\n";
  file << "########################################################\n";
  file << "\n";
  file << "METABOLISM_TIMESTEP  " << _metabolism_timestep << "\n";
  file << "\n";
  file << "ESSENTIAL_METABOLITES_TOXICITY_THRESHOLD      " << _essential_metabolites_toxicity_threshold << "\n";
  file << "NON_ESSENTIAL_METABOLITES_TOXICITY_THRESHOLD  " << _non_essential_metabolites_toxicity_threshold << "\n";
  file << "\n";
  file << "INITIAL_METABOLITES_AMOUNT_IN_CELLS  " << _initial_metabolites_amount_in_cells << "\n";
  file << "\n";
  file << "MAXIMUM_REACTION_SIZE  " << _maximum_reaction_size << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# ENERGY\n";
  file << "########################################################\n";
  file << "\n";
  file << "ENERGY_TRANSCRIPTION_COST  " << _energy_transcription_cost << "\n";
  file << "ENERGY_DEGRADATION_COST    " << _energy_degradation_cost << "\n";
  file << "ENERGY_ENZYMATIC_COST      " << _energy_enzymatic_cost << "\n";
  file << "ENERGY_PUMPING_COST        " << _energy_pumping_cost << "\n";
  file << "\n";
  file << "ENERGY_DISSIPATION_RATE  " << _energy_dissipation_rate << "\n";
  file << "\n";
  file << "ENERGY_TOXICITY_THRESHOLD  " << _energy_toxicity_threshold << "\n";
  file << "\n";
  file << "INITIAL_ENERGY_AMOUNT_IN_CELLS  " << _initial_energy_amount_in_cells << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# CELL\n";
  file << "########################################################\n";
  file << "\n";
  file << "MEMBRANE_PERMEABILITY  " << _membrane_permeability << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# POPULATION\n";
  file << "########################################################\n";
  file << "\n";
  file << "DEATH_PROBABILITY  " << _death_probability << "\n";
  file << "MIGRATION_RATE     " << _migration_rate << "\n";
  file << "HGT_RATE           " << _hgt_rate << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# ENVIRONMENT\n";
  file << "########################################################\n";
  file << "\n";
  file << "ENVIRONMENT_INITIALIZATION_CYCLES    " << _environment_properties.number_of_init_cycles << "\n";
  file << "ENVIRONMENT_SPECIES_TAG_RANGE        " << _environment_properties.species_tag_range.min << " " << _environment_properties.species_tag_range.max << "\n";
  file << "ENVIRONMENT_CONCENTRATION_RANGE      " << _environment_properties.concentration_range.min << " " << _environment_properties.concentration_range.max << "\n";
  file << "ENVIRONMENT_NUMBER_OF_SPECIES_RANGE  " << _environment_properties.number_of_species_range.min << " " << _environment_properties.number_of_species_range.max << "\n";
  file << "\n";
  file << "ENVIRONMENT_INTERACTION_SCHEME   ";
  switch (_environment_properties.interaction_scheme)
  {
    case NO_INTERACTION:
      file << "NO_INTERACTION (NO_INTERACTION/INTERACTION)\n";
      break;
    case INTERACTION:
      file << "INTERACTION (NO_INTERACTION/INTERACTION)\n";
      break;
  }
  file << "ENVIRONMENT_RENEWAL_SCHEME       ";
  switch (_environment_properties.renewal_scheme)
  {
    case KEEP_MATTER:
      file << "KEEP_MATTER (KEEP_MATTER/CLEAR_MATTER)\n";
      break;
    case CLEAR_MATTER:
      file << "CLEAR_MATTER (KEEP_MATTER/CLEAR_MATTER)\n";
      break;
  }
  file << "ENVIRONMENT_VARIATION_SCHEME     ";
  switch (_environment_properties.variation_scheme)
  {
    case RANDOM_SCHEME:
      file << "RANDOM (RANDOM/PERIODIC/CYCLIC)\n";
      break;
    case PERIODIC_SCHEME:
      file << "PERIODIC (RANDOM/PERIODIC/CYCLIC)\n";
      break;
    case CYCLIC_SCHEME:
      file << "CYCLIC (RANDOM/PERIODIC/CYCLIC)\n";
      break;
  }
  file << "ENVIRONMENT_LOCALIZATION_SCHEME  ";
  switch (_environment_properties.localization_scheme)
  {
    case GLOBAL_LOCALIZATION:
      file << "GLOBAL (GLOBAL/RANDOM/SPOT/CENTER)\n";
      break;
    case RANDOM_LOCALIZATION:
      file << "RANDOM (GLOBAL/RANDOM/SPOT/CENTER)\n";
      break;
    case SPOT_LOCALIZATION:
      file << "SPOT (GLOBAL/RANDOM/SPOT/CENTER)\n";
      break;
    case CENTER_LOCALIZATION:
      file << "CENTER (GLOBAL/RANDOM/SPOT/CENTER)\n";
      break;
  }
  file << "ENVIRONMENT_METABOLIC_SCHEME     ";
  switch (_environment_properties.metabolic_scheme)
  {
    case UNIQUE_METABOLITE:
      file << "UNIQUE (UNIQUE/MULTIPLE/BOUNDARIES)\n";
      break;
    case MULTIPLE_METABOLITES:
      file << "MULTIPLE (UNIQUE/MULTIPLE/BOUNDARIES)\n";
      break;
    case BOUNDARIES:
      file << "BOUNDARIES (UNIQUE/MULTIPLE/BOUNDARIES)\n";
      break;
  }
  file << "\n";
  file << "ENVIRONMENT_INTRODUCTION_RATE      " << _environment_properties.introduction_rate << "\n";
  file << "ENVIRONMENT_DIFFUSION_COEFFICIENT  " << _environment_properties.diffusion_coefficient << "\n";
  file << "ENVIRONMENT_DEGRADATION_RATE       " << _environment_properties.degradation_rate << "\n";
  file << "\n";
  file << "\n";
  file.close();
}

/**
 * \brief    Draw initial metabolite tag in its range
 * \details  --
 * \param    void
 * \return   \e int
 */
int Parameters::draw_initial_metabolite_tag( void )
{
  return _simulation_prng->uniform(_metabolite_tag_initial_range.min, _metabolite_tag_initial_range.max);
}

/**
 * \brief    Draw initial binding site tag in its range
 * \details  --
 * \param    void
 * \return   \e int
 */
int Parameters::draw_initial_binding_site_tag( void )
{
  return _simulation_prng->uniform(_binding_site_tag_initial_range.min, _binding_site_tag_initial_range.max);
}

/**
 * \brief    Draw initial co-enzyme tag in its range
 * \details  --
 * \param    void
 * \return   \e int
 */
int Parameters::draw_initial_co_enzyme_tag( void )
{
  return _simulation_prng->uniform(_co_enzyme_tag_initial_range.min, _co_enzyme_tag_initial_range.max);
}

/**
 * \brief    Draw initial transcription factor tag in its range
 * \details  --
 * \param    void
 * \return   \e int
 */
int Parameters::draw_initial_transcription_factor_tag( void )
{
  return _simulation_prng->uniform(_transcription_factor_tag_initial_range.min, _transcription_factor_tag_initial_range.max);
}

/**
 * \brief    Draw environmental number of species in its range
 * \details  --
 * \param    void
 * \return   \e int
 */
int Parameters::draw_environment_number_of_species( void )
{
  return _simulation_prng->uniform(_environment_properties.number_of_species_range.min, _environment_properties.number_of_species_range.max);
}

/**
 * \brief    Draw environmental species tag in its range
 * \details  --
 * \param    void
 * \return   \e int
 */
int Parameters::draw_environment_species_tag( void )
{
  return _simulation_prng->uniform(_environment_properties.species_tag_range.min, _environment_properties.species_tag_range.max);
}

/**
 * \brief    Draw environmental concentration in its range
 * \details  --
 * \param    void
 * \return   \e double
 */
double Parameters::draw_environment_concentration( void )
{
  return _simulation_prng->uniform()*(_environment_properties.concentration_range.max-_environment_properties.concentration_range.min)+_environment_properties.concentration_range.min;
}

/**
 * \brief    Draw genetic unit's attributes at random
 * \details  --
 * \param    genetic_unit& unit
 * \param    genetic_unit_type type
 * \return   \e void
 */
void Parameters::draw_random_genetic_unit( genetic_unit& unit, genetic_unit_type type )
{
  /*------------------------------------------------------------------ Global attributes */
  
  unit.type              = type;
  unit.identifier        = 0;
  unit.parent_identifier = 0;
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  unit.s             = draw_initial_metabolite_tag();
  unit.s             = (unit.s > 0 ? unit.s : 1);
  unit.p             = draw_initial_metabolite_tag();
  unit.p             = (unit.p > 0 ? unit.p : 1);
  unit.kcat          = pow(10.0, _simulation_prng->uniform()*(KCAT_MAX_LOG-KCAT_MIN_LOG)+KCAT_MIN_LOG);
  unit.kcat          = unit.kcat*(_simulation_prng->uniform() < 0.5 ? -1.0 : 1.0);
  unit.kcat_km_ratio = pow(10.0, _simulation_prng->uniform()*(KCAT_KM_RATIO_MAX_LOG-KCAT_KM_RATIO_MIN_LOG)+KCAT_KM_RATIO_MIN_LOG);
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  unit.BS_tag         = draw_initial_binding_site_tag();
  unit.coE_tag        = draw_initial_co_enzyme_tag();
  unit.coE_tag        = (unit.coE_tag > 0 ? unit.coE_tag : 1);
  unit.free_activity  = (_simulation_prng->uniform() < 0.5 ? true : false);
  unit.bound_activity = (_simulation_prng->uniform() < 0.5 ? true : false);
  unit.binding_window = get_transcription_factor_binding_window();
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  unit.TF_tag = draw_initial_transcription_factor_tag();
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  unit.basal_expression_level = get_simulation_prng()->uniform();
}

/**
 * \brief    Load the initial genome from file
 * \details  --
 * \param    std::vector<genetic_unit>& sequence
 * \return   \e void
 */
void Parameters::load_genome_from_file( std::vector<genetic_unit>& sequence )
{
  std::ifstream file("initial_genome.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: initial_genome.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    sequence.clear();
    std::string line;
    std::string param_name;
    while(getline(file, line))
    {
      std::vector<std::string> words;
      if(parse_line(&words, line))
      {
        /*****************************/
        /* NON CODING UNIT           */
        /*****************************/
        if ( strcmp(words[0].c_str(), "NC") == 0 )
        {
          size_t nb_NC = 0;
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> nb_NC;
          for (size_t i = 0; i < nb_NC; i++)
          {
            genetic_unit unit;
            draw_random_genetic_unit(unit, NON_CODING);
            sequence.push_back(unit);
          }
        }
        /*****************************/
        /* ENZYME UNIT               */
        /*****************************/
        else if ( strcmp(words[0].c_str(), "E") == 0 )
        {
          genetic_unit unit;
          draw_random_genetic_unit(unit, ENZYME);
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> unit.s >> unit.p >> unit.kcat >> unit.kcat_km_ratio;
          assert(!(unit.s < 0 && unit.p < 0));
          assert(unit.kcat >= KCAT_MIN_LOG);
          assert(unit.kcat <= KCAT_MAX_LOG);
          assert(unit.kcat_km_ratio >= KCAT_KM_RATIO_MIN_LOG);
          assert(unit.kcat_km_ratio <= KCAT_KM_RATIO_MAX_LOG);
          if (unit.s < 0 && unit.p > 0)
          {
            unit.s = -unit.s;
            unit.kcat = pow(10.0, unit.kcat);
          }
          else if (unit.s > 0 && unit.p < 0)
          {
            unit.p = -unit.p;
            unit.kcat = -pow(10.0, unit.kcat);
          }
          unit.kcat_km_ratio = pow(10.0, unit.kcat_km_ratio);
          sequence.push_back(unit);
        }
        /*****************************/
        /* TRANSCRIPTION FACTOR UNIT */
        /*****************************/
        else if ( strcmp(words[0].c_str(), "TF") == 0 )
        {
          genetic_unit unit;
          draw_random_genetic_unit(unit, TRANSCRIPTION_FACTOR);
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          std::string buffer1;
          std::string buffer2;
          flux >> param_name >> unit.BS_tag >> unit.coE_tag >> buffer1 >> buffer2 >> unit.kcat >> unit.kcat_km_ratio;
          assert(unit.kcat >= KCAT_MIN_LOG);
          assert(unit.kcat <= KCAT_MAX_LOG);
          assert(unit.kcat_km_ratio >= KCAT_KM_RATIO_MIN_LOG);
          assert(unit.kcat_km_ratio <= KCAT_KM_RATIO_MAX_LOG);
          unit.kcat = pow(10.0, unit.kcat);
          unit.kcat_km_ratio = pow(10.0, unit.kcat_km_ratio);
          if (strcmp(words[3].c_str(), "false") == 0)
          {
            unit.free_activity = false;
          }
          else if (strcmp(words[3].c_str(), "true") == 0)
          {
            unit.free_activity = true;
          }
          else
          {
            std::cout << "Error : TF wrong value at line:\n " << line.c_str() << ".\n\n";
            exit(EXIT_FAILURE);
          }
          if (strcmp(words[4].c_str(), "false") == 0)
          {
            unit.bound_activity = false;
          }
          else if (strcmp(words[4].c_str(), "true") == 0)
          {
            unit.bound_activity = true;
          }
          else
          {
            std::cout << "Error : TF wrong value at line:\n " << line.c_str() << ".\n\n";
            exit(EXIT_FAILURE);
          }
          sequence.push_back(unit);
        }
        /*****************************/
        /* BINDING SITE UNIT         */
        /*****************************/
        else if ( strcmp(words[0].c_str(), "BS") == 0 )
        {
          genetic_unit unit;
          draw_random_genetic_unit(unit, BINDING_SITE);
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> unit.TF_tag;
          sequence.push_back(unit);
        }
        /*****************************/
        /* PROMOTER UNIT             */
        /*****************************/
        else if ( strcmp(words[0].c_str(), "P") == 0 )
        {
          genetic_unit unit;
          draw_random_genetic_unit(unit, PROMOTER);
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> unit.basal_expression_level;
          assert(unit.basal_expression_level >= 0.0);
          assert(unit.basal_expression_level <= 1.0);
          sequence.push_back(unit);
        }
        else
        {
          std::cout << "Unknown parameter '" << words[0] << "'. Exit.\n\n";
          exit(EXIT_FAILURE);
        }
      }
    }
    file.close();
  }
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Parse a line in words
 * \details  --
 * \param    std::vector<std::string>* words
 * \param    std::string* line
 * \return   \e void
 */
bool Parameters::parse_line( std::vector<std::string>* words, std::string line )
{
  if (line.size() == 0)
  {
    return false;
  }
  if (line[0] == '#' || line[0] == '\n')
  {
    return false;
  }
  words->clear();
  bool new_word   = false;
  bool word_saved = false;
  int pos       = -1;
  int length    = -1;
  for (int i = 0; i < (int)line.length(); i++)
  {
    /* if find a new character, start word saving */
    if (line[i] != ' ' && line[i] != '\n' && !new_word)
    {
      new_word   = true;
      word_saved = false;
      pos      = i;
      length   = 0;
    }
    /* else stop it */
    else if (line[i] == ' ' || line[i] == '\n')
    {
      new_word = false;
    }
    /* if a new word is found, save word length */
    if (new_word)
    {
      length++;
    }
    if (!new_word && !word_saved && pos >= 0 && length > 0)
    {
      words->push_back(line.substr(pos, length));
      word_saved = true;
    }
  }
  if (new_word && !word_saved && pos >= 0 && length > 0)
  {
    words->push_back(line.substr(pos, length));
  }
  return true;
}

/**
 * \brief    Build the list of the prime numbers until maximum value
 * \details  Use the 'crible' method
 * \param    size_t maximum
 * \return   \e void
 */
void Parameters::build_prime_numbers_list( size_t maximum )
{
  _prime_numbers = NULL;
  _prime_numbers = new int[maximum];
  for (size_t i = 1; i <= maximum; i++)
  {
    _prime_numbers[i-1] = 1;
  }
  _prime_numbers[0] = 0;
  for (size_t i = 2; i <= maximum; i++)
  {
    size_t multiple = 2*i;
    while (multiple <= maximum)
    {
      _prime_numbers[multiple-1] = 0;
      multiple += i;
    }
  }
}

/**
 * \brief    Initialize the prngs seeds based on the simulation prng
 * \details  --
 * \param    void
 * \return   \e void
 */
void Parameters::initialize_prngs( void )
{
  /*----------------------*/
  /* 1) Free prngs memory */
  /*----------------------*/
  if (_simulation_prng != NULL)
  {
    delete _simulation_prng;
    _simulation_prng = NULL;
  }
  if (_population_prng != NULL)
  {
    for (size_t pos = 0; pos < _width*_height; pos++)
    {
      delete _population_prng[pos];
      _population_prng[pos] = NULL;
    }
  }
  delete[] _population_prng;
  if (_environment_prng != NULL)
  {
    delete _environment_prng;
    _environment_prng = NULL;
  }
  
  /*----------------------*/
  /* 2) Create prngs      */
  /*----------------------*/
  _simulation_prng = new Prng();
  _simulation_prng->set_seed(_seed);
  _population_prng = NULL;
  _population_prng = new Prng*[_width*_height];
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    _population_prng[pos] = new Prng();
    _population_prng[pos]->set_seed((unsigned long int)_simulation_prng->uniform(1, MAXIMUM_SEED));
  }
  _environment_prng = new Prng();
  _environment_prng->set_seed((unsigned long int)_simulation_prng->uniform(1, MAXIMUM_SEED));
}
