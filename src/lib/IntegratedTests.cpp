
/**
 * \file      IntegratedTests.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     IntegratedTests class definition
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

#include "IntegratedTests.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters1
 * \param    Parameters* parameters2
 * \return   \e void
 */
IntegratedTests::IntegratedTests( Parameters* parameters1, Parameters* parameters2 )
{
  _parameters1 = parameters1;
  _parameters2 = parameters2;
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
IntegratedTests::~IntegratedTests( void )
{
  /* NOTHING TODO */
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Run integrated tests
 * \details  --
 * \param    std::string parameters_filename
 * \param    size_t number_of_tests
 * \param    size_t number_of_steps
 * \param    bool random_seed
 * \param    bool random_parameters
 * \return   \e void
 */
void IntegratedTests::run_integrated_tests( size_t number_of_tests, size_t number_of_steps, bool random_seed, bool random_parameters )
{
  printf("\n************ Run integrated tests ************\n");
  
  /*-------------------------------------------------------------*/
  /* 1) Modify backup steps                                      */
  /*-------------------------------------------------------------*/
  _parameters1->set_simulation_backup_step(500);
  _parameters2->set_simulation_backup_step(500);
  
  /*-------------------------------------------------------------*/
  /* 2) Create subfolders in order to run simulations separately */
  /*-------------------------------------------------------------*/
  mkdir("simulation1", 0777);
  mkdir("simulation2", 0777);
  
  /*-------------------------------------------------------------*/
  /* 3) Run tests                                                */
  /*-------------------------------------------------------------*/
  Prng* prng = new Prng();
  prng->set_seed((unsigned long int)time(NULL));
  
  for (size_t test = 0; test < number_of_tests; test++)
  {
    printf("\n-----------------------\n> Launching test %lu\n-----------------------\n", test+1);
    unsigned long int seed = _parameters1->get_seed();
    if (random_seed)
    {
      seed = prng->uniform(1, MAXIMUM_SEED);
    }
    if (random_parameters)
    {
      generate_random_parameters(_parameters1, prng);
      delete _parameters2;
      _parameters2 = new Parameters(*_parameters1);
    }
    test_simulation_creation(seed);
    
    size_t backup_step = 0;
    size_t step        = 1;
    bool   run         = true;
    while (backup_step < 500*number_of_steps && run)
    {
      printf("\n    > Backup step %lu (step %lu)\n", backup_step, step);
      run = test_simulation_execution(backup_step, 500);
      backup_step += 500;
      step++;
    }
    if (!run)
    {
      printf("    > Population extincted.\n");
    }
  }
  
  /*-------------------------------------------------------------*/
  /* 4) Remove subfolders                                        */
  /*-------------------------------------------------------------*/
  int success = 0;
  success = system("rm -rf simulation1");
  success = system("rm -rf simulation2");
  (void)success;
  
  printf("\n************ End of integrated tests ************\n\n");
}

/**
 * \brief    Create a simulation
 * \details  Create a simulation with the seed given in parameter
 * \param    Parameters* parameters
 * \param    unsigned long int seed
 * \return   \e Simulation*
 */
Simulation* IntegratedTests::create_simulation( Parameters* parameters, unsigned long int seed )
{
  /*----------------------------------*/
  /* 1) Create simulation folders     */
  /*----------------------------------*/
  int success = 0;
  success = system("rm -rf statistics");
  mkdir("statistics", 0777);
  success = system("rm -rf environment");
  mkdir("environment", 0777);
  success = system("rm -rf population");
  mkdir("population", 0777);
  success = system("rm -rf trophic_network");
  mkdir("trophic_network", 0777);
  success = system("rm -rf tree");
  mkdir("tree", 0777);
  success = system("rm -rf parameters");
  mkdir("parameters", 0777);
  success = system("rm -rf prng");
  mkdir("prng", 0777);
  success = system("rm -rf figures");
  mkdir("figures", 0777);
  (void)success;
  
  /*----------------------------------*/
  /* 2) Initialize PRNG seed          */
  /*----------------------------------*/
  parameters->set_seed(seed);
  
  /*----------------------------------*/
  /* 3) Load and initialize           */
  /*    simulation                    */
  /*----------------------------------*/
  Simulation* simulation = new Simulation(parameters);
  simulation->initialize(NEW_SIMULATION);
  
  /*----------------------------------*/
  /* 4) Save simulation               */
  /*----------------------------------*/
  simulation->save_simulation();
  
  /*----------------------------------*/
  /* 5) Return the simulation         */
  /*----------------------------------*/
  return simulation;
}

/**
 * \brief    Test simulation creation
 * \details  --
 * \param    unsigned long int seed
 * \return   \e Simulation*
 */
void IntegratedTests::test_simulation_creation( unsigned long int seed )
{
  /*
   A) Create both simulations (simulation1 and simulation2)
   B) Compare simulations
   */
  
  /*----------------------------------------------------------*/
  /* A) Create both simulations (simulation1 and simulation2) */
  /*----------------------------------------------------------*/
  printf("    > Create simulations with seed %ld ...\n", seed);
  int success = 0;
  success = chdir("./simulation1/");
  Simulation* simulation1 = create_simulation(_parameters1, seed);
  _parameters1->write("parameters.txt");
  success = chdir("../simulation2/");
  Simulation* simulation2 = create_simulation(_parameters2, seed);
  _parameters2->write("parameters.txt");
  success = chdir("..");
  (void)success;
  
  /*----------------------------------------------------------*/
  /* B) Test simulations equality                             */
  /*----------------------------------------------------------*/
  Simulation_isEqualTo(simulation1, simulation2);
  
  delete simulation1;
  simulation1 = NULL;
  delete simulation2;
  simulation2 = NULL;
}

/**
 * \brief    Test simulation creation
 * \details  Return true if the population survived, false else
 * \param    size_t backup_time
 * \param    size_t simulation_time
 * \return   \e bool
 */
bool IntegratedTests::test_simulation_execution( size_t backup_time, size_t simulation_time )
{
  /*
   A) Reload simulations from backup and compare
   B) Run both simulations for 500 timesteps and compare at each timestep
   C) Compare at the end of the run
   */
  
  /*-------------------------------------------------------------------------*/
  /* A) Run both simulations for 500 timesteps (simulation1 and simulation2) */
  /*-------------------------------------------------------------------------*/
  printf("      Load simulations from backup ...\n");
  int success = 0;
  success = chdir("./simulation1/");
  Simulation* simulation1 = new Simulation(_parameters1, backup_time, true);
  simulation1->initialize( FROM_BACKUP );
  success = chdir("../simulation2/");
  Simulation* simulation2 = new Simulation(_parameters2, backup_time, true);
  simulation2->initialize( FROM_BACKUP );
  success = chdir("..");
  
  Simulation_isEqualTo(simulation1, simulation2);
  
  printf("      Run simulations for 500 timesteps ...\n");
  double t    = (double)backup_time;
  double TIME = (double)backup_time+simulation_time;
  bool   run  = true;
  
  while (t < TIME && run)
  {
    t += 1.0;
    if ((int)t%100 == 0)
    {
      printf("      Elapsed time %d\n", (int)t);
    }
    success = chdir("./simulation1/");
    simulation1->update();
    success = chdir("../simulation2/");
    simulation2->update();
    success = chdir("..");
    
    if (_parameters1->get_simulation_backup_step() > 0 && simulation1->get_population()->get_time() % _parameters1->get_simulation_backup_step() == 0)
    {
      success = chdir("./simulation1/");
      simulation1->save_simulation();
      success = chdir("../simulation2/");
      simulation2->save_simulation();
      success = chdir("..");
      assert(simulation1->get_population()->get_time() == simulation2->get_population()->get_time());
    }
    
    if (simulation1->get_population()->get_population_size() == 0)
    {
      assert(simulation1->get_population()->get_population_size() == simulation2->get_population()->get_population_size());
      run = false;
    }
    else
    {
      success = chdir("./simulation1/");
      simulation1->get_statistics()->write_stats();
      success = chdir("../simulation2/");
      simulation2->get_statistics()->write_stats();
      success = chdir("..");
    }
    if (run)
    {
      Simulation_isEqualTo(simulation1, simulation2);
    }
  }
  
  success = chdir("./simulation1/");
  simulation1->get_statistics()->close_files();
  success = chdir("../simulation2/");
  simulation2->get_statistics()->close_files();
  success = chdir("..");
  (void)success;
  
  Simulation_isEqualTo(simulation1, simulation2);
  
  delete simulation1;
  simulation1 = NULL;
  delete simulation2;
  simulation2 = NULL;
  
  return run;
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Test Parameters struct equality
 * \details  --
 * \param    Parameters* obj1
 * \param    Parameters* obj2
 * \return   \e void
 */
void IntegratedTests::Parameters_isEqualTo( Parameters* obj1, Parameters* obj2 )
{
  /*------------------------------------------------------------------ prng seed */
  
  assert(obj1->get_seed() == obj2->get_seed());
  
  /*------------------------------------------------------------------ parallel computing */
  
  assert(obj1->get_parallel_computing() == obj2->get_parallel_computing());
  
  /*------------------------------------------------------------------ simulation schemes */
  
  assert(obj1->get_energy_costs_scheme() == obj2->get_energy_costs_scheme());
  assert(obj1->get_membrane_permeability_scheme() == obj2->get_membrane_permeability_scheme());
  assert(obj1->get_metabolic_inheritance_scheme() == obj2->get_metabolic_inheritance_scheme());
  assert(obj1->get_enzymatic_inheritance_scheme() == obj2->get_enzymatic_inheritance_scheme());
  assert(obj1->get_co_enzyme_activity_scheme() == obj2->get_co_enzyme_activity_scheme());
  assert(obj1->get_score_scheme() == obj2->get_score_scheme());
  assert(obj1->get_selection_threshold() == obj2->get_selection_threshold());
  
  /*------------------------------------------------------------------ space */
  
  assert(obj1->get_width() == obj2->get_width());
  assert(obj1->get_height() == obj2->get_height());
  
  /*------------------------------------------------------------------ output */
  
  assert(obj1->get_simulation_backup_step() == obj2->get_simulation_backup_step());
  assert(obj1->get_figures_generation_step() == obj2->get_figures_generation_step());
  
  /*------------------------------------------------------------------ genome */
  
  assert(obj1->get_metabolite_tag_initial_range()->min == obj2->get_metabolite_tag_initial_range()->min);
  assert(obj1->get_metabolite_tag_initial_range()->max == obj2->get_metabolite_tag_initial_range()->max);
  
  assert(obj1->get_binding_site_tag_initial_range()->min == obj2->get_binding_site_tag_initial_range()->min);
  assert(obj1->get_binding_site_tag_initial_range()->max == obj2->get_binding_site_tag_initial_range()->max);
  
  assert(obj1->get_co_enzyme_tag_initial_range()->min == obj2->get_co_enzyme_tag_initial_range()->min);
  assert(obj1->get_co_enzyme_tag_initial_range()->max == obj2->get_co_enzyme_tag_initial_range()->max);
  
  assert(obj1->get_transcription_factor_tag_initial_range()->min == obj2->get_transcription_factor_tag_initial_range()->min);
  assert(obj1->get_transcription_factor_tag_initial_range()->max == obj2->get_transcription_factor_tag_initial_range()->max);
  
  assert(obj1->get_transcription_factor_binding_window() == obj2->get_transcription_factor_binding_window());
  
  assert(obj1->get_initial_number_of_NC_units() == obj2->get_initial_number_of_NC_units());
  assert(obj1->get_initial_number_of_E_units() == obj2->get_initial_number_of_E_units());
  assert(obj1->get_initial_number_of_TF_units() == obj2->get_initial_number_of_TF_units());
  assert(obj1->get_initial_number_of_BS_units() == obj2->get_initial_number_of_BS_units());
  assert(obj1->get_initial_number_of_P_units() == obj2->get_initial_number_of_P_units());
  
  assert(obj1->get_point_mutation_rate() == obj2->get_point_mutation_rate());
  assert(obj1->get_duplication_rate() == obj2->get_duplication_rate());
  assert(obj1->get_deletion_rate() == obj2->get_deletion_rate());
  assert(obj1->get_translocation_rate() == obj2->get_translocation_rate());
  assert(obj1->get_inversion_rate() == obj2->get_inversion_rate());
  assert(obj1->get_transition_rate() == obj2->get_transition_rate());
  assert(obj1->get_breakpoint_rate() == obj2->get_breakpoint_rate());
  
  assert(obj1->get_substrate_tag_mutation_size() == obj2->get_substrate_tag_mutation_size());
  assert(obj1->get_product_tag_mutation_size() == obj2->get_product_tag_mutation_size());
  assert(obj1->get_kcat_mutation_size() == obj2->get_kcat_mutation_size());
  assert(obj1->get_kcat_km_ratio_mutation_size() == obj2->get_kcat_km_ratio_mutation_size());
  assert(obj1->get_binding_site_tag_mutation_size() == obj2->get_binding_site_tag_mutation_size());
  assert(obj1->get_co_enzyme_tag_mutation_size() == obj2->get_co_enzyme_tag_mutation_size());
  assert(obj1->get_transcription_factor_tag_mutation_size() == obj2->get_transcription_factor_tag_mutation_size());
  assert(obj1->get_basal_expression_level_mutation_size() == obj2->get_basal_expression_level_mutation_size());
  
  /*------------------------------------------------------------------ genetic regulation network */
  
  assert(obj1->get_genetic_regulation_network_timestep() == obj2->get_genetic_regulation_network_timestep());
  
  assert(obj1->get_hill_function_theta() == obj2->get_hill_function_theta());
  assert(obj1->get_hill_function_n() == obj2->get_hill_function_n());
  assert(obj1->get_protein_degradation_rate() == obj2->get_protein_degradation_rate());
  
  /*------------------------------------------------------------------ metabolic network */
  
  assert(obj1->get_metabolism_timestep() == obj2->get_metabolism_timestep());
  
  assert(obj1->get_essential_metabolites_toxicity_threshold() == obj2->get_essential_metabolites_toxicity_threshold());
  assert(obj1->get_non_essential_metabolites_toxicity_threshold() == obj2->get_non_essential_metabolites_toxicity_threshold());
  
  assert(obj1->get_initial_metabolites_amount_in_cells() == obj2->get_initial_metabolites_amount_in_cells());
  
  assert(obj1->get_energy_toxicity_threshold() == obj2->get_energy_toxicity_threshold());
  
  /*------------------------------------------------------------------ energy */
  
  assert(obj1->get_energy_transcription_cost() == obj2->get_energy_transcription_cost());
  assert(obj1->get_energy_degradation_cost() == obj2->get_energy_degradation_cost());
  assert(obj1->get_energy_enzymatic_cost() == obj2->get_energy_enzymatic_cost());
  assert(obj1->get_energy_pumping_cost() == obj2->get_energy_pumping_cost());
  
  assert(obj1->get_energy_dissipation_rate() == obj2->get_energy_dissipation_rate());
  
  assert(obj1->get_energy_toxicity_threshold() == obj2->get_energy_toxicity_threshold());
  
  /*------------------------------------------------------------------ cell */
  
  assert(obj1->get_membrane_permeability() == obj2->get_membrane_permeability());
  
  /*------------------------------------------------------------------ population */
  
  assert(obj1->get_death_probability() == obj2->get_death_probability());
  assert(obj1->get_migration_rate() == obj2->get_migration_rate());
  assert(obj1->get_hgt_rate() == obj2->get_hgt_rate());
  
  /*------------------------------------------------------------------ environment */
  
  assert(obj1->get_environment_properties()->number_of_init_cycles == obj2->get_environment_properties()->number_of_init_cycles);
  
  assert(obj1->get_environment_properties()->species_tag_range.min == obj2->get_environment_properties()->species_tag_range.min);
  assert(obj1->get_environment_properties()->species_tag_range.max == obj2->get_environment_properties()->species_tag_range.max);
  
  assert(obj1->get_environment_properties()->concentration_range.min == obj2->get_environment_properties()->concentration_range.min);
  assert(obj1->get_environment_properties()->concentration_range.max == obj2->get_environment_properties()->concentration_range.max);
  
  assert(obj1->get_environment_properties()->number_of_species_range.min == obj2->get_environment_properties()->number_of_species_range.min);
  assert(obj1->get_environment_properties()->number_of_species_range.max == obj2->get_environment_properties()->number_of_species_range.max);
  
  assert(obj1->get_environment_properties()->interaction_scheme == obj2->get_environment_properties()->interaction_scheme);
  assert(obj1->get_environment_properties()->renewal_scheme == obj2->get_environment_properties()->renewal_scheme);
  assert(obj1->get_environment_properties()->variation_scheme == obj2->get_environment_properties()->variation_scheme);
  assert(obj1->get_environment_properties()->localization_scheme == obj2->get_environment_properties()->localization_scheme);
  assert(obj1->get_environment_properties()->metabolic_scheme == obj2->get_environment_properties()->metabolic_scheme);
  
  assert(obj1->get_environment_properties()->introduction_rate == obj2->get_environment_properties()->introduction_rate);
  assert(obj1->get_environment_properties()->diffusion_coefficient == obj2->get_environment_properties()->diffusion_coefficient);
  assert(obj1->get_environment_properties()->degradation_rate == obj2->get_environment_properties()->degradation_rate);
  
  /*------------------------------------------------------------------ prime numbers */
  
  for (size_t i = 0; i < 5000; i++)
  {
    assert(obj1->get_prime_numbers()[i] == obj2->get_prime_numbers()[i]);
  }
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test genetic unit struct equality
 * \details  --
 * \param    genetic_unit* obj1
 * \param    genetic_unit* obj2
 * \return   \e void
 */
void IntegratedTests::genetic_unit_isEqualTo( genetic_unit* obj1, genetic_unit* obj2 )
{
  /*------------------------------------------------------------------ Global attributes */
  
  assert(obj1->type == obj2->type);
  assert(obj1->identifier == obj2->identifier);
  assert(obj1->parent_identifier == obj2->parent_identifier);
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  assert(obj1->s == obj2->s);
  assert(obj1->p == obj2->p);
  assert(obj1->kcat == obj2->kcat);
  assert(obj1->kcat_km_ratio == obj2->kcat_km_ratio);
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  assert(obj1->BS_tag == obj2->BS_tag);
  assert(obj1->coE_tag == obj2->coE_tag);
  assert(obj1->free_activity == obj2->free_activity);
  assert(obj1->bound_activity == obj2->bound_activity);
  assert(obj1->binding_window == obj2->binding_window);
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  assert(obj1->TF_tag == obj2->TF_tag);
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  assert(obj1->basal_expression_level == obj2->basal_expression_level);
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  assert(obj1->functional == obj2->functional);
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test MutationVector equality
 * \details  --
 * \param    MutationVector* obj1
 * \param    MutationVector* obj2
 * \return   \e void
 */
void IntegratedTests::MutationVector_isEqualTo( MutationVector* obj1, MutationVector* obj2 )
{
  genetic_unit_isEqualTo(obj1->get_dX(), obj2->get_dX());
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test MutationEvent equality
 * \details  --
 * \param    MutationEvent* obj1
 * \param    MutationEvent* obj2
 * \return   \e void
 */
void IntegratedTests::MutationEvent_isEqualTo( MutationEvent* obj1, MutationEvent* obj2 )
{
  assert(obj1->get_mutation_type() == obj2->get_mutation_type());
  assert(obj1->get_point_mutation_location() == obj2->get_point_mutation_location());
  assert(obj1->get_hgt_insert() == obj2->get_hgt_insert());
  assert(obj1->get_nb_NC() == obj2->get_nb_NC());
  assert(obj1->get_nb_E() == obj2->get_nb_E());
  assert(obj1->get_nb_TF() == obj2->get_nb_TF());
  assert(obj1->get_nb_BS() == obj2->get_nb_BS());
  assert(obj1->get_nb_P() == obj2->get_nb_P());
  assert(obj1->get_src_breakpoint1() == obj2->get_src_breakpoint1());
  assert(obj1->get_src_breakpoint2() == obj2->get_src_breakpoint2());
  assert(obj1->get_tgt_breakpoint() == obj2->get_tgt_breakpoint());
  assert(obj1->get_size() == obj2->get_size());
  if (obj1->get_mutation_type() != POINT_MUTATION)
  {
    assert(obj1->get_mutation_vector() == NULL && obj2->get_mutation_vector() == NULL);
  }
  else
  {
    MutationVector_isEqualTo(obj1->get_mutation_vector(), obj2->get_mutation_vector());
  }
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test ReplicationReport equality
 * \details  --
 * \param    ReplicationReport* obj1
 * \param    ReplicationReport* obj2
 * \return   \e void
 */
void IntegratedTests::ReplicationReport_isEqualTo( ReplicationReport* obj1, ReplicationReport* obj2 )
{
  /*------------------------------------------------------------------ Genome structure */
  
  assert(obj1->get_old_genome_size() == obj2->get_old_genome_size());
  assert(obj1->get_new_genome_size() == obj2->get_new_genome_size());
  assert(obj1->get_genome_functional_size() == obj2->get_genome_functional_size());
  assert(obj1->get_genome_nb_NC() == obj2->get_genome_nb_NC());
  assert(obj1->get_genome_nb_E() == obj2->get_genome_nb_E());
  assert(obj1->get_genome_nb_TF() == obj2->get_genome_nb_TF());
  assert(obj1->get_genome_nb_BS() == obj2->get_genome_nb_BS());
  assert(obj1->get_genome_nb_P() == obj2->get_genome_nb_P());
  assert(obj1->get_genome_nb_inner_enzymes() == obj2->get_genome_nb_inner_enzymes());
  assert(obj1->get_genome_nb_inflow_pumps() == obj2->get_genome_nb_inflow_pumps());
  assert(obj1->get_genome_nb_outflow_pumps() == obj2->get_genome_nb_outflow_pumps());
  
  /*------------------------------------------------------------------ GRN data */
  
  assert(obj1->get_nb_functional_regions() == obj2->get_nb_functional_regions());
  assert(obj1->get_nb_enhancers() == obj2->get_nb_enhancers());
  assert(obj1->get_nb_operators() == obj2->get_nb_operators());
  assert(obj1->get_nb_E_regions() == obj2->get_nb_E_regions());
  assert(obj1->get_nb_TF_regions() == obj2->get_nb_TF_regions());
  assert(obj1->get_nb_mixed_regions() == obj2->get_nb_mixed_regions());
  assert(obj1->get_mean_functional_region_size() == obj2->get_mean_functional_region_size());
  assert(obj1->get_mean_E_region_size() == obj2->get_mean_E_region_size());
  assert(obj1->get_mean_TF_region_size() == obj2->get_mean_TF_region_size());
  assert(obj1->get_mean_mixed_region_size() == obj2->get_mean_mixed_region_size());
  assert(obj1->get_mean_enhancer_size() == obj2->get_mean_enhancer_size());
  assert(obj1->get_mean_operator_size() == obj2->get_mean_operator_size());
  assert(obj1->get_mean_operon_size() == obj2->get_mean_operon_size());
  assert(obj1->get_mean_E_operon_size() == obj2->get_mean_E_operon_size());
  assert(obj1->get_mean_TF_operon_size() == obj2->get_mean_TF_operon_size());
  assert(obj1->get_mean_mixed_operon_size() == obj2->get_mean_mixed_operon_size());
  
  /*------------------------------------------------------------------ Genetic redundancy */
  
  assert(obj1->get_mean_regulation_redundancy() == obj2->get_mean_regulation_redundancy());
  assert(obj1->get_mean_metabolic_redundancy() == obj2->get_mean_metabolic_redundancy());
  
  /*------------------------------------------------------------------ Inherited proteins structure */
  
  assert(obj1->get_inherited_size() == obj2->get_inherited_size());
  assert(obj1->get_inherited_nb_E() == obj2->get_inherited_nb_E());
  assert(obj1->get_inherited_nb_TF() == obj2->get_inherited_nb_TF());
  assert(obj1->get_inherited_nb_inner_enzymes() == obj2->get_inherited_nb_inner_enzymes());
  assert(obj1->get_inherited_nb_inflow_pumps() == obj2->get_inherited_nb_inflow_pumps());
  assert(obj1->get_inherited_nb_outflow_pumps() == obj2->get_inherited_nb_outflow_pumps());
  
  /*------------------------------------------------------------------ Phenotype */
  
  assert(obj1->get_id() == obj2->get_id());
  assert(obj1->get_parent_id() == obj2->get_parent_id());
  
  assert(obj1->get_generation() == obj2->get_generation());
  assert(obj1->get_x() == obj2->get_x());
  assert(obj1->get_y() == obj2->get_y());
  assert(obj1->get_number_of_updates() == obj2->get_number_of_updates());
  assert(obj1->get_number_of_divisions() == obj2->get_number_of_divisions());
  assert(obj1->get_birth_time() == obj2->get_birth_time());
  assert(obj1->get_death_time() == obj2->get_death_time());
  assert(obj1->get_lifespan() == obj2->get_lifespan());
  assert(obj1->get_toxicity() == obj2->get_toxicity());
  assert(obj1->get_inherited_TF_amount() == obj2->get_inherited_TF_amount());
  assert(obj1->get_inherited_E_amount() == obj2->get_inherited_E_amount());
  assert(obj1->get_TF_amount() == obj2->get_TF_amount());
  assert(obj1->get_E_amount() == obj2->get_E_amount());
  assert(obj1->get_inherited_metabolic_amount() == obj2->get_inherited_metabolic_amount());
  assert(obj1->get_min_metabolic_amount() == obj2->get_min_metabolic_amount());
  assert(obj1->get_metabolic_amount() == obj2->get_metabolic_amount());
  assert(obj1->get_max_metabolic_amount() == obj2->get_max_metabolic_amount());
  assert(obj1->get_metabolic_uptake() == obj2->get_metabolic_uptake());
  assert(obj1->get_metabolic_release() == obj2->get_metabolic_release());
  assert(obj1->get_min_energy() == obj2->get_min_energy());
  assert(obj1->get_mean_energy() == obj2->get_mean_energy());
  assert(obj1->get_max_energy() == obj2->get_max_energy());
  assert(obj1->get_min_score() == obj2->get_min_score());
  assert(obj1->get_mean_score() == obj2->get_mean_score());
  assert(obj1->get_max_score() == obj2->get_max_score());
  assert(obj1->get_metabolic_growth_rate() == obj2->get_metabolic_growth_rate());
  assert(obj1->get_Dmetabolic_growth_rate() == obj2->get_Dmetabolic_growth_rate());
  assert(obj1->get_grn_nb_nodes() == obj2->get_grn_nb_nodes());
  assert(obj1->get_grn_nb_edges() == obj2->get_grn_nb_edges());
  assert(obj1->get_metabolic_nb_nodes() == obj2->get_metabolic_nb_nodes());
  assert(obj1->get_metabolic_nb_edges() == obj2->get_metabolic_nb_edges());
  
  assert(obj1->get_trophic_group() == obj2->get_trophic_group());
  assert(obj1->get_trophic_level() == obj2->get_trophic_level());
  
  /*------------------------------------------------------------------ List of mutation events */
  
  assert(obj1->get_number_of_events() == obj2->get_number_of_events());
  assert(obj1->get_list_of_events()->size() == obj2->get_list_of_events()->size());
  for (size_t i = 0; i < obj1->get_list_of_events()->size(); i++)
  {
    MutationEvent_isEqualTo(obj1->get_list_of_events()->at(i), obj2->get_list_of_events()->at(i));
  }
  
  /*------------------------------------------------------------------ Point mutations data */
  
  assert(obj1->get_nb_point_mutations() == obj2->get_nb_point_mutations());
  
  assert(obj1->get_nb_NC_point_mutations() == obj2->get_nb_NC_point_mutations());
  assert(obj1->get_nb_E_point_mutations() == obj2->get_nb_E_point_mutations());
  assert(obj1->get_nb_TF_point_mutations() == obj2->get_nb_TF_point_mutations());
  assert(obj1->get_nb_BS_point_mutations() == obj2->get_nb_BS_point_mutations());
  assert(obj1->get_nb_P_point_mutations() == obj2->get_nb_P_point_mutations());
  
  assert(obj1->get_nb_NC_to_E_transitions() == obj2->get_nb_NC_to_E_transitions());
  assert(obj1->get_nb_NC_to_TF_transitions() == obj2->get_nb_NC_to_TF_transitions());
  assert(obj1->get_nb_NC_to_BS_transitions() == obj2->get_nb_NC_to_BS_transitions());
  assert(obj1->get_nb_NC_to_P_transitions() == obj2->get_nb_NC_to_P_transitions());
  
  assert(obj1->get_nb_E_to_NC_transitions() == obj2->get_nb_E_to_NC_transitions());
  assert(obj1->get_nb_E_to_TF_transitions() == obj2->get_nb_E_to_TF_transitions());
  assert(obj1->get_nb_E_to_BS_transitions() == obj2->get_nb_E_to_BS_transitions());
  assert(obj1->get_nb_E_to_P_transitions() == obj2->get_nb_E_to_P_transitions());
  
  assert(obj1->get_nb_TF_to_NC_transitions() == obj2->get_nb_TF_to_NC_transitions());
  assert(obj1->get_nb_TF_to_E_transitions() == obj2->get_nb_TF_to_E_transitions());
  assert(obj1->get_nb_TF_to_BS_transitions() == obj2->get_nb_TF_to_BS_transitions());
  assert(obj1->get_nb_TF_to_P_transitions() == obj2->get_nb_TF_to_P_transitions());
  
  assert(obj1->get_nb_BS_to_NC_transitions() == obj2->get_nb_BS_to_NC_transitions());
  assert(obj1->get_nb_BS_to_E_transitions() == obj2->get_nb_BS_to_E_transitions());
  assert(obj1->get_nb_BS_to_TF_transitions() == obj2->get_nb_BS_to_TF_transitions());
  assert(obj1->get_nb_BS_to_P_transitions() == obj2->get_nb_BS_to_P_transitions());
  
  assert(obj1->get_nb_P_to_NC_transitions() == obj2->get_nb_P_to_NC_transitions());
  assert(obj1->get_nb_P_to_E_transitions() == obj2->get_nb_P_to_E_transitions());
  assert(obj1->get_nb_P_to_TF_transitions() == obj2->get_nb_P_to_TF_transitions());
  assert(obj1->get_nb_P_to_BS_transitions() == obj2->get_nb_P_to_BS_transitions());
  
  assert(obj1->get_mean_s_mutation_size() == obj2->get_mean_s_mutation_size());
  assert(obj1->get_mean_p_mutation_size() == obj2->get_mean_p_mutation_size());
  assert(obj1->get_mean_kcat_mutation_size() == obj2->get_mean_kcat_mutation_size());
  assert(obj1->get_mean_kcat_km_ratio_mutation_size() == obj2->get_mean_kcat_km_ratio_mutation_size());
  
  assert(obj1->get_mean_BS_tag_mutation_size() == obj2->get_mean_BS_tag_mutation_size());
  assert(obj1->get_mean_coE_tag_mutation_size() == obj2->get_mean_coE_tag_mutation_size());
  
  assert(obj1->get_mean_TF_tag_mutation_size() == obj2->get_mean_TF_tag_mutation_size());
  
  assert(obj1->get_mean_basal_expression_level_mutation_size() == obj2->get_mean_basal_expression_level_mutation_size());
  
  /*------------------------------------------------------------------ HGT data */
  
  assert(obj1->get_nb_HGT() == obj2->get_nb_HGT());
  assert(obj1->get_mean_HGT_size() == obj2->get_mean_HGT_size());
  assert(obj1->get_nb_NC_HGT() == obj2->get_nb_NC_HGT());
  assert(obj1->get_nb_E_HGT() == obj2->get_nb_E_HGT());
  assert(obj1->get_nb_TF_HGT() == obj2->get_nb_TF_HGT());
  assert(obj1->get_nb_BS_HGT() == obj2->get_nb_BS_HGT());
  assert(obj1->get_nb_P_HGT() == obj2->get_nb_P_HGT());
  
  /*------------------------------------------------------------------ Rearrangements data */
  
  assert(obj1->get_nb_rearrangements() == obj2->get_nb_rearrangements());
  
  assert(obj1->get_nb_duplicated_NC() == obj2->get_nb_duplicated_NC());
  assert(obj1->get_nb_duplicated_E() == obj2->get_nb_duplicated_E());
  assert(obj1->get_nb_duplicated_TF() == obj2->get_nb_duplicated_TF());
  assert(obj1->get_nb_duplicated_BS() == obj2->get_nb_duplicated_BS());
  assert(obj1->get_nb_duplicated_P() == obj2->get_nb_duplicated_P());
  
  assert(obj1->get_nb_deleted_NC() == obj2->get_nb_deleted_NC());
  assert(obj1->get_nb_deleted_E() == obj2->get_nb_deleted_E());
  assert(obj1->get_nb_deleted_TF() == obj2->get_nb_deleted_TF());
  assert(obj1->get_nb_deleted_BS() == obj2->get_nb_deleted_BS());
  assert(obj1->get_nb_deleted_P() == obj2->get_nb_deleted_P());
  
  assert(obj1->get_nb_duplications() == obj2->get_nb_duplications());
  assert(obj1->get_nb_deletions() == obj2->get_nb_deletions());
  assert(obj1->get_nb_translocations() == obj2->get_nb_translocations());
  assert(obj1->get_nb_inversions() == obj2->get_nb_inversions());
  
  assert(obj1->get_mean_rearrangement_size() == obj2->get_mean_rearrangement_size());
  assert(obj1->get_mean_duplication_size() == obj2->get_mean_duplication_size());
  assert(obj1->get_mean_deletion_size() == obj2->get_mean_deletion_size());
  assert(obj1->get_mean_translocation_size() == obj2->get_mean_translocation_size());
  assert(obj1->get_mean_inversion_size() == obj2->get_mean_inversion_size());
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test Genome equality
 * \details  --
 * \param    Genome* obj1
 * \param    Genome* obj2
 * \return   \e void
 */
void IntegratedTests::Genome_isEqualTo( Genome* obj1, Genome* obj2 )
{
  /*------------------------------------------------------------------ genetic sequence */
  
  assert(obj1->get_size() == obj2->get_size());
  for (size_t i = 0; i < obj1->get_size(); i++)
  {
    genetic_unit_isEqualTo(obj1->get_genetic_unit(i), obj2->get_genetic_unit(i));
  }
  assert(obj1->get_buffer_size() == obj2->get_buffer_size());
  assert(obj1->get_coding_size() == obj2->get_coding_size());
  assert(obj1->get_non_coding_size() == obj2->get_non_coding_size());
  
  /*------------------------------------------------------------------ concentration vector */
  
  for (size_t i = 0; i < obj1->get_size(); i++)
  {
    double conc1 = obj1->get_concentration_vector()[i];
    double conc2 = obj2->get_concentration_vector()[i];
    assert(conc1 == conc2);
    assert(conc1 >= 0.0);
    
    (void)conc1;
    (void)conc2;
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  assert(obj1->get_nb_NC() == obj2->get_nb_NC());
  assert(obj1->get_nb_E() == obj2->get_nb_E());
  assert(obj1->get_nb_TF() == obj2->get_nb_TF());
  assert(obj1->get_nb_BS() == obj2->get_nb_BS());
  assert(obj1->get_nb_P() == obj2->get_nb_P());
  assert(obj1->get_nb_inner_enzymes() == obj2->get_nb_inner_enzymes());
  assert(obj1->get_nb_inflow_pumps() == obj2->get_nb_inflow_pumps());
  assert(obj1->get_nb_outflow_pumps() == obj2->get_nb_outflow_pumps());
  
  /*------------------------------------------------------------------ genetic unit positions */
  
  for (size_t i = 0; i < obj1->get_nb_TF(); i++)
  {
    assert(obj1->get_TFi()[i] == obj2->get_TFi()[i]);
    assert(obj1->get_TFi()[i] < obj1->get_size());
    assert(obj1->get_genetic_unit(obj1->get_TFi()[i])->type == TRANSCRIPTION_FACTOR);
    assert(obj2->get_genetic_unit(obj2->get_TFi()[i])->type == TRANSCRIPTION_FACTOR);
  }
  for (size_t i = 0; i < obj1->get_nb_P(); i++)
  {
    assert(obj1->get_Pi()[i] == obj2->get_Pi()[i]);
    assert(obj1->get_Pi()[i] < obj1->get_size());
    assert(obj1->get_genetic_unit(obj1->get_Pi()[i])->type == PROMOTER);
    assert(obj2->get_genetic_unit(obj2->get_Pi()[i])->type == PROMOTER);
  }
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test InheritedProteins equality
 * \details  --
 * \param    InheritedProteins* obj1
 * \param    InheritedProteins* obj2
 * \return   \e void
 */
void IntegratedTests::InheritedProteins_isEqualTo( InheritedProteins* obj1, InheritedProteins* obj2 )
{
  /*------------------------------------------------------------------ genetic sequence */
  
  assert(obj1->get_size() == obj2->get_size());
  for (size_t i = 0; i < obj1->get_size(); i++)
  {
    genetic_unit_isEqualTo(obj1->get_genetic_unit(i), obj2->get_genetic_unit(i));
  }
  assert(obj1->get_buffer_size() == obj2->get_buffer_size());
  
  /*------------------------------------------------------------------ concentration vector */
  
  for (size_t i = 0; i < obj1->get_size(); i++)
  {
    assert(obj1->get_concentration_vector()[i] == obj2->get_concentration_vector()[i]);
    assert(obj1->get_concentration_vector()[i] >= 0.0);
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  assert(obj1->get_nb_E() == obj2->get_nb_E());
  assert(obj1->get_nb_TF() == obj2->get_nb_TF());
  assert(obj1->get_nb_inner_enzymes() == obj2->get_nb_inner_enzymes());
  assert(obj1->get_nb_inflow_pumps() == obj2->get_nb_inflow_pumps());
  assert(obj1->get_nb_outflow_pumps() == obj2->get_nb_outflow_pumps());
  
  /*------------------------------------------------------------------ genetic unit positions */
  
  for (size_t i = 0; i < obj1->get_nb_E(); i++)
  {
    assert(obj1->get_Ei()[i] == obj2->get_Ei()[i]);
    assert(obj1->get_Ei()[i] < obj1->get_size());
    assert(obj1->get_genetic_unit(obj1->get_Ei()[i])->type == ENZYME);
    assert(obj2->get_genetic_unit(obj2->get_Ei()[i])->type == ENZYME);
  }
  for (size_t i = 0; i < obj1->get_nb_TF(); i++)
  {
    assert(obj1->get_TFi()[i] == obj2->get_TFi()[i]);
    assert(obj1->get_TFi()[i] < obj1->get_size());
    assert(obj1->get_genetic_unit(obj1->get_TFi()[i])->type == TRANSCRIPTION_FACTOR);
    assert(obj2->get_genetic_unit(obj2->get_TFi()[i])->type == TRANSCRIPTION_FACTOR);
  }
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test reaction_list equality
 * \details  --
 * \param    reaction_list* obj1
 * \param    reaction_list* obj2
 * \return   \e void
 */
void IntegratedTests::reaction_list_isEqualTo( reaction_list* obj1, reaction_list* obj2 )
{
  /*------------------------------------------------------------------ list of genetic regulation network equations */
  
  assert(obj1->grn_N == obj2->grn_N);
  
  /* For each promoter */
  size_t enhancer_TF_index     = 0;
  size_t operator_TF_index     = 0;
  size_t regulated_genes_index = 0;
  for (size_t i = 0; i < obj1->grn_N; i++)
  {
    assert(obj1->grn_promoter[i] == obj2->grn_promoter[i]);
    assert(obj1->grn_Nenhancer[i] == obj2->grn_Nenhancer[i]);
    assert(obj1->grn_Noperator[i] == obj2->grn_Noperator[i]);
    assert(obj1->grn_Ngenes[i] == obj2->grn_Ngenes[i]);
    assert(obj1->grn_enhancer_nb_BS[i] == obj2->grn_enhancer_nb_BS[i]);
    assert(obj1->grn_operator_nb_BS[i] == obj2->grn_operator_nb_BS[i]);
    
    /* For each TF linked to the enhancer */
    for (size_t j = 0; j < obj1->grn_Nenhancer[i]; j++)
    {
      assert(obj1->grn_enhancer_TF_list[enhancer_TF_index] == obj2->grn_enhancer_TF_list[enhancer_TF_index]);
      assert(obj1->grn_enhancer_affinity_list[enhancer_TF_index] == obj2->grn_enhancer_affinity_list[enhancer_TF_index]);
      assert(obj1->grn_enhancer_coe_list[enhancer_TF_index] == obj2->grn_enhancer_coe_list[enhancer_TF_index]);
      assert(obj1->grn_enhancer_coe_type[enhancer_TF_index] == obj2->grn_enhancer_coe_type[enhancer_TF_index]);
      enhancer_TF_index++;
    }
    
    /* For each TF linked to the operator */
    for (size_t j = 0; j < obj1->grn_Noperator[i]; j++)
    {
      assert(obj1->grn_operator_TF_list[operator_TF_index] == obj2->grn_operator_TF_list[operator_TF_index]);
      assert(obj1->grn_operator_affinity_list[operator_TF_index] == obj2->grn_operator_affinity_list[operator_TF_index]);
      assert(obj1->grn_operator_coe_list[operator_TF_index] == obj2->grn_operator_coe_list[operator_TF_index]);
      assert(obj1->grn_operator_coe_type[operator_TF_index] == obj2->grn_operator_coe_type[operator_TF_index]);
      operator_TF_index++;
    }
    /* For each regulated gene */
    for (size_t j = 0; j < obj1->grn_Ngenes[i]; j++)
    {
      assert(obj1->grn_regulated_genes[regulated_genes_index] == obj2->grn_regulated_genes[regulated_genes_index]);
      regulated_genes_index++;
    }
    assert(obj1->grn_beta[i] == obj2->grn_beta[i]);
  }
  assert(obj1->grn_hill_n == obj2->grn_hill_n);
  assert(obj1->grn_theta_n == obj2->grn_theta_n);
  assert(obj1->grn_protein_degradation_rate == obj2->grn_protein_degradation_rate);
  assert(obj1->grn_energy_transcription_cost == obj2->grn_energy_transcription_cost);
  assert(obj1->grn_energy_degradation_cost == obj2->grn_energy_degradation_cost);
  assert(obj1->grn_energy_dissipation_rate == obj2->grn_energy_dissipation_rate);
  
  /*------------------------------------------------------------------ list of metabolic equations */
  
  assert(obj1->metabolic_N == obj2->metabolic_N);
  for (size_t i = 0; i < obj1->metabolic_N; i++)
  {
    assert(obj1->metabolic_type[i] == obj2->metabolic_type[i]);
    assert(obj1->metabolic_s[i] == obj2->metabolic_s[i]);
    assert(obj1->metabolic_p[i] == obj2->metabolic_p[i]);
    assert(obj1->metabolic_km[i] == obj2->metabolic_km[i]);
    assert(obj1->metabolic_kcat[i] == obj2->metabolic_kcat[i]);
    assert(obj1->metabolic_delta_g[i] == obj2->metabolic_delta_g[i]);
    assert(obj1->metabolic_e[i] == obj2->metabolic_e[i]);
  }
  
  /*------------------------------------------------------------------ size of each state vector subspace */
  
  assert(obj1->Ngenome == obj2->Ngenome);
  assert(obj1->Ninherited == obj2->Ninherited);
  assert(obj1->Ncell == obj2->Ncell);
  assert(obj1->Nenv == obj2->Nenv);
  assert(obj1->N == obj2->N);
  
  /*------------------------------------------------------------------ modeling schemes */
  
  assert(obj1->energy_constraints == obj2->energy_constraints);
  assert(obj1->membrane_permeability == obj2->membrane_permeability);
  assert(obj1->metabolic_inheritance == obj2->metabolic_inheritance);
  assert(obj1->enzymatic_inheritance == obj2->enzymatic_inheritance);
  assert(obj1->co_enzyme_activity == obj2->co_enzyme_activity);
  
  /*------------------------------------------------------------------ membrane permeability */
  
  assert(obj1->kmp == obj2->kmp);
  
  /*------------------------------------------------------------------ metabolic exchanges with the environment */
  
  assert(obj1->metabolic_uptake == obj2->metabolic_uptake);
  assert(obj1->metabolic_release == obj2->metabolic_release);
  
  /*------------------------------------------------------------------ timestep ratio */
  
  assert(obj1->timestep_ratio == obj2->timestep_ratio);
  
  /*------------------------------------------------------------------ environment interaction scheme */
  
  assert(obj1->interacts_with_environment == obj2->interacts_with_environment);
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test ODE equality
 * \details  --w
 * \param    ODE* obj1
 * \param    ODE* obj2
 * \return   \e void
 */
void IntegratedTests::ODE_isEqualTo( ODE* obj1, ODE* obj2 )
{
  if (obj1 == NULL)
  {
    assert(obj2 == NULL);
  }
  else
  {
    reaction_list_isEqualTo(obj1->get_reaction_list(), obj2->get_reaction_list());
  }
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test SpeciesList equality
 * \details  --
 * \param    SpeciesList* obj1
 * \param    SpeciesList* obj2
 * \return   \e void
 */
void IntegratedTests::SpeciesList_isEqualTo( SpeciesList* obj1, SpeciesList* obj2 )
{
  assert(obj1->get_size() == obj2->get_size());
  assert(obj1->get_amount() == obj2->get_amount());
  for (size_t i = 0; i < obj1->get_size(); i++)
  {
    assert(obj1->get_X()[i] == obj2->get_X()[i]);
    assert(obj1->get((int)i+1) == obj2->get((int)i+1));
    assert(obj1->get_X()[i] == obj2->get((int)i+1));
    assert(obj1->get((int)i+1) == obj2->get_X()[i]);
    assert(obj1->get_X()[i] >= 0.0);
  }
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test Cell equality
 * \details  --
 * \param    Cell* obj1
 * \param    Cell* obj2
 * \return   \e void
 */
void IntegratedTests::Cell_isEqualTo( Cell* obj1, Cell* obj2 )
{
  /*------------------------------------------------------------------ main cell classes */
  
  ReplicationReport_isEqualTo(obj1->get_replication_report(), obj2->get_replication_report());
  Genome_isEqualTo(obj1->get_genome(), obj2->get_genome());
  if (obj1->get_inherited_proteins() != NULL)
  {
    InheritedProteins_isEqualTo(obj1->get_inherited_proteins(), obj2->get_inherited_proteins());
  }
  SpeciesList_isEqualTo(obj1->get_inherited_species_list(), obj2->get_inherited_species_list());
  SpeciesList_isEqualTo(obj1->get_species_list(), obj2->get_species_list());
  
  assert(obj1->get_inherited_species_list() == obj1->get_species_list());
  
  /*------------------------------------------------------------------ ODE system */
  
  ODE_isEqualTo(obj1->get_ode(), obj2->get_ode());
  
  /*------------------------------------------------------------------ main cell variables */
  
  assert(obj1->get_id() == obj2->get_id());
  assert(obj1->get_parent_id() == obj2->get_parent_id());
  assert(obj1->get_generation() == obj2->get_generation());
  assert(obj1->get_energy() == obj2->get_energy());
  assert(obj1->get_active() == obj2->get_active());
  assert(obj1->get_alive() == obj2->get_alive());
  assert(obj1->get_x() == obj2->get_x());
  assert(obj1->get_y() == obj2->get_y());
  assert(obj1->get_score() == obj2->get_score());
  assert(obj1->get_number_of_updates() == obj2->get_number_of_updates());
  assert(obj1->get_number_of_divisions() == obj2->get_number_of_divisions());
  assert(obj1->isActive() == obj2->isActive());
  assert(obj1->isAlive() == obj2->isAlive());
  assert(obj1->isTagged() == obj2->isTagged());
  
  /*------------------------------------------------------------------ mutation rates */
  
  assert(obj1->get_point_mutation_rate() == obj2->get_point_mutation_rate());
  assert(obj1->get_duplication_rate() == obj2->get_duplication_rate());
  assert(obj1->get_deletion_rate() == obj2->get_deletion_rate());
  assert(obj1->get_translocation_rate() == obj2->get_translocation_rate());
  assert(obj1->get_inversion_rate() == obj2->get_inversion_rate());
  assert(obj1->get_transition_rate() == obj2->get_transition_rate());
  assert(obj1->get_breakpoint_rate() == obj2->get_breakpoint_rate());
  assert(obj1->get_substrate_tag_mutation_size() == obj2->get_substrate_tag_mutation_size());
  assert(obj1->get_product_tag_mutation_size() == obj2->get_product_tag_mutation_size());
  assert(obj1->get_kcat_mutation_size() == obj2->get_kcat_mutation_size());
  assert(obj1->get_kcat_km_ratio_mutation_size() == obj2->get_kcat_km_ratio_mutation_size());
  assert(obj1->get_binding_site_tag_mutation_size() == obj2->get_binding_site_tag_mutation_size());
  assert(obj1->get_co_enzyme_tag_mutation_size() == obj2->get_co_enzyme_tag_mutation_size());
  assert(obj1->get_transcription_factor_tag_mutation_size() == obj2->get_transcription_factor_tag_mutation_size());
  assert(obj1->get_basal_expression_level_mutation_size() == obj2->get_basal_expression_level_mutation_size());
  
  /*------------------------------------------------------------------ time variables */
  
  assert(obj1->get_birth_time() == obj2->get_birth_time());
  assert(obj1->get_death_time() == obj2->get_death_time());
  assert(obj1->get_lifespan() == obj2->get_lifespan());
  
  /*------------------------------------------------------------------ global phenotypic variables */
  
  assert(obj1->get_toxicity() == obj2->get_toxicity());
  assert(obj1->get_inherited_TF_amount() == obj2->get_inherited_TF_amount());
  assert(obj1->get_inherited_E_amount() == obj2->get_inherited_E_amount());
  assert(obj1->get_TF_amount() == obj2->get_TF_amount());
  assert(obj1->get_E_amount() == obj2->get_E_amount());
  assert(obj1->get_inherited_metabolic_amount() == obj2->get_inherited_metabolic_amount());
  assert(obj1->get_min_metabolic_amount() == obj2->get_min_metabolic_amount());
  assert(obj1->get_max_metabolic_amount() == obj2->get_max_metabolic_amount());
  assert(obj1->get_metabolic_uptake() == obj2->get_metabolic_uptake());
  assert(obj1->get_metabolic_release() == obj2->get_metabolic_release());
  assert(obj1->get_min_energy() == obj2->get_min_energy());
  assert(obj1->get_mean_energy() == obj2->get_mean_energy());
  assert(obj1->get_max_energy() == obj2->get_max_energy());
  assert(obj1->get_min_score() == obj2->get_min_score());
  assert(obj1->get_mean_score() == obj2->get_mean_score());
  assert(obj1->get_max_score() == obj2->get_max_score());
  assert(obj1->get_metabolic_growth_rate() == obj2->get_metabolic_growth_rate());
  assert(obj1->get_Dmetabolic_growth_rate() == obj2->get_Dmetabolic_growth_rate());
  assert(obj1->get_grn_nb_nodes() == obj2->get_grn_nb_nodes());
  assert(obj1->get_grn_nb_edges() == obj2->get_grn_nb_edges());
  assert(obj1->get_metabolic_nb_nodes() == obj2->get_metabolic_nb_nodes());
  assert(obj1->get_metabolic_nb_edges() == obj2->get_metabolic_nb_edges());
  assert(obj1->get_inflowing_pumps()->size() == obj2->get_inflowing_pumps()->size());
  for (size_t i = 0; i < obj1->get_inflowing_pumps()->size(); i++)
  {
    assert(obj1->get_inflowing_pumps()->at(i) == obj2->get_inflowing_pumps()->at(i));
  }
  assert(obj1->get_outflowing_pumps()->size() == obj2->get_outflowing_pumps()->size());
  for (size_t i = 0; i < obj1->get_outflowing_pumps()->size(); i++)
  {
    assert(obj1->get_outflowing_pumps()->at(i) == obj2->get_outflowing_pumps()->at(i));
  }
  assert(obj1->get_trophic_group() == obj2->get_trophic_group());
  assert(obj1->get_trophic_level() == obj2->get_trophic_level());
  
  /*------------------------------------------------------------------ cell color */
  
  assert(obj1->get_red_color() == obj2->get_red_color());
  assert(obj1->get_green_color() == obj2->get_green_color());
  assert(obj1->get_blue_color() == obj2->get_blue_color());
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test Population equality
 * \details  --
 * \param    Population* obj1
 * \param    Population* obj2
 * \return   \e void
 */
void IntegratedTests::Population_isEqualTo( Population* obj1, Population* obj2 )
{
  assert(obj1->get_width() == obj2->get_width());
  assert(obj1->get_height() == obj2->get_height());
  for (size_t i = 0; i < obj1->get_width()*obj1->get_height(); i++)
  {
    Cell_isEqualTo(obj1->get_cell(i), obj2->get_cell(i));
  }
  for (size_t i = 0; i < obj1->get_width(); i++)
  {
    for (size_t j = 0; j < obj1->get_height(); j++)
    {
      Cell_isEqualTo(obj1->get_cell(i, j), obj2->get_cell(i, j));
    }
  }
  Cell_isEqualTo(obj1->get_best_cell(), obj2->get_best_cell());
  assert(obj1->get_time() == obj2->get_time());
  assert(obj1->get_population_size() == obj2->get_population_size());
  assert(obj1->get_growth_rate() == obj2->get_growth_rate());
  assert(obj1->get_best_id() == obj2->get_best_id());
  assert(obj1->get_best_position() == obj2->get_best_position());
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test Environment equality
 * \details  --
 * \param    Environment* obj1
 * \param    Environment* obj2
 * \return   \e void
 */
void IntegratedTests::Environment_isEqualTo( Environment* obj1, Environment* obj2 )
{
  assert(obj1->get_width() == obj2->get_width());
  assert(obj1->get_height() == obj2->get_height());
  assert(obj1->get_size() == obj2->get_size());
  assert(obj1->get_total_amount() == obj2->get_total_amount());
  assert(obj1->get_min_amount() == obj2->get_min_amount());
  assert(obj1->get_max_amount() == obj2->get_max_amount());
  assert(obj1->get_inflowing_amount() == obj2->get_inflowing_amount());
  assert(obj1->get_outflowing_amount() == obj2->get_outflowing_amount());
  SpeciesList_isEqualTo(obj1->get_X(), obj2->get_X());
  
  for (size_t i = 0; i < obj1->get_width(); i++)
  {
    for (size_t j = 0; j < obj1->get_height(); j++)
    {
      SpeciesList_isEqualTo(obj1->get_species_list(i, j), obj2->get_species_list(i, j));
      for (int k = 0; k < (int)obj1->get_size(); k++)
      {
        assert(obj1->get(i, j, k+1) == obj2->get(i, j, k+1));
      }
    }
  }
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test TrophicGroup equality
 * \details  --
 * \param    TrophicGroup* obj1
 * \param    TrophicGroup* obj2
 * \return   \e void
 */
void IntegratedTests::TrophicGroup_isEqualTo( TrophicGroup* obj1, TrophicGroup* obj2 )
{
  /*------------------------------------------------------------------ trophic group properties */
  
  assert(obj1->get_identifier() == obj2->get_identifier());
  assert(obj1->_trophic_profile.compare(obj2->_trophic_profile) == 0);
  assert(obj1->_production_profile.compare(obj2->_production_profile) == 0);
  assert(obj1->_uptake_profile.compare(obj2->_uptake_profile) == 0);
  assert(obj1->_release_profile.compare(obj2->_release_profile) == 0);
  assert(obj1->get_necrophagy_links()->size() == obj2->get_necrophagy_links()->size());
  for (size_t i = 0; i < obj1->get_necrophagy_links()->size(); i++)
  {
    assert(obj1->get_necrophagy_links()->at(i) == obj2->get_necrophagy_links()->at(i));
  }
  assert(obj1->get_active_release_links()->size() == obj2->get_active_release_links()->size());
  for (size_t i = 0; i < obj1->get_active_release_links()->size(); i++)
  {
    assert(obj1->get_active_release_links()->at(i) == obj2->get_active_release_links()->at(i));
  }
  assert(obj1->get_appearance_time() == obj2->get_appearance_time());
  assert(obj1->get_lifespan() == obj2->get_lifespan());
  assert(obj1->get_number_of_cells() == obj2->get_number_of_cells());
  assert(obj1->get_trophic_level() == obj2->get_trophic_level());
  
  /*------------------------------------------------------------------ trophic group statistics */
  
  /* PHENOTYPE */
  assert(obj1->get_mean_generations() == obj2->get_mean_generations());
  assert(obj1->get_mean_inherited_TF_amount() == obj2->get_mean_inherited_TF_amount());
  assert(obj1->get_mean_inherited_E_amount() == obj2->get_mean_inherited_E_amount());
  assert(obj1->get_mean_TF_amount() == obj2->get_mean_TF_amount());
  assert(obj1->get_mean_E_amount() == obj2->get_mean_E_amount());
  assert(obj1->get_mean_inherited_metabolic_amount() == obj2->get_mean_inherited_metabolic_amount());
  assert(obj1->get_mean_metabolic_amount() == obj2->get_mean_metabolic_amount());
  assert(obj1->get_mean_energy() == obj2->get_mean_energy());
  assert(obj1->get_mean_score() == obj2->get_mean_score());
  assert(obj1->get_mean_lifespan() == obj2->get_mean_lifespan());
  assert(obj1->get_mean_number_of_divisions() == obj2->get_mean_number_of_divisions());
  assert(obj1->get_mean_toxicity() == obj2->get_mean_toxicity());
  assert(obj1->get_mean_metabolic_uptake() == obj2->get_mean_metabolic_uptake());
  assert(obj1->get_mean_metabolic_release() == obj2->get_mean_metabolic_release());
  assert(obj1->get_mean_metabolic_growth_rate() == obj2->get_mean_metabolic_growth_rate());
  assert(obj1->get_mean_Dmetabolic_growth_rate() == obj2->get_mean_Dmetabolic_growth_rate());
  assert(obj1->get_mean_grn_nb_nodes() == obj2->get_mean_grn_nb_nodes());
  assert(obj1->get_mean_grn_nb_edges() == obj2->get_mean_grn_nb_edges());
  assert(obj1->get_mean_metabolic_nb_nodes() == obj2->get_mean_metabolic_nb_nodes());
  assert(obj1->get_mean_metabolic_nb_edges() == obj2->get_mean_metabolic_nb_edges());
  assert(obj1->get_mean_regulation_redundancy() == obj2->get_mean_regulation_redundancy());
  assert(obj1->get_mean_metabolic_redundancy() == obj2->get_mean_metabolic_redundancy());
  
  /* GENOME STRUCTURE */
  assert(obj1->get_mean_genome_size() == obj2->get_mean_genome_size());
  assert(obj1->get_mean_functional_size() == obj2->get_mean_functional_size());
  assert(obj1->get_mean_genome_nb_NC() == obj2->get_mean_genome_nb_NC());
  assert(obj1->get_mean_genome_nb_E() == obj2->get_mean_genome_nb_E());
  assert(obj1->get_mean_genome_nb_TF() == obj2->get_mean_genome_nb_TF());
  assert(obj1->get_mean_genome_nb_BS() == obj2->get_mean_genome_nb_BS());
  assert(obj1->get_mean_genome_nb_P() == obj2->get_mean_genome_nb_P());
  assert(obj1->get_mean_genome_nb_inner_enzymes() == obj2->get_mean_genome_nb_inner_enzymes());
  assert(obj1->get_mean_genome_nb_inflow_pumps() == obj2->get_mean_genome_nb_inflow_pumps());
  assert(obj1->get_mean_genome_nb_outflow_pumps() == obj2->get_mean_genome_nb_outflow_pumps());
  assert(obj1->get_mean_genome_nb_functional_regions() == obj2->get_mean_genome_nb_functional_regions());
  assert(obj1->get_mean_genome_nb_enhancers() == obj2->get_mean_genome_nb_enhancers());
  assert(obj1->get_mean_genome_nb_operators() == obj2->get_mean_genome_nb_operators());
  assert(obj1->get_mean_genome_nb_E_regions() == obj2->get_mean_genome_nb_E_regions());
  assert(obj1->get_mean_genome_nb_TF_regions() == obj2->get_mean_genome_nb_TF_regions());
  assert(obj1->get_mean_genome_nb_mixed_regions() == obj2->get_mean_genome_nb_mixed_regions());
  assert(obj1->get_mean_genome_functional_region_size() == obj2->get_mean_genome_functional_region_size());
  assert(obj1->get_mean_genome_E_region_size() == obj2->get_mean_genome_E_region_size());
  assert(obj1->get_mean_genome_TF_region_size() == obj2->get_mean_genome_TF_region_size());
  assert(obj1->get_mean_genome_mixed_region_size() == obj2->get_mean_genome_mixed_region_size());
  assert(obj1->get_mean_genome_enhancer_size() == obj2->get_mean_genome_enhancer_size());
  assert(obj1->get_mean_genome_operator_size() == obj2->get_mean_genome_operator_size());
  assert(obj1->get_mean_genome_operon_size() == obj2->get_mean_genome_operon_size());
  assert(obj1->get_mean_genome_E_operon_size() == obj2->get_mean_genome_E_operon_size());
  assert(obj1->get_mean_genome_TF_operon_size() == obj2->get_mean_genome_TF_operon_size());
  assert(obj1->get_mean_genome_mixed_operon_size() == obj2->get_mean_genome_mixed_operon_size());
  
  /* INHERITED STRUCTURE */
  assert(obj1->get_mean_inherited_size() == obj2->get_mean_inherited_size());
  assert(obj1->get_mean_inherited_nb_E() == obj2->get_mean_inherited_nb_E());
  assert(obj1->get_mean_inherited_nb_TF() == obj2->get_mean_inherited_nb_TF());
  assert(obj1->get_mean_inherited_nb_inner_enzymes() == obj2->get_mean_inherited_nb_inner_enzymes());
  assert(obj1->get_mean_inherited_nb_inflow_pumps() == obj2->get_mean_inherited_nb_inflow_pumps());
  assert(obj1->get_mean_inherited_nb_outflow_pumps() == obj2->get_mean_inherited_nb_outflow_pumps());
  
  /*------------------------------------------------------------------ trophic group color */
  
  assert(obj1->get_red_color() == obj2->get_red_color());
  assert(obj1->get_green_color() == obj2->get_green_color());
  assert(obj1->get_blue_color() == obj2->get_blue_color());
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test TrophicNetwork equality
 * \details  --
 * \param    TrophicNetwork* obj1
 * \param    TrophicNetwork* obj2
 * \return   \e void
 */
void IntegratedTests::TrophicNetwork_isEqualTo( TrophicNetwork* obj1, TrophicNetwork* obj2 )
{
  /*------------------------------------------------------------------ Trophic network attributes */
  
  assert(obj1->get_current_id() == obj2->get_current_id());
  assert(obj1->get_number_of_groups() == obj2->get_number_of_groups());
  TrophicGroup* g1 = obj1->get_first_group();
  TrophicGroup* g2 = obj2->get_first_group();
  while (g1 != NULL)
  {
    TrophicGroup_isEqualTo(g1, g2);
    TrophicGroup_isEqualTo(obj1->get_group(g1->get_identifier()), obj2->get_group(g2->get_identifier()));
    g1 = obj1->get_next_group();
    g2 = obj2->get_next_group();
  }
  assert(g1 == NULL);
  assert(g2 == NULL);
  
  /*------------------------------------------------------------------ Trophic network statistics */
  
  assert(obj1->get_nb_level_1_groups() == obj2->get_nb_level_1_groups());
  assert(obj1->get_nb_level_2_groups() == obj2->get_nb_level_2_groups());
  assert(obj1->get_nb_no_level_groups() == obj2->get_nb_no_level_groups());
  
  assert(obj1->get_nb_level_0_cells() == obj2->get_nb_level_0_cells());
  assert(obj1->get_nb_level_1_cells() == obj2->get_nb_level_1_cells());
  assert(obj1->get_nb_level_2_cells() == obj2->get_nb_level_2_cells());
  assert(obj1->get_nb_no_level_cells() == obj2->get_nb_no_level_cells());
  
  assert(obj1->get_nb_group_appearances() == obj2->get_nb_group_appearances());
  assert(obj1->get_nb_group_extinctions() == obj2->get_nb_group_extinctions());
  
  assert(obj1->get_mean_group_lifespan() == obj2->get_mean_group_lifespan());
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test Node equality
 * \details  --
 * \param    Node* obj1
 * \param    Node* obj2
 * \return   \e void
 */
void IntegratedTests::Node_isEqualTo( Node* obj1, Node* obj2 )
{
  assert(obj1->get_id() == obj2->get_id());
  assert(obj1->get_number_of_children() == obj2->get_number_of_children());
  assert(obj1->get_node_state() == obj2->get_node_state());
  assert(obj1->isMasterRoot() == obj2->isMasterRoot());
  assert(obj1->isRoot() == obj2->isRoot());
  assert(obj1->isNormal() == obj2->isNormal());
  assert(obj1->isDead() == obj2->isDead());
  assert(obj1->isAlive() == obj2->isAlive());
  assert(obj1->isTagged() == obj2->isTagged());
  if (obj1->get_alive_cell() == NULL)
  {
    assert(obj2->get_alive_cell() == NULL);
  }
  else
  {
    Cell_isEqualTo(obj1->get_alive_cell(), obj2->get_alive_cell());
  }
  if (obj1->get_replication_report() == NULL)
  {
    assert(obj2->get_replication_report() == NULL);
  }
  else
  {
    ReplicationReport_isEqualTo(obj1->get_replication_report(), obj2->get_replication_report());
  }
  if (obj1->get_parent() == NULL)
  {
    assert(obj2->get_parent() == NULL);
  }
  else
  {
    assert(obj1->get_parent()->get_id() == obj2->get_parent()->get_id());
  }
  if (obj1->get_number_of_children() > 0)
  {
    for (size_t i = 0; i < obj1->get_number_of_children(); i++)
    {
      assert(obj1->get_child(i)->get_id() == obj2->get_child(i)->get_id());
    }
  }
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test Tree equality
 * \details  --
 * \param    Tree* obj1
 * \param    Tree* obj2
 * \return   \e void
 */
void IntegratedTests::Tree_isEqualTo( Tree* obj1, Tree* obj2 )
{
  assert(obj1->get_current_id() == obj2->get_current_id());
  assert(obj1->get_number_of_nodes() == obj2->get_number_of_nodes());
  Node* node1 = obj1->get_first_node();
  Node* node2 = obj2->get_first_node();
  while (node1 != NULL || node2 != NULL)
  {
    assert(node1 != NULL);
    assert(node2 != NULL);
    Node_isEqualTo(node1, node2);
    Node_isEqualTo(obj1->get_node(node1->get_id()), obj2->get_node(node1->get_id()));
    node1 = obj1->get_next_node();
    node2 = obj2->get_next_node();
  }
  std::vector<unsigned long long int> alive_nodes1;
  std::vector<unsigned long long int> alive_nodes2;
  obj1->get_alive_nodes(&alive_nodes1);
  obj2->get_alive_nodes(&alive_nodes2);
  assert(alive_nodes1.size() == alive_nodes2.size());
  for (size_t i = 0; i < alive_nodes1.size(); i++)
  {
    assert(alive_nodes1[i] == alive_nodes2[i]);
    Node_isEqualTo(obj1->get_node(alive_nodes1[i]), obj2->get_node(alive_nodes2[i]));
  }
  if (obj1->get_best_alive_node() == NULL)
  {
    assert(obj2->get_best_alive_node() == NULL);
  }
  else
  {
    Node_isEqualTo(obj1->get_best_alive_node(), obj2->get_best_alive_node());
  }
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test Statistics equality
 * \details  --
 * \param    Statistics* obj1
 * \param    Statistics* obj2
 * \return   \e void
 */
void IntegratedTests::Statistics_isEqualTo( Statistics* obj1, Statistics* obj2 )
{
  /*------------------------------------------------------------------ MEAN statistical variables */
  
  /* PHENOTYPE */
  assert(obj1->_mean_generations == obj2->_mean_generations);
  assert(obj1->_mean_inherited_TF_amount == obj2->_mean_inherited_TF_amount);
  assert(obj1->_mean_inherited_E_amount == obj2->_mean_inherited_E_amount);
  assert(obj1->_mean_TF_amount == obj2->_mean_TF_amount);
  assert(obj1->_mean_E_amount == obj2->_mean_E_amount);
  assert(obj1->_mean_inherited_metabolic_amount == obj2->_mean_inherited_metabolic_amount);
  assert(obj1->_mean_metabolic_amount == obj2->_mean_metabolic_amount);
  assert(obj1->_mean_energy == obj2->_mean_energy);
  assert(obj1->_mean_score == obj2->_mean_score);
  assert(obj1->_mean_lifespan == obj2->_mean_lifespan);
  assert(obj1->_mean_number_of_divisions == obj2->_mean_number_of_divisions);
  assert(obj1->_mean_toxicity == obj2->_mean_toxicity);
  assert(obj1->_mean_metabolic_uptake == obj2->_mean_metabolic_uptake);
  assert(obj1->_mean_metabolic_release == obj2->_mean_metabolic_release);
  assert(obj1->_mean_metabolic_growth_rate == obj2->_mean_metabolic_growth_rate);
  assert(obj1->_mean_Dmetabolic_growth_rate == obj2->_mean_Dmetabolic_growth_rate);
  assert(obj1->_mean_grn_nb_nodes == obj2->_mean_grn_nb_nodes);
  assert(obj1->_mean_grn_nb_edges == obj2->_mean_grn_nb_edges);
  assert(obj1->_mean_metabolic_nb_nodes == obj2->_mean_metabolic_nb_nodes);
  assert(obj1->_mean_metabolic_nb_edges == obj2->_mean_metabolic_nb_edges);
  assert(obj1->_mean_regulation_redundancy == obj2->_mean_regulation_redundancy);
  assert(obj1->_mean_metabolic_redundancy == obj2->_mean_metabolic_redundancy);
  
  /* GENOME STRUCTURE */
  assert(obj1->_mean_genome_size == obj2->_mean_genome_size);
  assert(obj1->_mean_functional_size == obj2->_mean_functional_size);
  assert(obj1->_mean_genome_nb_NC == obj2->_mean_genome_nb_NC);
  assert(obj1->_mean_genome_nb_E == obj2->_mean_genome_nb_E);
  assert(obj1->_mean_genome_nb_TF == obj2->_mean_genome_nb_TF);
  assert(obj1->_mean_genome_nb_BS == obj2->_mean_genome_nb_BS);
  assert(obj1->_mean_genome_nb_P == obj2->_mean_genome_nb_P);
  assert(obj1->_mean_genome_nb_inner_enzymes == obj2->_mean_genome_nb_inner_enzymes);
  assert(obj1->_mean_genome_nb_inflow_pumps == obj2->_mean_genome_nb_inflow_pumps);
  assert(obj1->_mean_genome_nb_outflow_pumps == obj2->_mean_genome_nb_outflow_pumps);
  assert(obj1->_mean_genome_nb_functional_regions == obj2->_mean_genome_nb_functional_regions);
  assert(obj1->_mean_genome_nb_enhancers == obj2->_mean_genome_nb_enhancers);
  assert(obj1->_mean_genome_nb_operators == obj2->_mean_genome_nb_operators);
  assert(obj1->_mean_genome_nb_E_regions == obj2->_mean_genome_nb_E_regions);
  assert(obj1->_mean_genome_nb_TF_regions == obj2->_mean_genome_nb_TF_regions);
  assert(obj1->_mean_genome_nb_mixed_regions == obj2->_mean_genome_nb_mixed_regions);
  assert(obj1->_mean_genome_functional_region_size == obj2->_mean_genome_functional_region_size);
  assert(obj1->_mean_genome_E_region_size == obj2->_mean_genome_E_region_size);
  assert(obj1->_mean_genome_TF_region_size == obj2->_mean_genome_TF_region_size);
  assert(obj1->_mean_genome_mixed_region_size == obj2->_mean_genome_mixed_region_size);
  assert(obj1->_mean_genome_enhancer_size == obj2->_mean_genome_enhancer_size);
  assert(obj1->_mean_genome_operator_size == obj2->_mean_genome_operator_size);
  assert(obj1->_mean_genome_operon_size == obj2->_mean_genome_operon_size);
  assert(obj1->_mean_genome_E_operon_size == obj2->_mean_genome_E_operon_size);
  assert(obj1->_mean_genome_TF_operon_size == obj2->_mean_genome_TF_operon_size);
  assert(obj1->_mean_genome_mixed_operon_size == obj2->_mean_genome_mixed_operon_size);
  
  /* INHERITED STRUCTURE */
  assert(obj1->_mean_inherited_size == obj2->_mean_inherited_size);
  assert(obj1->_mean_inherited_nb_E == obj2->_mean_inherited_nb_E);
  assert(obj1->_mean_inherited_nb_TF == obj2->_mean_inherited_nb_TF);
  assert(obj1->_mean_inherited_nb_inner_enzymes == obj2->_mean_inherited_nb_inner_enzymes);
  assert(obj1->_mean_inherited_nb_inflow_pumps == obj2->_mean_inherited_nb_inflow_pumps);
  assert(obj1->_mean_inherited_nb_outflow_pumps == obj2->_mean_inherited_nb_outflow_pumps);
  
  /*------------------------------------------------------------------ VARIANCE statistical variables */
  
  /* PHENOTYPE */
  assert(obj1->_var_generations == obj2->_var_generations);
  assert(obj1->_var_inherited_TF_amount == obj2->_var_inherited_TF_amount);
  assert(obj1->_var_inherited_E_amount == obj2->_var_inherited_E_amount);
  assert(obj1->_var_TF_amount == obj2->_var_TF_amount);
  assert(obj1->_var_E_amount == obj2->_var_E_amount);
  assert(obj1->_var_inherited_metabolic_amount == obj2->_var_inherited_metabolic_amount);
  assert(obj1->_var_metabolic_amount == obj2->_var_metabolic_amount);
  assert(obj1->_var_energy == obj2->_var_energy);
  assert(obj1->_var_score == obj2->_var_score);
  assert(obj1->_var_lifespan == obj2->_var_lifespan);
  assert(obj1->_var_number_of_divisions == obj2->_var_number_of_divisions);
  assert(obj1->_var_toxicity == obj2->_var_toxicity);
  assert(obj1->_var_metabolic_uptake == obj2->_var_metabolic_uptake);
  assert(obj1->_var_metabolic_release == obj2->_var_metabolic_release);
  assert(obj1->_var_metabolic_growth_rate == obj2->_var_metabolic_growth_rate);
  assert(obj1->_var_Dmetabolic_growth_rate == obj2->_var_Dmetabolic_growth_rate);
  assert(obj1->_var_grn_nb_nodes == obj2->_var_grn_nb_nodes);
  assert(obj1->_var_grn_nb_edges == obj2->_var_grn_nb_edges);
  assert(obj1->_var_metabolic_nb_nodes == obj2->_var_metabolic_nb_nodes);
  assert(obj1->_var_metabolic_nb_edges == obj2->_var_metabolic_nb_edges);
  assert(obj1->_var_regulation_redundancy == obj2->_var_regulation_redundancy);
  assert(obj1->_var_metabolic_redundancy == obj2->_var_metabolic_redundancy);
  
  /* GENOME STRUCTURE */
  assert(obj1->_var_genome_size == obj2->_var_genome_size);
  assert(obj1->_var_functional_size == obj2->_var_functional_size);
  assert(obj1->_var_genome_nb_NC == obj2->_var_genome_nb_NC);
  assert(obj1->_var_genome_nb_E == obj2->_var_genome_nb_E);
  assert(obj1->_var_genome_nb_TF == obj2->_var_genome_nb_TF);
  assert(obj1->_var_genome_nb_BS == obj2->_var_genome_nb_BS);
  assert(obj1->_var_genome_nb_P == obj2->_var_genome_nb_P);
  assert(obj1->_var_genome_nb_inner_enzymes == obj2->_var_genome_nb_inner_enzymes);
  assert(obj1->_var_genome_nb_inflow_pumps == obj2->_var_genome_nb_inflow_pumps);
  assert(obj1->_var_genome_nb_outflow_pumps == obj2->_var_genome_nb_outflow_pumps);
  assert(obj1->_var_genome_nb_functional_regions == obj2->_var_genome_nb_functional_regions);
  assert(obj1->_var_genome_nb_enhancers == obj2->_var_genome_nb_enhancers);
  assert(obj1->_var_genome_nb_operators == obj2->_var_genome_nb_operators);
  assert(obj1->_var_genome_nb_E_regions == obj2->_var_genome_nb_E_regions);
  assert(obj1->_var_genome_nb_TF_regions == obj2->_var_genome_nb_TF_regions);
  assert(obj1->_var_genome_nb_mixed_regions == obj2->_var_genome_nb_mixed_regions);
  assert(obj1->_var_genome_functional_region_size == obj2->_var_genome_functional_region_size);
  assert(obj1->_var_genome_E_region_size == obj2->_var_genome_E_region_size);
  assert(obj1->_var_genome_TF_region_size == obj2->_var_genome_TF_region_size);
  assert(obj1->_var_genome_mixed_region_size == obj2->_var_genome_mixed_region_size);
  assert(obj1->_var_genome_enhancer_size == obj2->_var_genome_enhancer_size);
  assert(obj1->_var_genome_operator_size == obj2->_var_genome_operator_size);
  assert(obj1->_var_genome_operon_size == obj2->_var_genome_operon_size);
  assert(obj1->_var_genome_E_operon_size == obj2->_var_genome_E_operon_size);
  assert(obj1->_var_genome_TF_operon_size == obj2->_var_genome_TF_operon_size);
  assert(obj1->_var_genome_mixed_operon_size == obj2->_var_genome_mixed_operon_size);
  
  /* INHERITED STRUCTURE */
  assert(obj1->_var_inherited_size == obj2->_var_inherited_size);
  assert(obj1->_var_inherited_nb_E == obj2->_var_inherited_nb_E);
  assert(obj1->_var_inherited_nb_TF == obj2->_var_inherited_nb_TF);
  assert(obj1->_var_inherited_nb_inner_enzymes == obj2->_var_inherited_nb_inner_enzymes);
  assert(obj1->_var_inherited_nb_inflow_pumps == obj2->_var_inherited_nb_inflow_pumps);
  assert(obj1->_var_inherited_nb_outflow_pumps == obj2->_var_inherited_nb_outflow_pumps);
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test Simulation equality
 * \details  --
 * \param    Simulation* obj1
 * \param    Simulation* obj2
 * \return   \e void
 */
void IntegratedTests::Simulation_isEqualTo( Simulation* obj1, Simulation* obj2 )
{
  Parameters_isEqualTo(obj1->get_parameters(), obj2->get_parameters());
  Population_isEqualTo(obj1->get_population(), obj2->get_population());
  Environment_isEqualTo(obj1->get_environment(), obj2->get_environment());
  TrophicNetwork_isEqualTo(obj1->get_trophic_network(), obj2->get_trophic_network());
  Tree_isEqualTo(obj1->get_lineage_tree(), obj2->get_lineage_tree());
  Tree_isEqualTo(obj1->get_phylogenetic_tree(), obj2->get_phylogenetic_tree());
  Statistics_isEqualTo(obj1->get_statistics(), obj2->get_statistics());
  assert(obj1->get_min_score() == obj2->get_min_score());
  assert(obj1->get_max_score() == obj2->get_max_score());
  assert(obj1->get_width() == obj2->get_width());
  assert(obj1->get_height() == obj2->get_height());
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test Prng equality
 * \details  --
 * \param    Prng* obj1
 * \param    Prng* obj2
 * \return   \e void
 */
void IntegratedTests::Prng_isEqualTo( Prng* obj1, Prng* obj2 )
{
  size_t test_size = 10000;
  
  /*-----------------------------*/
  /* 1) test uniform law         */
  /*-----------------------------*/
  for (size_t i = 0; i < test_size; i++)
  {
    assert(obj1->uniform() == obj2->uniform());
  }
  
  /*-----------------------------*/
  /* 2) test integer uniform law */
  /*-----------------------------*/
  for (size_t i = 0; i < test_size; i++)
  {
    assert(obj1->uniform(-10000, 10000) == obj2->uniform(-10000, 10000));
  }
  
  /*-----------------------------*/
  /* 3) test bernouilli law      */
  /*-----------------------------*/
  for (size_t i = 0; i < test_size; i++)
  {
    assert(obj1->bernouilli(0.5) == obj2->bernouilli(0.5));
  }
  
  /*-----------------------------*/
  /* 4) test binomial law        */
  /*-----------------------------*/
  for (size_t i = 0; i < test_size; i++)
  {
    assert(obj1->binomial(1000, 0.5) == obj2->binomial(1000, 0.5));
  }
  
  /*-----------------------------*/
  /* 5) test multinomial law     */
  /*-----------------------------*/
  for (size_t i = 0; i < test_size; i++)
  {
    unsigned int draws1[10];
    unsigned int draws2[10];
    double probas[10] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    int N = 1000, K = 10;
    obj1->multinomial(draws1, probas, N, K);
    obj2->multinomial(draws2, probas, N, K);
    for (size_t j = 0; j < 10; j++)
    {
      assert(draws1[j] == draws2[j]);
    }
  }
  
  /*-----------------------------*/
  /* 6) test gaussian law        */
  /*-----------------------------*/
  for (size_t i = 0; i < test_size; i++)
  {
    assert(obj1->gaussian(0, 1) == obj2->gaussian(0, 1));
  }
  
  /*-----------------------------*/
  /* 7) test exponential law     */
  /*-----------------------------*/
  for (size_t i = 0; i < test_size; i++)
  {
    assert(obj1->exponential(10.0) == obj2->exponential(10.0));
  }
  
  /*-----------------------------*/
  /* 8) test exponential law     */
  /*-----------------------------*/
  for (size_t i = 0; i < test_size; i++)
  {
    double probas[10] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double sum = 5.0;
    int N = 10;
    assert(obj1->roulette_wheel(probas, sum, N) == obj2->roulette_wheel(probas, sum, N));
    
    (void)probas;
    (void)sum;
    (void)N;
  }
}

/**
 * \brief    Generate a random parameters set
 * \details  --
 * \param    Parameters* parameters
 * \return   \e void
 */
void IntegratedTests::generate_random_parameters( Parameters* parameters, Prng* prng )
{
  /*------------------------------------------------------------------ parallel computing */
  
  parameters->set_parallel_computing((prng->uniform() < 0.5 ? false : true));
  
  /*------------------------------------------------------------------ simulation schemes */
  
  parameters->set_energy_costs_scheme((prng->uniform() < 0.5 ? false : true));
  parameters->set_membrane_permeability_scheme((prng->uniform() < 0.5 ? false : true));
  parameters->set_metabolic_inheritance_scheme((prng->uniform() < 0.5 ? false : true));
  parameters->set_enzymatic_inheritance_scheme((prng->uniform() < 0.5 ? false : true));
  parameters->set_co_enzyme_activity_scheme((prng->uniform() < 0.5 ? false : true));
  int draw = prng->uniform(1, 3);
  switch (draw)
  {
    case 1:
      parameters->set_score_scheme(ESSENTIAL_METABOLITES_SUM);
      break;
    case 2:
      parameters->set_score_scheme(ESSENTIAL_METABOLITES_SUM_MINUS_DEVIATION);
      break;
    case 3:
      parameters->set_score_scheme(ESSENTIAL_METABOLITES_COMBINATORIAL_CONTRIBUTION);
      break;
  }
  parameters->set_selection_threshold(prng->uniform());
  
  /*------------------------------------------------------------------ genome */
  
  //parameters->set_cell_metabolic_range( const distribution_law* range_law );
  
  parameters->set_transcription_factor_binding_window(prng->uniform(0, 10));
  
  parameters->set_initial_number_of_NC_units(prng->uniform(0, 30));
  parameters->set_initial_number_of_E_units(prng->uniform(0, 30));
  parameters->set_initial_number_of_TF_units(prng->uniform(0, 30));
  parameters->set_initial_number_of_BS_units(prng->uniform(0, 30));
  parameters->set_initial_number_of_P_units(prng->uniform(0, 30));
  
  parameters->set_point_mutation_rate(pow(10.0, prng->uniform(-6, -2)));
  parameters->set_duplication_rate(pow(10.0, prng->uniform(-6, -2)));
  parameters->set_deletion_rate(pow(10.0, prng->uniform(-6, -2)));
  parameters->set_translocation_rate(pow(10.0, prng->uniform(-6, -2)));
  parameters->set_inversion_rate(pow(10.0, prng->uniform(-6, -2)));
  parameters->set_transition_rate(pow(10.0, prng->uniform(-6, -2)));
  parameters->set_breakpoint_rate(pow(10.0, prng->uniform(-6, -2)));
  
  parameters->set_substrate_tag_mutation_size(prng->uniform(0, 50));
  parameters->set_product_tag_mutation_size(prng->uniform(0, 50));
  parameters->set_kcat_mutation_size(prng->uniform());
  parameters->set_kcat_km_ratio_mutation_size(prng->uniform());
  parameters->set_binding_site_tag_mutation_size(prng->uniform(0, 50));
  parameters->set_co_enzyme_tag_mutation_size(prng->uniform(0, 50));
  parameters->set_transcription_factor_tag_mutation_size(prng->uniform(0, 50));
  parameters->set_basal_expression_level_mutation_size(prng->uniform());
  
  /*------------------------------------------------------------------ genetic regulation network */
  
  parameters->set_genetic_regulation_network_timestep((double)prng->uniform(1, 100));
  //obj->set_hill_function_theta(0.1234);
  //obj->set_hill_function_n(0.1234);
  parameters->set_protein_degradation_rate(prng->uniform());
  
  /*------------------------------------------------------------------ metabolic network */
  
  parameters->set_metabolism_timestep((double)prng->uniform(1, 100));
  
  parameters->set_essential_metabolites_toxicity_threshold(pow(10.0, prng->uniform(-1, 2)));
  parameters->set_non_essential_metabolites_toxicity_threshold(pow(10.0, prng->uniform(-1, 2)));
  
  //obj->set_initial_metabolites_amount_in_cells(0.1234);
  
  //parameters->set_energy_toxicity_threshold( double energy_toxicity_threshold );
  
  /*------------------------------------------------------------------ energy */
  
  //parameters->set_energy_transcription_cost(-0.001234);
  //parameters->set_energy_degradation_cost(-0.001234);
  //parameters->set_energy_enzymatic_cost(-0.001234);
  //parameters->set_energy_pumping_cost(-0.001234);
  
  //parameters->set_energy_dissipation_rate(0.1234);
  
  //parameters->set_energy_toxicity_threshold(12.34);
  
  /*------------------------------------------------------------------ cell */
  
  parameters->set_death_probability(pow(10.0, prng->uniform(-3, -1)));
  
  parameters->set_membrane_permeability(0.0);
  if (prng->uniform() < 0.5)
  {
    parameters->set_membrane_permeability(prng->uniform());
  }
  
  /*------------------------------------------------------------------ population */
  
  parameters->set_migration_rate(prng->uniform());
  parameters->set_hgt_rate(pow(10.0, prng->uniform(-6, -2)));
  
  /*------------------------------------------------------------------ environment */
  
  environment_properties* properties = parameters->get_environment_properties();
  
  double conc = pow(10.0, prng->uniform(-1, 2));
  properties->concentration_range.min = conc;
  properties->concentration_range.max = conc;
  properties->interaction_scheme = (prng->uniform() < 0.5 ? INTERACTION : NO_INTERACTION);
  properties->renewal_scheme     = (prng->uniform() < 0.5 ? KEEP_MATTER : CLEAR_MATTER);
  
  draw = prng->uniform(1, 3);
  switch (draw)
  {
    case 1:
      properties->variation_scheme = RANDOM_SCHEME;
      break;
    case 2:
      properties->variation_scheme = PERIODIC_SCHEME;
      break;
    case 3:
      properties->variation_scheme = CYCLIC_SCHEME;
      break;
  }
  draw = prng->uniform(1, 4);
  switch (draw)
  {
    case 1:
      properties->localization_scheme = GLOBAL_LOCALIZATION;
      break;
    case 2:
      properties->localization_scheme = RANDOM_LOCALIZATION;
      break;
    case 3:
      properties->localization_scheme = SPOT_LOCALIZATION;
      break;
    case 4:
      properties->localization_scheme = CENTER_LOCALIZATION;
      break;
  }
  properties->introduction_rate     = prng->uniform();
  properties->diffusion_coefficient = pow(10.0, prng->uniform(-5, -1));
  properties->degradation_rate      = pow(10.0, prng->uniform(-5, -1));
}
