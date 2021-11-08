
/**
 * \file      UnitaryTests.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     UnitaryTests class definition
 */

/****************************************************************************
 * Evo2Sim (Evolution of Evolution Simulator)
 * -------------------------------------------
 * Digital evolution model dedicated to
 * bacterial in silico experimental evolution.
 *
 * Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon
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

#include "UnitaryTests.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    void
 * \return   \e void
 */
UnitaryTests::UnitaryTests( Parameters* parameters )
{
  _parameters = parameters;
  initialize_mutation_rates();
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
UnitaryTests::~UnitaryTests( void )
{
  delete[] _mutation_rates;
  _mutation_rates = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Run unitary tests
 * \details  --
 * \param    std::string parameters_filename
 * \return   \e void
 */
void UnitaryTests::run_unitary_tests( std::string parameters_filename )
{
  printf("> Run unitary tests:\n");
  
  printf("  > Testing Parameters class ...");
  test_Parameters_class(parameters_filename);
  printf(" ok.\n");
  
  printf("  > Testing MutationVector class ...");
  test_MutationVector_class();
  printf(" ok.\n");
  
  printf("  > Testing MutationEvent class ...");
  test_MutationEvent_class();
  printf(" ok.\n");
  
  printf("  > Testing ReplicationReport class ...");
  test_ReplicationReport_class();
  printf(" ok.\n");
  
  printf("  > Testing Genome class ...");
  test_Genome_class();
  printf(" ok.\n");
  
  printf("  > Testing InheritedProteins class ...");
  test_InheritedProteins_class();
  printf(" ok.\n");
  
  printf("  > Testing SpeciesList class ...");
  test_SpeciesList_class();
  printf(" ok.\n");
  
  printf("  > Testing Prng class ...");
  test_Prng_class();
  printf(" ok.\n");
  
  printf("> End of unitary tests.\n\n");
}

/**
 * \brief    Test Parameters class
 * \details  --
 * \param    std::string filename
 * \return   \e void
 */
void UnitaryTests::test_Parameters_class( std::string filename )
{
  Parameters* parameters1 = new Parameters();
  parameters1->load_parameters_from_file(filename);
  
  /*------------------------------------------*/
  /* Test Parameters backup constructor       */
  /*------------------------------------------*/
  mkdir("parameters", 0777);
  mkdir("prng", 0777);
  parameters1->save(1000);
  
  Parameters* parameters2 = new Parameters(1000);
  
  Parameters_isEqualTo(parameters2, parameters1);
  
  delete parameters1;
  parameters1 = NULL;
  
  /*------------------------------------------*/
  /* Test Parameters copy constructor         */
  /*------------------------------------------*/
  modify_Parameters(parameters2);
  Parameters* parameters3 = new Parameters(*parameters2);
  
  Parameters_isEqualTo(parameters3, parameters2);
  
  delete parameters2;
  parameters2 = NULL;
  
  /*------------------------------------------*/
  /* Test Parameters file writing and loading */
  /*------------------------------------------*/
  parameters3->write("./parameters/tmp.txt");
  
  Parameters* parameters4 = new Parameters();
  parameters4->load_parameters_from_file("./parameters/tmp.txt");
  
  Parameters_isEqualTo(parameters4, parameters3);
  
  int success = 0;
  success = system("rm -rf parameters");
  success = system("rm -rf prng");
  (void)success;
  
  delete parameters3;
  parameters3 = NULL;
  delete parameters4;
  parameters4 = NULL;
}

/**
 * \brief    Test MutationVector class
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::test_MutationVector_class( void )
{
  Prng* prng = new Prng();
  prng->set_seed((unsigned long int)time(NULL));
  
  MutationVector* dX = new MutationVector();
  
  for (size_t i = 0; i < 100; i++)
  {
    while (!dX->draw(prng, _mutation_rates));
  }
  
  /*----------------------------------------*/
  /* Test struct native copy constructor    */
  /*----------------------------------------*/
  const genetic_unit unit = *dX->get_dX();
  genetic_unit unit1 = unit;
  genetic_unit_isEqualTo(&unit1, dX->get_dX());
  
  /*----------------------------------------*/
  /* Test MutationVector backup constructor */
  /*----------------------------------------*/
  gzFile backup_file = gzopen("tmp", "w");
  dX->save(backup_file);
  gzclose(backup_file);
  
  backup_file = gzopen("tmp", "r");
  MutationVector* dX2 = new MutationVector(backup_file);
  gzclose(backup_file);
  
  MutationVector_isEqualTo(dX2, dX);
  
  int success = system("rm -rf tmp");
  (void)success;
  
  delete dX;
  dX = NULL;
  
  /*----------------------------------------*/
  /* Test MutationVector copy constructor   */
  /*----------------------------------------*/
  MutationVector* dX3 = new MutationVector(*dX2);
  
  MutationVector_isEqualTo(dX3, dX2);
  
  delete dX2;
  dX2 = NULL;
  delete dX3;
  dX3 = NULL;
  delete prng;
  prng = NULL;
}

/**
 * \brief    Test MutationEvent class
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::test_MutationEvent_class( void )
{
  Prng* prng = new Prng();
  prng->set_seed((unsigned long int)time(NULL));
  
  MutationVector* dX = new MutationVector();
  for (size_t i = 0; i < 100; i++)
  {
    while (!dX->draw(prng, _mutation_rates));
  }
  MutationEvent* event1 = new MutationEvent(POINT_MUTATION, 100, dX);
  
  /*---------------------------------------*/
  /* Test MutationEvent backup constructor */
  /*---------------------------------------*/
  gzFile backup_file = gzopen("tmp", "w");
  event1->save(backup_file);
  gzclose(backup_file);
  
  backup_file = gzopen("tmp", "r");
  MutationEvent* event2 = new MutationEvent(backup_file);
  gzclose(backup_file);
  
  MutationEvent_isEqualTo(event2, event1);
  
  int success = system("rm -rf tmp");
  (void)success;
  
  delete event1;
  event1 = NULL;
  
  /*---------------------------------------*/
  /* Test MutationEvent copy constructor   */
  /*---------------------------------------*/
  MutationEvent* event3 = new MutationEvent(*event2);
  
  MutationEvent_isEqualTo(event3, event2);
  
  delete event2;
  event2 = NULL;
  delete event3;
  event3 = NULL;
  delete prng;
  prng = NULL;
}

/**
 * \brief    Test ReplicationReport class
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::test_ReplicationReport_class( void )
{
  Prng* prng = new Prng();
  prng->set_seed((unsigned long int)time(NULL));
  
  /* 1) create a mutation event ------*/
  MutationVector* dX = new MutationVector();
  for (size_t i = 0; i < 100; i++)
  {
    while (!dX->draw(prng, _mutation_rates));
  }
  MutationEvent* point_mutation_event = new MutationEvent(POINT_MUTATION, 100, dX);
  MutationEvent* hgt_event            = new MutationEvent(HGT, 50, 25, 5, 5, 5, 5, 5);
  MutationEvent* duplication_event    = new MutationEvent(DUPLICATION, 100, 200, 300, 250, 50, 50, 50, 50, 50);
  MutationEvent* deletion_event       = new MutationEvent(DELETION, 100, 200, 300, 250, 50, 50, 50, 50, 50);
  MutationEvent* translocation_event  = new MutationEvent(TRANSLOCATION, 100, 200, 300, 250, 50, 50, 50, 50, 50);
  MutationEvent* inversion_event      = new MutationEvent(INVERSION, 100, 200, 300, 250, 50, 50, 50, 50, 50);
  
  /* 2) create the replication report */
  ReplicationReport* report1 = new ReplicationReport();
  report1->add_mutation_event(point_mutation_event);
  report1->add_mutation_event(hgt_event);
  report1->add_mutation_event(duplication_event);
  report1->add_mutation_event(deletion_event);
  report1->add_mutation_event(translocation_event);
  report1->add_mutation_event(inversion_event);
  report1->compute_mean();
  
  /*-------------------------------------------*/
  /* Test ReplicationReport atributes          */
  /*-------------------------------------------*/
  
  /*------------------------------------------------------------------ global data */
  
  assert(report1->get_number_of_events() == 6);
  
  /*------------------------------------------------------------------ point mutations data */
  
  assert(report1->get_nb_point_mutations() == 1);
  
  size_t nb_NC_point_mutations   = 0;
  size_t nb_E_point_mutations    = 0;
  size_t nb_TF_point_mutations   = 0;
  size_t nb_BS_point_mutations   = 0;
  size_t nb_P_point_mutations    = 0;
  size_t nb_NC_to_E_transitions  = 0;
  size_t nb_NC_to_TF_transitions = 0;
  size_t nb_NC_to_BS_transitions = 0;
  size_t nb_NC_to_P_transitions  = 0;
  size_t nb_E_to_NC_transitions  = 0;
  size_t nb_E_to_TF_transitions  = 0;
  size_t nb_E_to_BS_transitions  = 0;
  size_t nb_E_to_P_transitions   = 0;
  size_t nb_TF_to_NC_transitions = 0;
  size_t nb_TF_to_E_transitions  = 0;
  size_t nb_TF_to_BS_transitions = 0;
  size_t nb_TF_to_P_transitions  = 0;
  size_t nb_BS_to_NC_transitions = 0;
  size_t nb_BS_to_E_transitions  = 0;
  size_t nb_BS_to_TF_transitions = 0;
  size_t nb_BS_to_P_transitions  = 0;
  size_t nb_P_to_NC_transitions  = 0;
  size_t nb_P_to_E_transitions   = 0;
  size_t nb_P_to_TF_transitions  = 0;
  size_t nb_P_to_BS_transitions  = 0;
  
  genetic_unit_type utype = point_mutation_event->get_mutation_vector()->get_dX()->type;
  switch (utype)
  {
    case NON_CODING:
      nb_NC_point_mutations   = 1;
      break;
    case ENZYME:
      nb_E_point_mutations    = 1;
      break;
    case TRANSCRIPTION_FACTOR:
      nb_TF_point_mutations   = 1;
      break;
    case BINDING_SITE:
      nb_BS_point_mutations   = 1;
      break;
    case PROMOTER:
      nb_P_point_mutations    = 1;
      break;
    case NC_TO_E_TRANSITION:
      nb_NC_to_E_transitions  = 1;
      break;
    case NC_TO_TF_TRANSITION:
      nb_NC_to_TF_transitions = 1;
      break;
    case NC_TO_BS_TRANSITION:
      nb_NC_to_BS_transitions = 1;
      break;
    case NC_TO_P_TRANSITION:
      nb_NC_to_P_transitions  = 1;
      break;
    case E_TO_NC_TRANSITION:
      nb_E_to_NC_transitions  = 1;
      break;
    case E_TO_TF_TRANSITION:
      nb_E_to_TF_transitions  = 1;
      break;
    case E_TO_BS_TRANSITION:
      nb_E_to_BS_transitions  = 1;
      break;
    case E_TO_P_TRANSITION:
      nb_E_to_P_transitions   = 1;
      break;
    case TF_TO_NC_TRANSITION:
      nb_TF_to_NC_transitions = 1;
      break;
    case TF_TO_E_TRANSITION:
      nb_TF_to_E_transitions  = 1;
      break;
    case TF_TO_BS_TRANSITION:
      nb_TF_to_BS_transitions = 1;
      break;
    case TF_TO_P_TRANSITION:
      nb_TF_to_P_transitions  = 1;
      break;
    case BS_TO_NC_TRANSITION:
      nb_BS_to_NC_transitions = 1;
      break;
    case BS_TO_E_TRANSITION:
      nb_BS_to_E_transitions  = 1;
      break;
    case BS_TO_TF_TRANSITION:
      nb_BS_to_TF_transitions = 1;
      break;
    case BS_TO_P_TRANSITION:
      nb_BS_to_P_transitions  = 1;
      break;
    case P_TO_NC_TRANSITION:
      nb_P_to_NC_transitions  = 1;
      break;
    case P_TO_E_TRANSITION:
      nb_P_to_E_transitions   = 1;
      break;
    case P_TO_TF_TRANSITION:
      nb_P_to_TF_transitions  = 1;
      break;
    case P_TO_BS_TRANSITION:
      nb_P_to_BS_transitions  = 1;
      break;
    default:
      break;
  }
  
  assert(report1->get_nb_NC_point_mutations() == nb_NC_point_mutations);
  assert(report1->get_nb_E_point_mutations() == nb_E_point_mutations);
  assert(report1->get_nb_TF_point_mutations() == nb_TF_point_mutations);
  assert(report1->get_nb_BS_point_mutations() == nb_BS_point_mutations);
  assert(report1->get_nb_P_point_mutations() == nb_P_point_mutations);
  
  assert(report1->get_nb_NC_to_E_transitions() == nb_NC_to_E_transitions);
  assert(report1->get_nb_NC_to_TF_transitions() == nb_NC_to_TF_transitions);
  assert(report1->get_nb_NC_to_BS_transitions() == nb_NC_to_BS_transitions);
  assert(report1->get_nb_NC_to_P_transitions() == nb_NC_to_P_transitions);
  
  assert(report1->get_nb_E_to_NC_transitions() == nb_E_to_NC_transitions);
  assert(report1->get_nb_E_to_TF_transitions() == nb_E_to_TF_transitions);
  assert(report1->get_nb_E_to_BS_transitions() == nb_E_to_BS_transitions);
  assert(report1->get_nb_E_to_P_transitions() == nb_E_to_P_transitions);
  
  assert(report1->get_nb_TF_to_NC_transitions() == nb_TF_to_NC_transitions);
  assert(report1->get_nb_TF_to_E_transitions() == nb_TF_to_E_transitions);
  assert(report1->get_nb_TF_to_BS_transitions() == nb_TF_to_BS_transitions);
  assert(report1->get_nb_TF_to_P_transitions() == nb_TF_to_P_transitions);
  
  assert(report1->get_nb_BS_to_NC_transitions() == nb_BS_to_NC_transitions);
  assert(report1->get_nb_BS_to_E_transitions() == nb_BS_to_E_transitions);
  assert(report1->get_nb_BS_to_TF_transitions() == nb_BS_to_TF_transitions);
  assert(report1->get_nb_BS_to_P_transitions() == nb_BS_to_P_transitions);
  
  assert(report1->get_nb_P_to_NC_transitions() == nb_P_to_NC_transitions);
  assert(report1->get_nb_P_to_E_transitions() == nb_P_to_E_transitions);
  assert(report1->get_nb_P_to_TF_transitions() == nb_P_to_TF_transitions);
  assert(report1->get_nb_P_to_BS_transitions() == nb_P_to_BS_transitions);
  
  if (point_mutation_event->get_mutation_vector()->get_dX()->type == ENZYME)
  {
    assert(report1->get_mean_s_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->s);
    assert(report1->get_mean_p_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->p);
    assert(report1->get_mean_kcat_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->kcat);
    assert(report1->get_mean_kcat_km_ratio_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->kcat_km_ratio);
    assert(report1->get_mean_BS_tag_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->BS_tag);
    assert(report1->get_mean_coE_tag_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->coE_tag);
    assert(report1->get_mean_TF_tag_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->TF_tag);
    assert(report1->get_mean_basal_expression_level_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->basal_expression_level);
  }
  
  /*------------------------------------------------------------------ HGT data */
  
  assert(report1->get_nb_HGT() == 1);
  assert(report1->get_mean_HGT_size() == 25.0);
  assert(report1->get_nb_NC_HGT() == 5);
  assert(report1->get_nb_E_HGT() == 5);
  assert(report1->get_nb_TF_HGT() == 5);
  assert(report1->get_nb_BS_HGT() == 5);
  assert(report1->get_nb_P_HGT() == 5);
  
  /*------------------------------------------------------------------ rearrangements data */
  
  assert(report1->get_nb_rearrangements() == 4);
  assert(report1->get_nb_duplicated_NC() == 50);
  assert(report1->get_nb_duplicated_E() == 50);
  assert(report1->get_nb_duplicated_TF() == 50);
  assert(report1->get_nb_duplicated_BS() == 50);
  assert(report1->get_nb_duplicated_P() == 50);
  assert(report1->get_nb_deleted_NC() == 50);
  assert(report1->get_nb_deleted_E() == 50);
  assert(report1->get_nb_deleted_TF() == 50);
  assert(report1->get_nb_deleted_BS() == 50);
  assert(report1->get_nb_deleted_P() == 50);
  assert(report1->get_nb_duplications() == 1);
  assert(report1->get_nb_deletions() == 1);
  assert(report1->get_nb_translocations() == 1);
  assert(report1->get_nb_inversions() == 1);
  assert(report1->get_mean_rearrangement_size() == 250.0);
  assert(report1->get_mean_duplication_size() == 250.0);
  assert(report1->get_mean_deletion_size() == 250.0);
  assert(report1->get_mean_translocation_size() == 250.0);
  assert(report1->get_mean_inversion_size() == 250.0);
  
  /*-------------------------------------------*/
  /* Test ReplicationReport backup constructor */
  /*-------------------------------------------*/
  gzFile backup_file = gzopen("tmp", "w");
  report1->save(backup_file);
  gzclose(backup_file);
  
  backup_file = gzopen("tmp", "r");
  ReplicationReport* report2 = new ReplicationReport(backup_file);
  gzclose(backup_file);
  
  ReplicationReport_isEqualTo(report2, report1);
  
  int success = system("rm -rf tmp");
  (void)success;
  
  delete report1;
  report1 = NULL;
  
  /*-------------------------------------------*/
  /* Test ReplicationReport copy constructor   */
  /*-------------------------------------------*/
  ReplicationReport* report3 = new ReplicationReport(*report2);
  
  ReplicationReport_isEqualTo(report3, report2);
  
  delete report2;
  report2 = NULL;
  delete report3;
  report3 = NULL;
  
  (void)nb_NC_point_mutations;
  (void)nb_E_point_mutations;
  (void)nb_TF_point_mutations;
  (void)nb_BS_point_mutations;
  (void)nb_P_point_mutations;
  (void)nb_NC_to_E_transitions;
  (void)nb_NC_to_TF_transitions;
  (void)nb_NC_to_BS_transitions;
  (void)nb_NC_to_P_transitions;
  (void)nb_E_to_NC_transitions;
  (void)nb_E_to_TF_transitions;
  (void)nb_E_to_BS_transitions;
  (void)nb_E_to_P_transitions;
  (void)nb_TF_to_NC_transitions;
  (void)nb_TF_to_E_transitions;
  (void)nb_TF_to_BS_transitions;
  (void)nb_TF_to_P_transitions;
  (void)nb_BS_to_NC_transitions;
  (void)nb_BS_to_E_transitions;
  (void)nb_BS_to_TF_transitions;
  (void)nb_BS_to_P_transitions;
  (void)nb_P_to_NC_transitions;
  (void)nb_P_to_E_transitions;
  (void)nb_P_to_TF_transitions;
  (void)nb_P_to_BS_transitions;
  
  delete prng;
  prng = NULL;
}

/**
 * \brief    Test Genome class
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::test_Genome_class( void )
{
  Prng* prng1 = new Prng();
  Prng* prng2 = new Prng();
  unsigned long int seed = (unsigned long int)time(NULL);
  prng1->set_seed(seed);
  prng2->set_seed(seed);
  
  /* 1) build the genome */
  ReplicationReport* rep_report1 = new ReplicationReport();
  ReplicationReport* rep_report2 = new ReplicationReport();
  
  Genome* genome1 = new Genome(_parameters, prng1, rep_report1);
  
  Parameters* param_copy  = new Parameters(*_parameters);
  Genome*     genome_copy = new Genome(param_copy, prng2, rep_report2);
  
  genetic_unit* unit = new genetic_unit;
  initialize_genetic_unit(unit);
  bool suitable = false;
  while (!suitable)
  {
    for (size_t i = 0; i < 1000; i++)
    {
      genome1->add_genetic_unit(unit);
      genome_copy->add_genetic_unit(unit);
    }
    for (size_t i = 0; i < 10; i++)
    {
      genome1->mutate(_mutation_rates);
      genome_copy->mutate(_mutation_rates);
    }
    if (genome1->get_size() > 0)
    {
      suitable = true;
    }
  }
  delete unit;
  unit = NULL;
  
  /* 2) test random driven equality */
  ReplicationReport_isEqualTo(rep_report1, rep_report2);
  Genome_isEqualTo(genome_copy, genome1);
  delete genome_copy;
  genome_copy = NULL;
  delete param_copy;
  param_copy = NULL;
  
  /*---------------------------------------*/
  /* Test Genome backup constructor        */
  /*---------------------------------------*/
  gzFile backup_file = gzopen("tmp", "w");
  genome1->save(backup_file);
  gzclose(backup_file);
  
  backup_file = gzopen("tmp", "r");
  Genome* genome2 = new Genome(_parameters, prng2, rep_report1, backup_file);
  gzclose(backup_file);
  
  Genome_isEqualTo(genome2, genome1);
  
  int success = system("rm -rf tmp");
  (void)success;
  
  delete genome1;
  genome1 = NULL;
  
  /*---------------------------------------*/
  /* Test Genome copy constructor          */
  /*---------------------------------------*/
  Genome* genome3 = new Genome(*genome2, prng2, rep_report1);
  
  Genome_isEqualTo(genome3, genome2);
  
  delete genome2;
  genome2 = NULL;
  delete genome3;
  genome3 = NULL;
  delete prng1;
  prng1 = NULL;
  delete prng2;
  prng2 = NULL;
}

/**
 * \brief    Test InheritedProteins class
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::test_InheritedProteins_class( void )
{
  /* 1) build the inherited proteins class */
  InheritedProteins* inhprot1 = new InheritedProteins(_parameters);
  
  genetic_unit* unit = new genetic_unit;
  initialize_genetic_unit(unit);
  for (size_t i = 0; i < 1000; i++)
  {
    if (_parameters->get_simulation_prng()->uniform() < 0.5)
    {
      unit->type = ENZYME;
    }
    else
    {
      unit->type = TRANSCRIPTION_FACTOR;
    }
    inhprot1->add_genetic_unit(unit);
  }
  delete unit;
  unit = NULL;
  inhprot1->initialize_concentration_vector();
  inhprot1->build_index_list();
  
  /*---------------------------------------*/
  /* Test Genome backup constructor        */
  /*---------------------------------------*/
  gzFile backup_file = gzopen("tmp", "w");
  inhprot1->save(backup_file);
  gzclose(backup_file);
  
  backup_file = gzopen("tmp", "r");
  InheritedProteins* inhprot2 = new InheritedProteins(_parameters, backup_file);
  gzclose(backup_file);
  
  InheritedProteins_isEqualTo(inhprot2, inhprot1);
  
  int success = system("rm -rf tmp");
  (void)success;
  
  delete inhprot1;
  inhprot1 = NULL;
  
  /*---------------------------------------*/
  /* Test Genome copy constructor          */
  /*---------------------------------------*/
  InheritedProteins* inhprot3 = new InheritedProteins(*inhprot2);
  
  InheritedProteins_isEqualTo(inhprot3, inhprot2);
  
  delete inhprot2;
  inhprot2 = NULL;
  delete inhprot3;
  inhprot3 = NULL;
}

/**
 * \brief    Test SpeciesList class
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::test_SpeciesList_class( void )
{
  SpeciesList* spl1 = new SpeciesList();
  spl1->set(1, 10.0);
  spl1->set(2, 20.0);
  spl1->set(3, 30.0);
  spl1->set(4, 40.0);
  spl1->set(5, 50.0);
  
  /*---------------------------------------*/
  /* Test SpeciesList backup constructor   */
  /*---------------------------------------*/
  gzFile backup_file = gzopen("tmp", "w");
  spl1->save(backup_file);
  gzclose(backup_file);
  
  backup_file = gzopen("tmp", "r");
  SpeciesList* spl2 = new SpeciesList(backup_file);
  gzclose(backup_file);
  
  SpeciesList_isEqualTo(spl1, spl2);
  
  int success = system("rm -rf tmp");
  (void)success;
  
  delete spl1;
  spl1 = NULL;
  
  /*---------------------------------------*/
  /* Test SpeciesList copy constructor     */
  /*---------------------------------------*/
  SpeciesList* spl3 = new SpeciesList(*spl2);
  SpeciesList_isEqualTo(spl2, spl3);
  delete spl2;
  spl2 = NULL;
  delete spl3;
  spl3 = NULL;
}

/**
 * \brief    Test Prng class
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::test_Prng_class( void )
{
  mkdir("prng", 0777);
  unsigned long int seed = (unsigned long int)time(NULL);
  Prng* prng1 = new Prng();
  prng1->set_seed(seed);
  
  /*---------------------------------------*/
  /* Test Prng seed constructor            */
  /*---------------------------------------*/
  Prng* prng2 = new Prng();
  prng2->set_seed(seed);
  Prng_isEqualTo(prng1, prng2);
  delete prng1;
  prng1 = NULL;
  
  /*---------------------------------------*/
  /* Test Prng backup constructor          */
  /*---------------------------------------*/
  std::stringstream prng_file_name;
  prng_file_name << "./prng/prng_" << 10000;
  FILE * prng_file;
  prng_file = fopen(prng_file_name.str().c_str(), "w");
  prng2->save(prng_file);
  fclose(prng_file);
  
  prng_file = fopen(prng_file_name.str().c_str(), "r");
  Prng* prng3 = new Prng(prng_file);
  fclose(prng_file);
  
  Prng_isEqualTo(prng2, prng3);
  delete prng2;
  prng2 = NULL;
  
  /*---------------------------------------*/
  /* Test Prng copy constructor            */
  /*---------------------------------------*/
  Prng* prng4 = new Prng(*prng3);
  Prng_isEqualTo(prng3, prng4);
  delete prng3;
  prng3 = NULL;
  delete prng4;
  prng4 = NULL;
  int success = system("rm -rf prng");
  (void)success;
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Initialize the mutation rates vector
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::initialize_mutation_rates( void )
{
  _mutation_rates = new double[NUMBER_OF_MUTATION_RATES];
  _mutation_rates[POINT_MUTATION_RATE]                    = _parameters->get_point_mutation_rate();
  _mutation_rates[DUPLICATION_RATE]                       = _parameters->get_duplication_rate();
  _mutation_rates[DELETION_RATE]                          = _parameters->get_deletion_rate();
  _mutation_rates[TRANSLOCATION_RATE]                     = _parameters->get_translocation_rate();
  _mutation_rates[INVERSION_RATE]                         = _parameters->get_inversion_rate();
  _mutation_rates[TRANSITION_RATE]                        = _parameters->get_transition_rate();
  _mutation_rates[BREAKPOINT_RATE]                        = _parameters->get_breakpoint_rate();
  _mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE]            = _parameters->get_substrate_tag_mutation_size();
  _mutation_rates[PRODUCT_TAG_MUTATION_SIZE]              = _parameters->get_product_tag_mutation_size();
  _mutation_rates[KCAT_MUTATION_SIZE]                     = _parameters->get_kcat_mutation_size();
  _mutation_rates[KCAT_KM_RATIO_MUTATION_SIZE]            = _parameters->get_kcat_km_ratio_mutation_size();
  _mutation_rates[BINDING_SITE_TAG_MUTATION_SIZE]         = _parameters->get_binding_site_tag_mutation_size();
  _mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE]            = _parameters->get_co_enzyme_tag_mutation_size();
  _mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE] = _parameters->get_transcription_factor_tag_mutation_size();
  _mutation_rates[BASAL_EXPRESSION_LEVEL_MUTATION_SIZE]   = _parameters->get_basal_expression_level_mutation_size();
}

/**
 * \brief    Initialize a default genetic unit
 * \details  --
 * \param    genetic_unit* obj
 * \return   \e void
 */
void UnitaryTests::initialize_genetic_unit( genetic_unit* obj )
{
  /*------------------------------------------------------------------ Global attributes */
  
  obj->type              = NON_CODING;
  obj->identifier        = 0;
  obj->parent_identifier = 0;
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  obj->s             = 1;
  obj->p             = 1;
  obj->kcat          = pow(10.0, KCAT_MIN_LOG);
  obj->kcat_km_ratio = pow(10.0, KCAT_KM_RATIO_MIN_LOG);
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  obj->BS_tag         = 0;
  obj->coE_tag        = 1;
  obj->free_activity  = false;
  obj->bound_activity = false;
  obj->binding_window = 0;
  
  /*------------------------------------------------------------------ Bidnign site type (BS) attributes */
  
  obj->TF_tag = 0;
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  obj->basal_expression_level = 0.0;
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  obj->functional = false;
}

/**
 * \brief    Test Parameters struct equality
 * \details  --
 * \param    Parameters* obj1
 * \param    Parameters* obj2
 * \return   \e void
 */
void UnitaryTests::Parameters_isEqualTo( Parameters* obj1, Parameters* obj2 )
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
  
  /*------------------------------------------------------------------ prng */
  
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
void UnitaryTests::genetic_unit_isEqualTo( genetic_unit* obj1, genetic_unit* obj2 )
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
  
  /*------------------------------------------------------------------ Bidnign site type (BS) attributes */
  
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
void UnitaryTests::MutationVector_isEqualTo( MutationVector* obj1, MutationVector* obj2 )
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
void UnitaryTests::MutationEvent_isEqualTo( MutationEvent* obj1, MutationEvent* obj2 )
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
void UnitaryTests::ReplicationReport_isEqualTo( ReplicationReport* obj1, ReplicationReport* obj2 )
{
  /*------------------------------------------------------------------ Genome size */
  
  assert(obj1->get_old_genome_size() == obj2->get_old_genome_size());
  assert(obj1->get_new_genome_size() == obj2->get_new_genome_size());
  assert(obj1->get_genome_functional_size() == obj2->get_genome_functional_size());
  
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
  
  /*------------------------------------------------------------------ rearrangements data */
  
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
void UnitaryTests::Genome_isEqualTo( Genome* obj1, Genome* obj2 )
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
    assert(obj1->get_concentration_vector()[i] == obj2->get_concentration_vector()[i]);
    assert(obj1->get_concentration_vector()[i] >= 0.0);
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
void UnitaryTests::InheritedProteins_isEqualTo( InheritedProteins* obj1, InheritedProteins* obj2 )
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
 * \brief    Test SpeciesList equality
 * \details  --
 * \param    SpeciesList* obj1
 * \param    SpeciesList* obj2
 * \return   \e void
 */
void UnitaryTests::SpeciesList_isEqualTo( SpeciesList* obj1, SpeciesList* obj2 )
{
  assert(obj1->get_size() == obj2->get_size());
  assert(obj1->get_amount() == obj2->get_amount());
  for (int i = 0; i < (int)obj1->get_size(); i++)
  {
    assert(obj1->get(i+1) == obj2->get(i+1));
    assert(obj1->get_X()[i] == obj2->get_X()[i]);
    assert(obj1->get(i+1) == obj1->get_X()[i]);
    assert(obj2->get(i+1) == obj2->get_X()[i]);
  }
  
  for (int i = 1; i <= 20; i++)
  {
    obj1->add(i, 10.0);
    obj2->add(i, 10.0);
  }
  assert(obj1->get_size() == obj2->get_size());
  assert(obj1->get_amount() == obj2->get_amount());
  for (int i = 0; i < (int)obj1->get_size(); i++)
  {
    assert(obj1->get(i+1) == obj2->get(i+1));
    assert(obj1->get_X()[i] == obj2->get_X()[i]);
    assert(obj1->get(i+1) == obj1->get_X()[i]);
    assert(obj2->get(i+1) == obj2->get_X()[i]);
  }
  
  for (int i = 1; i <= 20; i++)
  {
    obj1->remove(i, 10.0);
    obj2->remove(i, 10.0);
  }
  assert(obj1->get_size() == obj2->get_size());
  assert(obj1->get_amount() == obj2->get_amount());
  for (int i = 0; i < (int)obj1->get_size(); i++)
  {
    assert(obj1->get(i+1) == obj2->get(i+1));
    assert(obj1->get_X()[i] == obj2->get_X()[i]);
    assert(obj1->get(i+1) == obj1->get_X()[i]);
    assert(obj2->get(i+1) == obj2->get_X()[i]);
  }
  
  size_t size = obj1->get_size();
  obj1->increase_size(1000);
  obj2->increase_size(1000);
  obj1->decrease_size(size);
  obj2->decrease_size(size);
  assert(obj1->get_size() == obj2->get_size());
  assert(obj1->get_amount() == obj2->get_amount());
  for (int i = 0; i < (int)obj1->get_size(); i++)
  {
    assert(obj1->get(i+1) == obj2->get(i+1));
    assert(obj1->get_X()[i] == obj2->get_X()[i]);
    assert(obj1->get(i+1) == obj1->get_X()[i]);
    assert(obj2->get(i+1) == obj2->get_X()[i]);
  }
  
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
void UnitaryTests::Prng_isEqualTo( Prng* prng1, Prng* prng2 )
{
  size_t test_size = 10000;
  
  /* 1) test uniform law */
  for (size_t i = 0; i < test_size; i++)
  {
    assert(prng1->uniform() == prng2->uniform());
  }
  /* 2) test integer uniform law */
  for (size_t i = 0; i < test_size; i++)
  {
    assert(prng1->uniform(-10000, 10000) == prng2->uniform(-10000, 10000));
  }
  /* 3) test bernouilli law */
  for (size_t i = 0; i < test_size; i++)
  {
    assert(prng1->bernouilli(0.5) == prng2->bernouilli(0.5));
  }
  /* 4) test binomial law */
  for (size_t i = 0; i < test_size; i++)
  {
    assert(prng1->binomial(1000, 0.5) == prng2->binomial(1000, 0.5));
  }
  /* 5) test multinomial law */
  for (size_t i = 0; i < test_size; i++)
  {
    unsigned int draws1[10];
    unsigned int draws2[10];
    double probas[10] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    int N = 1000, K = 10;
    prng1->multinomial(draws1, probas, N, K);
    prng2->multinomial(draws2, probas, N, K);
    for (size_t j = 0; j < 10; j++)
    {
      assert(draws1[j] == draws2[j]);
    }
  }
  /* 6) test gaussian law */
  for (size_t i = 0; i < test_size; i++)
  {
    assert(prng1->gaussian(0, 1) == prng2->gaussian(0, 1));
  }
  /* 7) test exponential law */
  for (size_t i = 0; i < test_size; i++)
  {
    assert(prng1->exponential(10.0) == prng2->exponential(10.0));
  }
  /* 8) test exponential law */
  for (size_t i = 0; i < test_size; i++)
  {
    double probas[10] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double sum = 5.0;
    int N = 10;
    assert(prng1->roulette_wheel(probas, sum, N) == prng2->roulette_wheel(probas, sum, N));
    
    (void)probas;
    (void)sum;
    (void)N;
  }
}

/**
 * \brief    modify Parameters attributes
 * \details  --
 * \param    Parameters* obj
 * \return   \e void
 */
void UnitaryTests::modify_Parameters( Parameters* obj )
{
  variable_range range;
  range.min = 1;
  range.max = 20;
  
  /*------------------------------------------------------------------ prng seed */
  
  obj->set_seed(101010);
  
  /*------------------------------------------------------------------ parallel computing */
  
  obj->set_parallel_computing(true);
  
  /*------------------------------------------------------------------ simulation schemes */
  
  obj->set_energy_costs_scheme(true);
  obj->set_membrane_permeability_scheme(true);
  obj->set_metabolic_inheritance_scheme(true);
  obj->set_enzymatic_inheritance_scheme(true);
  obj->set_co_enzyme_activity_scheme(true);
  obj->set_score_scheme(ESSENTIAL_METABOLITES_COMBINATORIAL_CONTRIBUTION);
  obj->set_selection_threshold(0.1234);
  
  /*------------------------------------------------------------------ space */
  
  //obj->set_width(101);
  //obj->set_height(101);
  
  /*------------------------------------------------------------------ output */
  
  obj->set_simulation_backup_step(101010);
  obj->set_figures_generation_step(101010);
  
  /*------------------------------------------------------------------ genome */
  
  obj->set_metabolite_tag_initial_range(&range);
  obj->set_binding_site_tag_initial_range(&range);
  obj->set_co_enzyme_tag_initial_range(&range);
  obj->set_transcription_factor_tag_initial_range(&range);
  
  obj->set_transcription_factor_binding_window(1234);
  
  obj->set_initial_number_of_NC_units(101);
  obj->set_initial_number_of_E_units(101);
  obj->set_initial_number_of_TF_units(101);
  obj->set_initial_number_of_BS_units(101);
  obj->set_initial_number_of_P_units(101);
  
  obj->set_point_mutation_rate(0.1234);
  obj->set_duplication_rate(0.1234);
  obj->set_deletion_rate(0.1234);
  obj->set_translocation_rate(0.1234);
  obj->set_inversion_rate(0.1234);
  obj->set_transition_rate(0.1234);
  obj->set_breakpoint_rate(0.1234);
  
  obj->set_substrate_tag_mutation_size(101);
  obj->set_product_tag_mutation_size(101);
  obj->set_kcat_mutation_size(0.1234);
  obj->set_kcat_km_ratio_mutation_size(0.1234);
  obj->set_binding_site_tag_mutation_size(101);
  obj->set_co_enzyme_tag_mutation_size(101);
  obj->set_transcription_factor_tag_mutation_size(101);
  obj->set_basal_expression_level_mutation_size(0.1234);
  
  /*------------------------------------------------------------------ genetic regulation network */
  
  obj->set_genetic_regulation_network_timestep(0.1234);
  obj->set_hill_function_theta(0.1234);
  obj->set_hill_function_n(0.1234);
  obj->set_protein_degradation_rate(0.1234);
  
  /*------------------------------------------------------------------ metabolic network */
  
  obj->set_metabolism_timestep(0.1234);
  
  obj->set_essential_metabolites_toxicity_threshold(0.1234);
  obj->set_non_essential_metabolites_toxicity_threshold(0.1234);
  
  obj->set_initial_metabolites_amount_in_cells(0.1234);
  
  obj->set_energy_toxicity_threshold(0.1234);
  
  /*------------------------------------------------------------------ energy */
  
  obj->set_energy_transcription_cost(0.001234);
  obj->set_energy_degradation_cost(0.001234);
  obj->set_energy_enzymatic_cost(0.001234);
  obj->set_energy_pumping_cost(0.001234);
  
  obj->set_energy_dissipation_rate(0.1234);
  
  obj->set_energy_toxicity_threshold(12.34);
  
  /*------------------------------------------------------------------ cell */
  
  obj->set_membrane_permeability(0.1234);
  
  /*------------------------------------------------------------------ population */
  
  obj->set_death_probability(0.1234);
  obj->set_migration_rate(0.1234);
  obj->set_hgt_rate(0.1234);
  
  /*------------------------------------------------------------------ environment */
  
  environment_properties properties;
  properties.number_of_init_cycles   = 12;
  properties.species_tag_range       = range;
  properties.concentration_range     = range;
  properties.number_of_species_range = range;
  properties.interaction_scheme      = NO_INTERACTION;
  properties.renewal_scheme          = CLEAR_MATTER;
  properties.variation_scheme        = RANDOM_SCHEME;
  properties.localization_scheme     = SPOT_LOCALIZATION;
  properties.metabolic_scheme        = MULTIPLE_METABOLITES;
  properties.introduction_rate       = 0.1234;
  properties.diffusion_coefficient   = 0.1234;
  properties.degradation_rate        = 0.1234;
  obj->set_environment_properties(&properties);
}
