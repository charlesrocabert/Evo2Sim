
/**
 * \file      IntegratedTests.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     IntegratedTests class declaration
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

#ifndef __EVOEVO__IntegratedTests__
#define __EVOEVO__IntegratedTests__

#include <iostream>
#include <cstring>
#include <sys/stat.h>
#include <assert.h>

#include "Macros.h"
#include "Enums.h"
#include "Structs.h"
#include "Parameters.h"
#include "MutationVector.h"
#include "MutationEvent.h"
#include "ReplicationReport.h"
#include "Genome.h"
#include "InheritedProteins.h"
#include "ODE.h"
#include "SpeciesList.h"
#include "Cell.h"
#include "Population.h"
#include "Environment.h"
#include "TrophicGroup.h"
#include "TrophicNetwork.h"
#include "Simulation.h"
#include "Node.h"
#include "Tree.h"
#include "Statistics.h"
#include "Prng.h"


class IntegratedTests
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  IntegratedTests( void ) = delete;
  IntegratedTests( Parameters* _parameters1, Parameters* _parameters2 );
  IntegratedTests( const IntegratedTests& test ) = delete;
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~IntegratedTests( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  IntegratedTests& operator=(const IntegratedTests&) = delete;
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void        run_integrated_tests( size_t number_of_tests, size_t number_of_steps, bool random_seed, bool random_parameters );
  Simulation* create_simulation( Parameters* parameters, unsigned long int seed );
  void        test_simulation_creation( unsigned long int seed );
  bool        test_simulation_execution( size_t backup_time, size_t simulation_time );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  void Parameters_isEqualTo( Parameters* obj1, Parameters* obj2 );
  void genetic_unit_isEqualTo( genetic_unit* obj1, genetic_unit* obj2 );
  void MutationVector_isEqualTo( MutationVector* obj1, MutationVector* obj2 );
  void MutationEvent_isEqualTo( MutationEvent* obj1, MutationEvent* obj2 );
  void ReplicationReport_isEqualTo( ReplicationReport* obj1, ReplicationReport* obj2 );
  void Genome_isEqualTo( Genome* obj1, Genome* obj2 );
  void InheritedProteins_isEqualTo( InheritedProteins* obj1, InheritedProteins* obj2 );
  void reaction_list_isEqualTo( reaction_list* obj1, reaction_list* obj2 );
  void ODE_isEqualTo( ODE* obj1, ODE* obj2 );
  void SpeciesList_isEqualTo( SpeciesList* obj1, SpeciesList* obj2 );
  void Cell_isEqualTo( Cell* obj1, Cell* obj2 );
  void Population_isEqualTo( Population* obj1, Population* obj2 );
  void Environment_isEqualTo( Environment* obj1, Environment* obj2 );
  void TrophicGroup_isEqualTo( TrophicGroup* obj1, TrophicGroup* obj2 );
  void TrophicNetwork_isEqualTo( TrophicNetwork* obj1, TrophicNetwork* obj2 );
  void Node_isEqualTo( Node* obj1, Node* obj2 );
  void Tree_isEqualTo( Tree* obj1, Tree* obj2 );
  void Statistics_isEqualTo( Statistics* obj1, Statistics* obj2 );
  void Simulation_isEqualTo( Simulation* obj1, Simulation* obj2 );
  void Prng_isEqualTo( Prng* obj1, Prng* obj2 );
  void generate_random_parameters( Parameters* parameters, Prng* prng );
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  Parameters* _parameters1; /*!< Parameters for the first simulation  */
  Parameters* _parameters2; /*!< Parameters for the second simulation */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/*----------------------------
 * SETTERS
 *----------------------------*/


#endif /* defined(__EVOEVO__IntegratedTests__) */
