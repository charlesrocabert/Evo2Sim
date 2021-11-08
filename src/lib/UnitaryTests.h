
/**
 * \file      UnitaryTests.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     UnitaryTests class declaration
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

#ifndef __Evo2Sim__UnitaryTests__
#define __Evo2Sim__UnitaryTests__

#include <iostream>
#include <cstring>
#include <sys/stat.h>
#include <assert.h>

#include "Macros.h"
#include "Enums.h"
#include "Structs.h"
#include "Prng.h"
#include "Parameters.h"
#include "MutationVector.h"
#include "MutationEvent.h"
#include "ReplicationReport.h"
#include "Genome.h"
#include "InheritedProteins.h"
#include "SpeciesList.h"


class UnitaryTests
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  UnitaryTests( void ) = delete;
  UnitaryTests( Parameters* parameters );
  UnitaryTests( const UnitaryTests& test ) = delete;
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~UnitaryTests( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  UnitaryTests& operator=(const UnitaryTests&) = delete;
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void run_unitary_tests( std::string parameters_filename );
  
  void test_Parameters_class( std::string filename );
  void test_MutationVector_class( void );
  void test_MutationEvent_class( void );
  void test_ReplicationReport_class( void );
  void test_Genome_class( void );
  void test_InheritedProteins_class( void );
  void test_SpeciesList_class( void );
  
  void test_Prng_class( void );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  void initialize_mutation_rates( void );
  void initialize_genetic_unit( genetic_unit* obj );
  
  void Parameters_isEqualTo( Parameters* obj1, Parameters* obj2 );
  void genetic_unit_isEqualTo( genetic_unit* obj1, genetic_unit* obj2 );
  void MutationVector_isEqualTo( MutationVector* obj1, MutationVector* obj2 );
  void MutationEvent_isEqualTo( MutationEvent* obj1, MutationEvent* obj2 );
  void ReplicationReport_isEqualTo( ReplicationReport* obj1, ReplicationReport* obj2 );
  void Genome_isEqualTo( Genome* obj1, Genome* obj2 );
  void InheritedProteins_isEqualTo( InheritedProteins* obj1, InheritedProteins* obj2 );
  void SpeciesList_isEqualTo( SpeciesList* obj1, SpeciesList* obj2 );
  
  void Prng_isEqualTo( Prng* prng1, Prng* prng2 );
  
  void modify_Parameters( Parameters* obj );
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  Parameters* _parameters;
  double*     _mutation_rates;
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/*----------------------------
 * SETTERS
 *----------------------------*/


#endif /* defined(__Evo2Sim__UnitaryTests__) */
