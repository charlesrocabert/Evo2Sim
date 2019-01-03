
/**
 * \file      Simulation.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Simulation class declaration
 */

/****************************************************************************
 * Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon
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

#ifndef __Evo2Sim__Simulation__
#define __Evo2Sim__Simulation__

#include <iostream>
#include <sstream>
#include <vector>
#include <cstring>
#include <cmath>
#include <zlib.h>
#include <string.h>
#include <assert.h>
#include <tbb/tbb.h>

#include "Macros.h"
#include "Enums.h"
#include "Structs.h"
#include "Parameters.h"
#include "Population.h"
#include "Environment.h"
#include "TrophicNetwork.h"
#include "Tree.h"
#include "Statistics.h"
#include "Prng.h"


class Simulation
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  Simulation( void ) = delete;
  Simulation( Parameters* parameters );
  Simulation( Parameters* parameters, size_t backup_time, bool clean_statistic_files );
  Simulation( const Simulation& simulation ) = delete;
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~Simulation( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  inline Parameters*     get_parameters( void );
  inline Population*     get_population( void );
  inline Environment*    get_environment( void );
  inline TrophicNetwork* get_trophic_network( void );
  inline Tree*           get_lineage_tree( void );
  inline Tree*           get_phylogenetic_tree( void );
  inline Statistics*     get_statistics( void );
  inline double          get_min_score( void ) const;
  inline double          get_max_score( void ) const;
  inline size_t          get_width( void ) const;
  inline size_t          get_height( void ) const;
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  Simulation& operator=(const Simulation&) = delete;
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void initialize( simulation_state state );
  void initialize_from_evolved_population( void );
  void update( void );
  void load_simulation( size_t backup_time );
  void save_simulation( void );
  void add_random_genetic_units( size_t cell_pos, size_t NC_type, size_t E_type, size_t TF_type, size_t BS_type, size_t P_type, bool shuffle );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  void initialize_environment( void );
  void initialize_population( simulation_state state );
  void initialize_population_from_evolved_population( void );
  void draw_random_genetic_unit( genetic_unit& unit, genetic_unit_type type );
  void update_environment( void );
  void update_population( void );
  void update_trophic_network( void );
  void update_trees( void );
  
  void save_parameters( void );
  void save_population( void );
  void save_environment( void );
  void save_trophic_network( void );
  void save_lineage_tree( void );
  void save_phylogenetic_tree( void );
  
  void mix( void );
  
  void kill_cell( size_t i );
  void mutate_cell( size_t i );
  void update_cell( size_t i );
  void divide_cell( size_t i, size_t child_position );
  bool find_empty_cases( size_t i, size_t& child_position );
  
#ifdef DEBUG
  void test_lineage_tree_structure( void );
  void test_phylogenetic_tree_structure( void );
#endif
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  Parameters*     _parameters;        /*!< Simulation parameters           */
  Population*     _population;        /*!< Population                      */
  Environment*    _environment;       /*!< Environment                     */
  TrophicNetwork* _trophic_network;   /*!< Trophic network                 */
  Tree*           _lineage_tree;      /*!< Lineage tree                    */
  Tree*           _phylogenetic_tree; /*!< Phylogenetic tree               */
  Statistics*     _statistics;        /*!< Statistics                      */
  double          _min_score;         /*!< Minimum score in the population */
  double          _max_score;         /*!< Maximum score in the population */
  size_t          _width;             /*!< Grid width                      */
  size_t          _height;            /*!< Grid height                     */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/**
 * \brief    Get simulation parameters
 * \details  --
 * \param    void
 * \return   \e Parameters*
 */
inline Parameters* Simulation::get_parameters( void )
{
  return _parameters;
}

/**
 * \brief    Get population
 * \details  --
 * \param    void
 * \return   \e Population*
 */
inline Population* Simulation::get_population( void )
{
  return _population;
}

/**
 * \brief    Get environment
 * \details  --
 * \param    void
 * \return   \e Environment*
 */
inline Environment* Simulation::get_environment( void )
{
  return _environment;
}

/**
 * \brief    Get trophic network
 * \details  --
 * \param    void
 * \return   \e TrophicNetwork*
 */
inline TrophicNetwork* Simulation::get_trophic_network( void )
{
  return _trophic_network;
}

/**
 * \brief    Get lineage tree
 * \details  --
 * \param    void
 * \return   \e Tree*
 */
inline Tree* Simulation::get_lineage_tree( void )
{
  return _lineage_tree;
}

/**
 * \brief    Get phylogenetic tree
 * \details  --
 * \param    void
 * \return   \e Tree*
 */
inline Tree* Simulation::get_phylogenetic_tree( void )
{
  return _phylogenetic_tree;
}

/**
 * \brief    Get statistics
 * \details  --
 * \param    Statistics*
 * \return   \e void
 */
inline Statistics* Simulation::get_statistics( void )
{
  return _statistics;
}

/**
 * \brief    Get minimum score in the population
 * \details  --
 * \param    void
 * \return   \e int
 */
inline double Simulation::get_min_score( void ) const
{
  return _min_score;
}

/**
 * \brief    Get maximum score in the population
 * \details  --
 * \param    void
 * \return   \e int
 */
inline double Simulation::get_max_score( void ) const
{
  return _max_score;
}

/**
 * \brief    Get the grid width
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Simulation::get_width( void ) const
{
  return _width;
}

/**
 * \brief    Get the grid height
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Simulation::get_height( void ) const
{
  return _height;
}

/*----------------------------
 * SETTERS
 *----------------------------*/

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/


#endif /* defined(__Evo2Sim__Simulation__) */
