
/**
 * \file      TrophicNetwork.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      14-10-2015
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     TrophicNetwork class declaration
 */

/****************************************************************************
 * Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * E-mail: charles.rocabert@gmail.com
 * Web: http://www.evoevo.eu/
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

#ifndef __EVOEVO__TrophicNetwork__
#define __EVOEVO__TrophicNetwork__

#include <iostream>
#include <unordered_map>
#include <zlib.h>
#include <assert.h>

#include "Macros.h"
#include "Structs.h"
#include "Enums.h"
#include "Parameters.h"
#include "TrophicGroup.h"
#include "Population.h"
#include "Environment.h"


class TrophicNetwork
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  TrophicNetwork( void ) = delete;
  TrophicNetwork( Parameters* parameters, Population* population, Environment* environment );
  TrophicNetwork( Parameters* parameters, Population* population, Environment* environment, gzFile backup_file );
  TrophicNetwork( const TrophicNetwork& trophic_network ) = delete;
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~TrophicNetwork( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  
  /*------------------------------------------------------------------ Trophic network attributes */
  
  inline unsigned long long int get_current_id( void ) const;
  inline size_t                 get_number_of_groups( void ) const;
  inline TrophicGroup*          get_group( unsigned long long int identifier );
  inline TrophicGroup*          get_first_group( void );
  inline TrophicGroup*          get_next_group( void );
  
  /*------------------------------------------------------------------ Trophic network statistics */
  
  inline size_t get_nb_level_0_groups( void ) const;
  inline size_t get_nb_level_1_groups( void ) const;
  inline size_t get_nb_level_2_groups( void ) const;
  inline size_t get_nb_no_level_groups( void ) const;
  
  inline size_t get_nb_level_0_cells( void ) const;
  inline size_t get_nb_level_1_cells( void ) const;
  inline size_t get_nb_level_2_cells( void ) const;
  inline size_t get_nb_no_level_cells( void ) const;
  
  inline size_t get_nb_group_appearances( void ) const;
  inline size_t get_nb_group_extinctions( void ) const;
  
  inline double get_mean_group_lifespan( void ) const;
  
  /*------------------------------------------------------------------ Identifier management */
  
  inline unsigned long long int get_new_id( void );
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  TrophicNetwork& operator=(const TrophicNetwork&) = delete;
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void save( gzFile backup_file );
  void initialize_trophic_network( void );
  void load_population( void );
  void write_trophic_network( std::string node_filename, std::string edge_filename );
  void write_trophic_network( double time, std::ofstream& file );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  bool groupExists( std::string trophic_profile, unsigned long long int &group_id );
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  
  /*------------------------------------------------------------------ Trophic network attributes */
  
  Parameters*                                                         _parameters;  /*!< Simulation parameters      */
  Population*                                                         _population;  /*!< Population                 */
  Environment*                                                        _environment; /*!< Environment                */
  unsigned long long int                                              _current_id;  /*!< Current group id           */
  std::unordered_map<unsigned long long int, TrophicGroup*>           _group_map;   /*!< Trophic group map          */
  std::unordered_map<unsigned long long int, TrophicGroup*>::iterator _iterator;    /*!< Trophic group map iterator */
  
  /*------------------------------------------------------------------ Trophic network statistics */
  
  size_t _nb_level_0_groups;  /*!< Number of groups of LEVEL 0 level  */
  size_t _nb_level_1_groups;  /*!< Number of groups of LEVEL 1 level  */
  size_t _nb_level_2_groups;  /*!< Number of groups of LEVEL 2 level  */
  size_t _nb_no_level_groups; /*!< Number of groups of NO LEVEL level */
  
  size_t _nb_level_0_cells;  /*!< Number of cells of LEVEL 0 level  */
  size_t _nb_level_1_cells;  /*!< Number of cells of LEVEL 1 level  */
  size_t _nb_level_2_cells;  /*!< Number of cells of LEVEL 2 level  */
  size_t _nb_no_level_cells; /*!< Number of cells of NO LEVEL level */
  
  size_t _nb_group_appearances; /*!< Number of group appearances per simulation timestep */
  size_t _nb_group_extinctions; /*!< Number of group extinctions per simulation timestep */
  
  double _mean_group_lifespan; /*!< Mean group lifespan */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/
/**
 * \brief    Get current group identifier
 * \details  --
 * \param    void
 * \return   \e unsigned long long int
 */
inline unsigned long long int TrophicNetwork::get_current_id( void ) const
{
  return _current_id;
}

/**
 * \brief    Get the number of groups of the network
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_number_of_groups( void ) const
{
  return _group_map.size();
}

/**
 * \brief    Get the node by its identifier
 * \details  Return NULL if the node do not exist
 * \param    unsigned long long int identifier
 * \return   \e TrophicGroup*
 */
inline TrophicGroup* TrophicNetwork::get_group( unsigned long long int identifier )
{
  if (_group_map.find(identifier) != _group_map.end())
  {
    return _group_map[identifier];
  }
  return NULL;
}

/**
 * \brief    Get the first group of the network
 * \details  Return NULL if the network is empty
 * \param    void
 * \return   \e TrophicGroup*
 */
inline TrophicGroup* TrophicNetwork::get_first_group( void )
{
  _iterator = _group_map.begin();
  if (_iterator != _group_map.end())
  {
    return _iterator->second;
  }
  return NULL;
}

/**
 * \brief    Get the next group
 * \details  Return NULL if the end of the network is reached
 * \param    void
 * \return   \e TrophicGroup*
 */
inline TrophicGroup* TrophicNetwork::get_next_group( void )
{
  _iterator++;
  if (_iterator != _group_map.end())
  {
    return _iterator->second;
  }
  return NULL;
}

/*------------------------------------------------------------------ Trophic network statistics */

/**
 * \brief    Get number of groups of level LEVEL 0
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_level_0_groups( void ) const
{
  return _nb_level_0_groups;
}

/**
 * \brief    Get number of groups of level LEVEL 1
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_level_1_groups( void ) const
{
  return _nb_level_1_groups;
}

/**
 * \brief    Get number of groups of level LEVEL 2
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_level_2_groups( void ) const
{
  return _nb_level_2_groups;
}

/**
 * \brief    Get number of groups of level NO LEVEL
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_no_level_groups( void ) const
{
  return _nb_no_level_groups;
}

/**
 * \brief    Get number of cells of level LEVEL 0
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_level_0_cells( void ) const
{
  return _nb_level_0_cells;
}

/**
 * \brief    Get number of cells of level LEVEL 1
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_level_1_cells( void ) const
{
  return _nb_level_1_cells;
}

/**
 * \brief    Get number of cells of level LEVEL 2
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_level_2_cells( void ) const
{
  return _nb_level_2_cells;
}

/**
 * \brief    Get number of cells of level NO LEVEL
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_no_level_cells( void ) const
{
  return _nb_no_level_cells;
}

/**
 * \brief    Get number of group appearances
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_group_appearances( void ) const
{
  return _nb_group_appearances;
}

/**
 * \brief    Get number of group extinctions
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_group_extinctions( void ) const
{
  return _nb_group_extinctions;
}

/**
 * \brief    Get mean group lifespan
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline double TrophicNetwork::get_mean_group_lifespan( void ) const
{
  return _mean_group_lifespan;
}

/*------------------------------------------------------------------ Identifier management */

/**
 * \brief    Get a new identifier
 * \details  --
 * \param    void
 * \return   \e unsigned long long int
 */
inline unsigned long long int TrophicNetwork::get_new_id( void )
{
  unsigned long long int new_id = _current_id;
  _current_id++;
  return new_id;
}

/*----------------------------
 * SETTERS
 *----------------------------*/


#endif /* defined(__EVOEVO__TrophicNetwork__) */
