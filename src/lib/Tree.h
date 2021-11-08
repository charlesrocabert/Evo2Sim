
/**
 * \file      Tree.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      10-09-2015
 * \copyright Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Tree class declaration
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

#ifndef __Evo2Sim__Tree__
#define __Evo2Sim__Tree__

#include <iostream>
#include <vector>
#include <unordered_map>
#include <zlib.h>
#include <cstring>
#include <cmath>

#include "Macros.h"
#include "Structs.h"
#include "Enums.h"
#include "Parameters.h"
#include "ReplicationReport.h"
#include "Population.h"
#include "Node.h"
#include "TrophicGroup.h"


class Tree
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  Tree( void ) = delete;
  Tree( Parameters* parameters );
  Tree( Parameters* parameters, Population* population, gzFile backup_file );
  Tree( const Tree& tree ) = delete;
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~Tree( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  inline unsigned long long int get_current_id( void ) const;
  inline size_t get_number_of_nodes( void ) const;
  inline Node*  get_node( unsigned long long int identifier );
  inline Node*  get_first_node( void );
  inline Node*  get_next_node( void );
  inline void   get_alive_nodes( std::vector<unsigned long long int>* alive_nodes );
  inline Node*  get_best_alive_node( void );
  inline Node*  get_common_ancestor( void );
  inline double get_common_ancestor_age( void );
  inline unsigned long long int get_node_id_by_alive_cell_id( unsigned long long int alive_cell_identifier );
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  Tree& operator=(const Tree&) = delete;
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void save( gzFile backup_file );
  void add_root( Cell* cell );
  void add_division( Cell* parent, Cell* child1, Cell* child2 );
  void freeze_node( unsigned long long int cell_identifier, size_t death_time );
  void delete_node( unsigned long long int node_identifier );
  void clean_cell_map( void );
  void prune();
  void shorten();
  
  void write_tree( std::string filename );
  void write_newick_tree( std::string filename );
  
  void write_lineage_statistics( std::string filename, unsigned long long int identifier );
  void write_phylogeny_statistics( std::string filename );
  void write_trophic_data( std::string filename );
  
  void compute_AI_score_on_SL( size_t backup_time );
  void compute_common_ancestor_SL_repartition( double& Sp_A, double& Sp_B );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  void inOrderNewick( Node* node, size_t parent_time, std::stringstream& output );
  void tag_tree();
  void untag_tree();
  void tag_offspring( Node* node, std::vector<Node*>* tagged_nodes );
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  Parameters*                                                 _parameters; /*!< Simulation parameters */
  unsigned long long int                                      _current_id; /*!< Current node id       */
  std::unordered_map<unsigned long long int, Node*>           _node_map;   /*!< Tree nodes map        */
  std::unordered_map<unsigned long long int, Node*>           _cell_map;   /*!< Cells map             */
  std::unordered_map<unsigned long long int, Node*>::iterator _iterator;   /*!< Tree map iterator     */
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/**
 * \brief    Get current node identifier
 * \details  --
 * \param    void
 * \return   \e unsigned long long int
 */
inline unsigned long long int Tree::get_current_id( void ) const
{
  return _current_id;
}

/**
 * \brief    Get the number of nodes of the tree
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Tree::get_number_of_nodes( void ) const
{
  return _node_map.size();
}

/**
 * \brief    Get the node by its identifier
 * \details  Return NULL if the node do not exist
 * \param    unsigned long long int identifier
 * \return   \e Node*
 */
inline Node* Tree::get_node( unsigned long long int identifier )
{
  if (_node_map.find(identifier) != _node_map.end())
  {
    return _node_map[identifier];
  }
  return NULL;
}

/**
 * \brief    Get the first node of the map
 * \details  Return NULL if the tree is empty
 * \param    void
 * \return   \e Node*
 */
inline Node* Tree::get_first_node( void )
{
  _iterator = _node_map.begin();
  if (_iterator != _node_map.end())
  {
    return _iterator->second;
  }
  return NULL;
}

/**
 * \brief    Get the next node
 * \details  Return NULL if the end of the tree is reached
 * \param    void
 * \return   \e Node*
 */
inline Node* Tree::get_next_node( void )
{
  _iterator++;
  if (_iterator != _node_map.end())
  {
    return _iterator->second;
  }
  return NULL;
}

/**
 * \brief    Get alive node identifiers list
 * \details  --
 * \param    std::vector<unsigned long long int>* alive_nodes
 * \return   \e void
 */
inline void Tree::get_alive_nodes( std::vector<unsigned long long int>* alive_nodes )
{
  alive_nodes->clear();
  for (_iterator = _cell_map.begin(); _iterator != _cell_map.end(); ++_iterator)
  {
    alive_nodes->push_back(_iterator->second->get_id());
  }
}

/**
 * \brief    Get best alive node
 * \details  --
 * \param    void
 * \return   \e Node*
 */
inline Node* Tree::get_best_alive_node( void )
{
  double best_score = 0.0;
  Node*  best_node  = NULL;
  for (_iterator = _node_map.begin(); _iterator != _node_map.end(); ++_iterator)
  {
    if (_iterator->second->isAlive())
    {
      if (best_score < _iterator->second->get_alive_cell()->get_score())
      {
        best_score = _iterator->second->get_alive_cell()->get_score();
        best_node  = _iterator->second;
      }
    }
  }
  return best_node;
}

/**
 * \brief    Get the common ancestor
 * \details  Return NULL if the population is extincted or if the tree is multirooted
 * \param    void
 * \return   \e Node*
 */
inline Node* Tree::get_common_ancestor( void )
{
  Node* master_root = _node_map[0];
  if (master_root->get_number_of_children() == 1)
  {
    return master_root->get_child(0);
  }
  return NULL;
}

/**
 * \brief    Get the common ancestor age
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Tree::get_common_ancestor_age( void )
{
  Node* master_root = _node_map[0];
  
  /*-----------------------------------------*/
  /* 1) If the population extincted          */
  /*-----------------------------------------*/
  if(master_root->get_number_of_children() == 0)
  {
    return 0.0;
  }
  
  /*-----------------------------------------*/
  /* 2) If there is a unique common ancestor */
  /*-----------------------------------------*/
  else if (master_root->get_number_of_children() == 1)
  {
    if (master_root->get_child(0)->isAlive())
    {
      return (double)master_root->get_child(0)->get_alive_cell()->get_birth_time();
    }
    else
    {
      return (double)master_root->get_child(0)->get_replication_report()->get_birth_time();
    }
  }
  
  /*-------------------------------------------*/
  /* 3) If there are multiple common ancestors */
  /*-------------------------------------------*/
  else
  {
    double mean = 0.0;
    for (size_t i = 0; i < master_root->get_number_of_children(); i++)
    {
      if (master_root->get_child(0)->isAlive())
      {
        mean += (double)master_root->get_child(0)->get_alive_cell()->get_birth_time();
      }
      else
      {
        mean += (double)master_root->get_child(0)->get_replication_report()->get_birth_time();
      }
    }
    return mean/master_root->get_number_of_children();
  }
}

/**
 * \brief    Get an alive node identifier by its linked alive cell identifier
 * \details  --
 * \param    void
 * \return   \e double
 */
inline unsigned long long int Tree::get_node_id_by_alive_cell_id( unsigned long long int alive_cell_identifier )
{
  for (_iterator = _node_map.begin(); _iterator != _node_map.end(); ++_iterator)
  {
    if (_iterator->second->isAlive())
    {
      if (_iterator->second->get_alive_cell()->get_id() == alive_cell_identifier)
      {
        return _iterator->first;
      }
    }
  }
  return 0;
}

/*----------------------------
 * SETTERS
 *----------------------------*/


#endif /* defined(__Evo2Sim__Tree__) */
