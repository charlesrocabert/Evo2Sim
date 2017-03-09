
/**
 * \file      Node.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      10-09-2015
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Node class declaration
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

#ifndef __EVOEVO__Node__
#define __EVOEVO__Node__

#include <iostream>
#include <vector>
#include <zlib.h>
#include <assert.h>

#include "Macros.h"
#include "Structs.h"
#include "Enums.h"
#include "Parameters.h"
#include "ReplicationReport.h"
#include "Cell.h"


class Node
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  Node( void ) = delete;
  Node( unsigned long long int identifier );
  Node( unsigned long long int identifier, Cell* alive_cell );
  Node( gzFile backup_file );
  Node( const Node& node ) = delete;
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~Node( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  inline unsigned long long int get_id( void ) const;
  inline Cell*                  get_alive_cell( void );
  inline ReplicationReport*     get_replication_report( void );
  inline Node*                  get_parent( void );
  inline Node*                  get_child( size_t pos );
  inline size_t                 get_number_of_children( void ) const;
  inline node_class             get_node_class( void ) const;
  inline node_state             get_node_state( void ) const;
  inline bool                   isTagged( void ) const;
  inline bool                   isMasterRoot( void ) const;
  inline bool                   isRoot( void ) const;
  inline bool                   isNormal( void ) const;
  inline bool                   isDead( void ) const;
  inline bool                   isAlive( void ) const;
  inline bool                   isAncestor( unsigned long long int ancestor_id );
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  Node& operator=(const Node&) = delete;
  
  inline void set_alive_cell( Cell* cell );
  inline void set_replication_report( ReplicationReport* replication_report );
  inline void set_parent( Node* node );
  inline void tag( void );
  inline void untag( void );
  inline void set_master_root( void );
  inline void set_root( void );
  inline void set_normal( void );
  inline void set_dead( size_t death_time );
  inline void set_alive( void );
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void save( gzFile backup_file );
  void add_child( Node* node );
  void remove_child( Node* node );
  void replace_children( Node* child_to_remove );
  void tag_lineage( void );
  void untag_lineage( void );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  unsigned long long int _identifier;         /*!< Node identifier                          */
  Cell*                  _alive_cell;         /*!< Alive cell pointer                       */
  ReplicationReport*     _replication_report; /*!< Cell's replication report                */
  Node*                  _parent;             /*!< Parent of the node                       */
  std::vector<Node*>     _children;           /*!< Children of the node                     */
  node_class             _node_class;         /*!< Node class (master root, root or normal) */
  node_state             _node_state;         /*!< Node state (dead or alive)               */
  bool                   _tagged;             /*!< Indicates if the node is tagged          */
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/**
 * \brief    Get node identifier
 * \details  --
 * \param    void
 * \return   \e unsigned long long int
 */
inline unsigned long long int Node::get_id( void ) const
{
  return _identifier;
}

/**
 * \brief    Get the alive cell pointer
 * \details  Return the cell if it is alive, or NULL if the cell is dead.
 * \param    void
 * \return   \e Cell*
 */
inline Cell* Node::get_alive_cell( void )
{
  return _alive_cell;
}

/**
 * \brief    Get the replication report linked to the node
 * \details  --
 * \param    void
 * \return   \e ReplicationReport*
 */
inline ReplicationReport* Node::get_replication_report( void )
{
  return _replication_report;
}

/**
 * \brief    Get the parent node
 * \details  --
 * \param    void
 * \return   \e Node*
 */
inline Node* Node::get_parent( void )
{
  return _parent;
}

/**
 * \brief    Get the child at position 'pos'
 * \details  --
 * \param    void
 * \return   \e Node*
 */
inline Node* Node::get_child( size_t pos )
{
  assert(pos < _children.size());
  return _children[pos];
}

/**
 * \brief    Get the number of children
 * \details  --
 * \param    void
 * \return   \e Node*
 */
inline size_t Node::get_number_of_children( void ) const
{
  return _children.size();
}

/**
 * \brief    Get the node class
 * \details  --
 * \param    void
 * \return   \e node_class
 */
inline node_class Node::get_node_class( void ) const
{
  return _node_class;
}

/**
 * \brief    Get the node state
 * \details  --
 * \param    void
 * \return   \e node_state
 */
inline node_state Node::get_node_state( void ) const
{
  return _node_state;
}

/**
 * \brief    Check if the node is tagged or not
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Node::isTagged( void ) const
{
  return _tagged;
}

/**
 * \brief    Check if the node is the master root
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Node::isMasterRoot( void ) const
{
  return (_node_class == MASTER_ROOT);
}

/**
 * \brief    Check if the node is a root or not
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Node::isRoot( void ) const
{
  return (_node_class == ROOT);
}

/**
 * \brief    Check if the node is normal or not
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Node::isNormal( void ) const
{
  return (_node_class == NORMAL);
}

/**
 * \brief    Check if the node is dead
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Node::isDead( void ) const
{
  return (_node_state == DEAD);
}

/**
 * \brief    Check if the node is alive
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Node::isAlive( void ) const
{
  return (_node_state == ALIVE);
}

/**
 * \brief    Check if the given identifier is an ancestor
 * \details  --
 * \param    unsigned long long int ancestor_id
 * \return   \e bool
 */
inline bool Node::isAncestor( unsigned long long int ancestor_id )
{
  Node* node = get_parent();
  while (node != NULL)
  {
    if (node->get_id() == ancestor_id)
    {
      return true;
    }
    node = node->get_parent();
  }
  return false;
}

/*----------------------------
 * SETTERS
 *----------------------------*/

/**
 * \brief    Set the alive cell's pointer (and its replication report pointer)
 * \details  --
 * \param    Cell* cell
 * \return   \e void
 */
inline void Node::set_alive_cell( Cell* cell )
{
  assert(_node_state == ALIVE);
  _alive_cell         = cell;
  _replication_report = cell->get_replication_report();
}

/**
 * \brief    Set the replication report pointer
 * \details  Only alive node should be modified
 * \param    ReplicationReport* replication_report
 * \return   \e void
 */
inline void Node::set_replication_report( ReplicationReport* replication_report )
{
  assert(_node_state == ALIVE);
  _replication_report = replication_report;
}

/**
 * \brief    Add a parent
 * \details  --
 * \param    Node* node
 * \return   \e void
 */
inline void Node::set_parent( Node* node )
{
  _parent = node;
}

/**
 * \brief    Tag the node
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Node::tag( void )
{
  _tagged = true;
}

/**
 * \brief    Untag the node
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Node::untag( void )
{
  _tagged = false;
}

/**
 * \brief    Set the node class as master root
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Node::set_master_root( void )
{
  _identifier         = 0;
  _alive_cell         = NULL;
  _replication_report = NULL;
  _parent             = NULL;
  _children.clear();
  _node_class = MASTER_ROOT;
  _node_state = DEAD;
  _tagged     = false;
}

/**
 * \brief    Set the node class as root
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Node::set_root( void )
{
  _node_class = ROOT;
}

/**
 * \brief    Set the node class as normal
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Node::set_normal( void )
{
  _node_class = NORMAL;
}

/**
 * \brief    Set the node state dead
 * \details  --
 * \param    size_t death_time
 * \return   \e void
 */
inline void Node::set_dead( size_t death_time )
{
  /*----------------------------------------*/
  /* 1) Copy and set the replication report */
  /*----------------------------------------*/
  _replication_report = new ReplicationReport(*_replication_report);
  _replication_report->set_death_time(death_time);
  
  /*----------------------------------------*/
  /* 2) Set the alive cell pointer to NULL  */
  /*----------------------------------------*/
  _alive_cell = NULL;
  
  /*----------------------------------------*/
  /* 3) set the new node state              */
  /*----------------------------------------*/
  _node_state = DEAD;
}


#endif /* defined(__EVOEVO__Node__) */
