
/**
 * \file      Node.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      10-09-2015
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Node class definition
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

#include "Node.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  Create a basic node
 * \param    unsigned long long int identifier
 * \return   \e void
 */
Node::Node( unsigned long long int identifier )
{
  _identifier         = identifier;
  _alive_cell         = NULL;
  _replication_report = NULL;
  _parent             = NULL;
  _children.clear();
  _node_class = MASTER_ROOT;
  _node_state = DEAD;
  _tagged     = false;
}

/**
 * \brief    Alive node constructor
 * \details  Create a new alive node with its related cell. New node must exclusively be linked to alive cells
 * \param    unsigned long long int identifier
 * \param    Cell* alive_cell
 * \return   \e void
 */
Node::Node( unsigned long long int identifier, Cell* alive_cell )
{
  assert(alive_cell->isAlive());
  _identifier         = identifier;
  _alive_cell         = alive_cell;
  _replication_report = alive_cell->get_replication_report();
  _parent             = NULL;
  _children.clear();
  _node_class = NORMAL;
  _node_state = ALIVE;
  _tagged     = false;
}

/**
 * \brief    Constructor from backup file
 * \details  Load Node class from backup file. Relationships are not loaded, but are managed in Tree class
 * \param    gzFile backup_file
 * \return   \e void
 */
Node::Node( gzFile backup_file )
{
  gzread( backup_file, &_identifier, sizeof(_identifier) );
  gzread( backup_file, &_node_class, sizeof(_node_class) );
  gzread( backup_file, &_node_state, sizeof(_node_state) );
  gzread( backup_file, &_tagged,     sizeof(_tagged) );
  _alive_cell         = NULL;
  _replication_report = NULL;
  if (_node_class != MASTER_ROOT && _node_state == DEAD)
  {
    _replication_report = new ReplicationReport(backup_file);
  }
  _parent = NULL;
  _children.clear();
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
Node::~Node( void )
{
  _alive_cell = NULL;
  if (_node_state == DEAD)
  {
    delete _replication_report;
  }
  _replication_report = NULL;
  _children.clear();
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Save in backup file
 * \details  Relationships are not saved, but are managed in Tree class
 * \param    gzFile backup_file
 * \return   \e void
 */
void Node::save( gzFile backup_file )
{
  gzwrite( backup_file, &_identifier, sizeof(_identifier) );
  gzwrite( backup_file, &_node_class, sizeof(_node_class) );
  gzwrite( backup_file, &_node_state, sizeof(_node_state) );
  gzwrite( backup_file, &_tagged, sizeof(_tagged) );
  if (_node_class != MASTER_ROOT && _node_state == DEAD)
  {
    assert(_replication_report != NULL);
    _replication_report->save(backup_file);
  }
}

/**
 * \brief    Add a child
 * \details  --
 * \param    Node* node
 * \return   \e void
 */
void Node::add_child( Node* node )
{
  for (size_t i = 0; i < _children.size(); i++)
  {
    assert(node->get_id() != _children[i]->get_id());
  }
  _children.push_back(node);
}

/**
 * \brief    Remove a child
 * \details  --
 * \param    Node* node
 * \return   \e void
 */
void Node::remove_child( Node* node )
{
  int pos = -1;
  for (size_t i = 0; i < _children.size(); i++)
  {
    if (node->get_id() == _children[i]->get_id() && pos == -1)
    {
      pos = (int)i;
    }
    else if (node->get_id() == _children[i]->get_id() && pos >= 0)
    {
      printf("Error in Node::remove_child(): multiple occurences of a node in children list. Exit.\n");
      exit(EXIT_FAILURE);
    }
  }
  if (pos == -1)
  {
    printf("Error in Node::remove_child(): node to remove do not exist. Exit.\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    _children.erase(_children.begin()+pos);
  }
}

/**
 * \brief    Replace this child by its own children
 * \details  --
 * \param    Node* node_to_remove
 * \return   \e void
 */
void Node::replace_children( Node* child_to_remove )
{
  /*---------------------------------------*/
  /* 1) remove the node from children list */
  /*---------------------------------------*/
  remove_child(child_to_remove);
  
  /*---------------------------------------*/
  /* 2) add children to the children list  */
  /*---------------------------------------*/
  for (size_t i = 0; i < child_to_remove->get_number_of_children(); i++)
  {
    add_child(child_to_remove->get_child(i));
  }
}

/**
 * \brief    Tag the lineage of this node
 * \details  --
 * \param    void
 * \return   \e void
 */
void Node::tag_lineage( void )
{
  _tagged = true;
  Node* node = _parent;
  while (node != NULL)
  {
    node->tag();
    node = node->get_parent();
    if (node != NULL)
    {
      if (node->isTagged())
      {
        node = NULL;
      }
    }
  }
}

/**
 * \brief    Untag the lineage of this node
 * \details  --
 * \param    void
 * \return   \e void
 */
void Node::untag_lineage( void )
{
  _tagged = false;
  Node* node = _parent;
  while (node != NULL)
  {
    node->untag();
    node = node->get_parent();
    if (node != NULL)
    {
      if (!node->isTagged())
      {
        node = NULL;
      }
    }
  }
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/
