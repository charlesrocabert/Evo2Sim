
/**
 * \file      Tree.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      10-09-2015
 * \copyright Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Tree class definition
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

#include "Tree.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  The tree starts with one node called the master root
 * \param    Parameters* parameters
 * \return   \e void
 */
Tree::Tree( Parameters* parameters )
{
  /*------------------------------------*/
  /* 1) Set simulation parameters       */
  /*------------------------------------*/
  _parameters = parameters;
  
  /*------------------------------------*/
  /* 2) Set the current identifier to 0 */
  /*------------------------------------*/
  _current_id = 0;
  
  /*------------------------------------*/
  /* 3) Initialize maps and iterator    */
  /*------------------------------------*/
  _node_map.clear();
  _cell_map.clear();
  _iterator = _node_map.begin();
  
  /*------------------------------------*/
  /* 4) Create the master root node and */
  /*    add it to the tree              */
  /*------------------------------------*/
  Node* master_root = new Node(_current_id);
  master_root->set_master_root();
  _node_map[_current_id] = master_root;
}

/**
 * \brief    Constructor from backup file
 * \details  Load Tree class from backup file. First, load nodes. Then, load tree relationships (parent and children information), and load alive cells in nodes.
 * \param    Parameters* parameters
 * \param    Population* pop
 * \param    gzFile backup_file
 * \return   \e void
 */
Tree::Tree( Parameters* parameters, Population* pop, gzFile backup_file )
{
  _parameters = parameters;
  
  /*-----------------------------------------*/
  /* 1) Recover the current identifier       */
  /*-----------------------------------------*/
  gzread( backup_file, &_current_id, sizeof(_current_id) );
  
  /*-----------------------------------------*/
  /* 2) Recover nodes map                    */
  /*-----------------------------------------*/
  _node_map.clear();
  size_t n = 0;
  gzread( backup_file, &n, sizeof(n) );
  for (size_t i = 0; i < n; i++)
  {
    Node* node = new Node(backup_file);
    _node_map[node->get_id()] = node;
  }
  
  /*-----------------------------------------*/
  /* 3) Recover nodes map relationships      */
  /*-----------------------------------------*/
  size_t number_of_edges = 0;
  gzread( backup_file, &number_of_edges, sizeof(number_of_edges) );
  for (size_t i = 0; i < number_of_edges; i++)
  {
    unsigned long long int first  = 0;
    unsigned long long int second = 0;
    gzread( backup_file, &first, sizeof(first) );
    gzread( backup_file, &second, sizeof(second) );
    assert(_node_map.find(first) != _node_map.end());
    assert(_node_map.find(second) != _node_map.end());
    _node_map[first]->add_child(_node_map[second]);
    _node_map[second]->set_parent(_node_map[first]);
  }
  
  /*-----------------------------------------*/
  /* 4) Recover cells map                    */
  /*-----------------------------------------*/
  _cell_map.clear();
  gzread( backup_file, &n, sizeof(n) );
  for (size_t i = 0; i < n; i++)
  {
    unsigned long long int first  = 0;
    unsigned long long int second = 0;
    gzread( backup_file, &first, sizeof(first) );
    gzread( backup_file, &second, sizeof(second) );
    assert(_node_map.find(second) != _node_map.end());
    _cell_map[first] = _node_map[second];
  }
  
  /*-----------------------------------------*/
  /* 5) Recover linked cells from population */
  /*-----------------------------------------*/
  for (_iterator = _cell_map.begin(); _iterator != _cell_map.end(); ++_iterator)
  {
    Cell* cell = pop->get_cell_by_id(_iterator->first);
    assert(cell != NULL);
    _iterator->second->set_alive_cell(cell); /* set_alive_cell also sets the replication report */
  }
  
  /*-----------------------------------------*/
  /* 6) Initialize iterator position         */
  /*-----------------------------------------*/
  _iterator = _node_map.begin();
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
Tree::~Tree( void )
{
  _iterator = _node_map.begin();
  for (_iterator = _node_map.begin(); _iterator != _node_map.end(); ++_iterator)
  {
    delete _iterator->second;
    _iterator->second = NULL;
  }
  _node_map.clear();
  _cell_map.clear();
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Save in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void Tree::save( gzFile backup_file )
{
  /*------------------------------------------------*/
  /* 1) Save the current identifier                 */
  /*------------------------------------------------*/
  gzwrite( backup_file, &_current_id, sizeof(_current_id) );
  
  /*------------------------------------------------*/
  /* 2) Save nodes map in the backup file           */
  /*------------------------------------------------*/
  size_t number_of_edges = 0;
  size_t n               = _node_map.size();
  gzwrite( backup_file, &n, sizeof(n) );
  for (_iterator = _node_map.begin(); _iterator != _node_map.end(); ++_iterator)
  {
    _iterator->second->save(backup_file);
    number_of_edges += _iterator->second->get_number_of_children();
  }
  
  /*------------------------------------------------*/
  /* 3) Save nodes map relationships in backup file */
  /*------------------------------------------------*/
  gzwrite( backup_file, &number_of_edges, sizeof(number_of_edges) );
  for (_iterator = _node_map.begin(); _iterator != _node_map.end(); ++_iterator)
  {
    unsigned long long int first = _iterator->first;
    for (size_t i = 0; i < _iterator->second->get_number_of_children(); i++)
    {
      unsigned long long int second = _iterator->second->get_child(i)->get_id();
      gzwrite( backup_file, &first, sizeof(first) );
      gzwrite( backup_file, &second, sizeof(second) );
    }
  }
  
  /*------------------------------------------------*/
  /* 4) Save cells map in the backup file           */
  /*------------------------------------------------*/
  n = _cell_map.size();
  gzwrite( backup_file, &n, sizeof(n) );
  
  for (_iterator = _cell_map.begin(); _iterator != _cell_map.end(); ++_iterator)
  {
    unsigned long long int first  = _iterator->first;
    unsigned long long int second = _iterator->second->get_id();
    gzwrite( backup_file, &first, sizeof(first) );
    gzwrite( backup_file, &second, sizeof(second) );
  }
}

/**
 * \brief    Add a root to the tree
 * \details  --
 * \param    Cell* cell
 * \return   \e void
 */
void Tree::add_root( Cell* cell )
{
  /*-----------------------------*/
  /* 1) Get the master root      */
  /*-----------------------------*/
  Node* master_root = _node_map[0];
  
  /*-----------------------------*/
  /* 2) Create the node          */
  /*-----------------------------*/
  _current_id++;
  Node* node = new Node(_current_id, cell);
  
  /*-----------------------------*/
  /* 3) Update nodes attributes  */
  /*-----------------------------*/
  node->set_root();
  node->set_parent(master_root);
  master_root->add_child(node);
  
  /*-----------------------------*/
  /* 4) Add the node to the tree */
  /*-----------------------------*/
  assert(_node_map.find(node->get_id()) == _node_map.end());
  _node_map[node->get_id()] = node;
  
  /*-----------------------------*/
  /* 5) Update the cells map     */
  /*-----------------------------*/
  _cell_map[cell->get_id()] = node;
}

/**
 * \brief    Add a division to the tree
 * \details  --
 * \param    Cell* parent
 * \param    Cell* child1
 * \param    Cell* child2
 * \return   \e void
 */
void Tree::add_division( Cell* parent, Cell* child1, Cell* child2 )
{
  /*----------------------------*/
  /* 1) Get parental node       */
  /*----------------------------*/
  assert(_cell_map.find(parent->get_id()) != _cell_map.end());
  Node* parent_node = _cell_map[parent->get_id()];
  
  /*----------------------------*/
  /* 2) Create new nodes        */
  /*----------------------------*/
  _current_id++;
  Node* node1 = new Node(_current_id, child1);
  _current_id++;
  Node* node2 = new Node(_current_id, child2);
  
  /*----------------------------*/
  /* 3) Update nodes attributes */
  /*----------------------------*/
  node1->set_parent(parent_node);
  node2->set_parent(parent_node);
  parent_node->add_child(node1);
  parent_node->add_child(node2);
  
  /*----------------------------*/
  /* 4) Add nodes to the tree   */
  /*----------------------------*/
  _node_map[node1->get_id()] = node1;
  _node_map[node2->get_id()] = node2;
  
  /*----------------------------*/
  /* 5) Update the cells map    */
  /*----------------------------*/
  _cell_map[child1->get_id()] = node1;
  _cell_map[child2->get_id()] = node2;
}

/**
 * \brief    Freeze the node
 * \details  Save the cell state by copying and killing it
 * \param    unsigned long long int cell_identifier
 * \param    node_state state
 * \param    size_t death_time
 * \return   \e void
 */
void Tree::freeze_node( unsigned long long int cell_identifier, size_t death_time )
{
  assert(_cell_map.find(cell_identifier) != _cell_map.end());
  _cell_map[cell_identifier]->set_dead(death_time);
}

/**
 * \brief    Delete a node and remove all links
 * \details  --
 * \param    unsigned long long int node_identifier
 * \return   \e void
 */
void Tree::delete_node( unsigned long long int node_identifier )
{
  assert(_node_map.find(node_identifier) != _node_map.end());
  Node* node = _node_map[node_identifier];
  assert(node->get_id() == node_identifier);
  assert(!node->isAlive());
  
  /*----------------------------------*/
  /* 1) Update parental children list */
  /*----------------------------------*/
  node->get_parent()->replace_children(node);
  
  /*-----------------------------------*/
  /* 2) Set the new parent of children */
  /*-----------------------------------*/
  for (size_t i = 0; i < node->get_parent()->get_number_of_children(); i++)
  {
    node->get_parent()->get_child(i)->set_parent(node->get_parent());
  }
  
  /*----------------------------------*/
  /* 3) Delete node                   */
  /*----------------------------------*/
  delete _node_map[node_identifier];
  _node_map[node_identifier] = NULL;
  _node_map.erase(node_identifier);
}

/**
 * \brief    Clean the cell map
 * \details  Remove dead nodes from the cell map
 * \param    void
 * \return   \e void
 */
void Tree::clean_cell_map( void )
{
  std::vector<unsigned long long int> to_remove;
  for (_iterator = _cell_map.begin(); _iterator != _cell_map.end(); ++_iterator)
  {
    if (!_iterator->second->isAlive())
    {
      to_remove.push_back(_iterator->first);
    }
    else
    {
      assert(_iterator->second->get_alive_cell()->isAlive());
    }
  }
  for (size_t i = 0; i < to_remove.size(); i++)
  {
    _cell_map[to_remove[i]] = NULL;
    _cell_map.erase(to_remove[i]);
  }
}

/**
 * \brief    Prune the tree
 * \details  Remove all dead branches
 * \param    void
 * \return   \e void
 */
void Tree::prune()
{
  untag_tree();
  
  /*-------------------------------------*/
  /* 1) Tag alive cells lineage          */
  /*-------------------------------------*/
  for (_iterator = _node_map.begin(); _iterator != _node_map.end(); ++_iterator)
  {
    assert(_iterator->first == _iterator->second->get_id());
    if (_iterator->second->isAlive())
    {
      _iterator->second->tag_lineage();
    }
  }
  
  /*-------------------------------------*/
  /* 2) Build the list of untagged nodes */
  /*-------------------------------------*/
  std::vector<unsigned long long int> remove_list;
  for (_iterator = _node_map.begin(); _iterator != _node_map.end(); ++_iterator)
  {
    if (!_iterator->second->isTagged() && !_iterator->second->isMasterRoot())
    {
      remove_list.push_back(_iterator->first);
    }
  }
  
  /*-------------------------------------*/
  /* 3) Delete untagged nodes            */
  /*-------------------------------------*/
  for (size_t i = 0; i < remove_list.size(); i++)
  {
    delete_node(remove_list[i]);
  }
  remove_list.clear();
  
  /*-------------------------------------*/
  /* 4) Set master root children as root */
  /*-------------------------------------*/
  Node* master_root = _node_map[0];
  for (size_t i = 0; i < master_root->get_number_of_children(); i++)
  {
    master_root->get_child(i)->set_root();
  }
}

/**
 * \brief    Shorten the tree
 * \details  Remove all the dead nodes that are not common ancestors
 * \param    void
 * \return   \e void
 */
void Tree::shorten()
{
  /*-------------------------------------*/
  /* 1) Select all intermediate nodes:   */
  /*    - not master root                */
  /*    - not alive                      */
  /*    - possessing exactly one child   */
  /*-------------------------------------*/
  std::vector<unsigned long long int> remove_list;
  remove_list.clear();
  for (_iterator = _node_map.begin(); _iterator != _node_map.end(); ++_iterator)
  {
    assert(_iterator->first == _iterator->second->get_id());
    if (!_iterator->second->isMasterRoot() && !_iterator->second->isAlive() && _iterator->second->get_number_of_children() == 1)
    {
      remove_list.push_back(_iterator->first);
    }
  }
  /*-------------------------------------*/
  /* 2) Delete nodes                     */
  /*-------------------------------------*/
  for (size_t i = 0; i < remove_list.size(); i++)
  {
    delete_node(remove_list[i]);
  }
  remove_list.clear();

#if DEBUG
  for (_iterator = _node_map.begin(); _iterator != _node_map.end(); ++_iterator)
  {
    if (!_iterator->second->isMasterRoot() && !_iterator->second->isAlive())
    {
      assert(_iterator->second->get_number_of_children() == 2);
    }
  }
#endif
  
  /*-------------------------------------*/
  /* 3) Set master root children as root */
  /*-------------------------------------*/
  Node* master_root = _node_map[0];
  for (size_t i = 0; i < master_root->get_number_of_children(); i++)
  {
    master_root->get_child(i)->set_root();
  }
}

/**
 * \brief    Write tree
 * \details  Write adjacency list in a file (.txt)
 * \param    std::string filename
 * \return   \e void
 */
void Tree::write_tree( std::string filename )
{
  std::ofstream file(filename.c_str(), std::ios::out | std::ios::trunc);
  for (_iterator = _node_map.begin(); _iterator != _node_map.end(); ++_iterator)
  {
    for (size_t i = 0; i < _iterator->second->get_number_of_children(); i++)
    {
      file << _iterator->second->get_id() << " " << _iterator->second->get_child(i)->get_id() << "\n";
    }
  }
  file.close();
}

/**
 * \brief    Write Newick tree
 * \details  Write tree in Newick format in a file (.phb)
 * \param    std::string filename
 * \return   \e void
 */
void Tree::write_newick_tree( std::string filename )
{
  std::ofstream file(filename.c_str(), std::ios::out | std::ios::trunc);
  
  for (size_t i = 0; i < _node_map[0]->get_number_of_children(); i++)
  {
    std::stringstream newick_tree;
    inOrderNewick(_node_map[0]->get_child(i), 0, newick_tree);
    newick_tree << ";\n";
    file << newick_tree.str();
  }
  file.close();
}

/**
 * \brief    Write lineage statistics from one specified node
 * \details  --
 * \param    std::string filename
 * \param    unsigned long long int identifier
 * \return   \e void
 */
void Tree::write_lineage_statistics( std::string filename, unsigned long long int identifier )
{
  /*---------------------------*/
  /* 1) Initialize data        */
  /*---------------------------*/
  std::stringstream genome_structure_filename;
  std::stringstream inherited_proteins_filename;
  std::stringstream phenotype_filename;
  std::stringstream fixed_mutations_filename;
  
  genome_structure_filename << "./statistics/genome_structure_" << filename << ".txt";
  inherited_proteins_filename << "./statistics/inherited_proteins_" << filename << ".txt";
  phenotype_filename << "./statistics/phenotype_" << filename << ".txt";
  fixed_mutations_filename << "./statistics/fixed_mutations_" << filename << ".txt";
  
  std::ofstream genome_structure_file(genome_structure_filename.str().c_str(), std::ios::out | std::ios::trunc);
  std::ofstream inherited_proteins_file(inherited_proteins_filename.str().c_str(), std::ios::out | std::ios::trunc);
  std::ofstream phenotype_file(phenotype_filename.str().c_str(), std::ios::out | std::ios::trunc);
  std::ofstream fixed_mutations_file(fixed_mutations_filename.str().c_str(), std::ios::out | std::ios::trunc);
  
  if (_node_map.find(identifier) == _node_map.end())
  {
    printf("Error in Tree::write_lineage_statistics(): call to a node that do not exist. Exit.\n");
    exit(EXIT_FAILURE);
  }
  
  /*---------------------------*/
  /* 2) Write headers          */
  /*---------------------------*/
  Node* node = _node_map[identifier];
  node->get_replication_report()->write_genome_structure_header(genome_structure_file);
  node->get_replication_report()->write_inherited_proteins_header(inherited_proteins_file);
  node->get_replication_report()->write_phenotype_header(phenotype_file);
  node->get_replication_report()->write_fixed_mutations_header(fixed_mutations_file);
  
  /*---------------------------*/
  /* 3) Climb the lineage tree */
  /*---------------------------*/
  while (!node->isMasterRoot())
  {
    node->get_replication_report()->write_genome_structure_data(genome_structure_file);
    node->get_replication_report()->write_inherited_proteins_data(inherited_proteins_file);
    node->get_replication_report()->write_phenotype_data(phenotype_file);
    node->get_replication_report()->write_fixed_mutations_data(fixed_mutations_file);
    node = node->get_parent();
  }
  
  /*---------------------------*/
  /* 4) Close files            */
  /*---------------------------*/
  genome_structure_file.close();
  inherited_proteins_file.close();
  phenotype_file.close();
  fixed_mutations_file.close();
}

/**
 * \brief    Write phylogeny statistics
 * \details  --
 * \param    std::string filename
 * \return   \e void
 */
void Tree::write_phylogeny_statistics( std::string filename )
{
  /*----------------------------------*/
  /* 1) Initialize data               */
  /*----------------------------------*/
  std::ofstream file(filename.c_str(), std::ios::out | std::ios::trunc);
  
  /*----------------------------------*/
  /* 2) Write headers                 */
  /*----------------------------------*/
  Node* node = get_first_node();
  node->get_replication_report()->write_replication_report_header(file);
  
  /*----------------------------------*/
  /* 3) Explore the phylogenetic tree */
  /*----------------------------------*/
  while (node != NULL)
  {
    if (!node->isMasterRoot())
    {
      node->get_replication_report()->write_replication_report_data(file);
    }
    node = get_next_node();
  }
  
  /*----------------------------------*/
  /* 4) Close files                   */
  /*----------------------------------*/
  file.close();
}

/**
 * \brief    Write trophic data
 * \details  --
 * \param    std::string filename
 * \return   \e void
 */
void Tree::write_trophic_data( std::string filename )
{
  /*----------------------------------*/
  /* 1) Initialize data               */
  /*----------------------------------*/
  std::ofstream file(filename.c_str(), std::ios::out | std::ios::trunc);
  
  /*----------------------------------*/
  /* 2) Write headers                 */
  /*----------------------------------*/
  Node* node = get_first_node();
  file << "id trophic_level trophic_group\n";
  
  /*----------------------------------*/
  /* 3) Explore the phylogenetic tree */
  /*----------------------------------*/
  while (node != NULL)
  {
    if (node->isAlive())
    {
      file << node->get_id() << " " << node->get_alive_cell()->get_trophic_level() << " " << node->get_alive_cell()->get_trophic_group() << "\n";
    }
    node = get_next_node();
  }
  
  /*----------------------------------*/
  /* 4) Close files                   */
  /*----------------------------------*/
  file.close();
}

/**
 * \brief    Compute the AI score
 * \details  Compute the AI score on the observed tree and then compute the distribution of AI scores on bootstraped trees
 * \param    size_t backup_time
 * \return   \e void
 */
void Tree::compute_AI_score_on_SL( size_t backup_time )
{
  std::stringstream filename;
  filename << "SL_stat_" << backup_time << ".txt";
  std::ofstream file(filename.str(), std::ios::out | std::ios::trunc);
  
  /*-------------------------------------------*/
  /* 1) Count the number of S and L cells      */
  /*-------------------------------------------*/
  double total_L = 0;
  double total_S = 0;
  for (_iterator = _node_map.begin(); _iterator != _node_map.end(); ++_iterator)
  {
    if (_iterator->second->isAlive())
    {
      trophic_level level = _iterator->second->get_alive_cell()->get_trophic_level();
      if (level == LEVEL_0 || level == LEVEL_1)
      {
        total_L += 1.0;
      }
      else if (level == LEVEL_2 || level == NO_LEVEL)
      {
        total_S += 1.0;
      }
    }
  }
  
  /*-------------------------------------------*/
  /* 2) Compute AI score                       */
  /*-------------------------------------------*/
  double AI = 0.0;
  std::unordered_map<unsigned long long int, Node*>::iterator it;
  for (it = _node_map.begin(); it != _node_map.end(); ++it)
  {
    if (!it->second->isAlive())
    {
      std::vector<Node*> tagged_nodes;
      tag_offspring(it->second, &tagged_nodes);
      double total  = 0.0;
      double Lcount = 0.0;
      for (size_t i = 0; i < tagged_nodes.size(); i++)
      {
        if (tagged_nodes[i]->isAlive())
        {
          total += 1.0;
          trophic_level level = tagged_nodes[i]->get_alive_cell()->get_trophic_level();
          if (level == LEVEL_0 || level == LEVEL_1)
          {
            Lcount += 1.0;
          }
        }
      }
      tagged_nodes.clear();
      if (total > 0.0)
      {
        double fi = (Lcount/total > 0.5 ? Lcount/total : 1.0-Lcount/total);
        AI += (1.0-fi)/pow(2.0, total-1.0);
      }
    }
  }
  double ca = get_common_ancestor_age();
  file << "backup depth AI L S\n";
  file << backup_time << " " << (double)backup_time-ca << " " << AI << " " << total_L << " " << total_S << "\n";
  file.close();
  
  /*-------------------------------------------*/
  /* 3) Compute AI score on bootstrapped trees */
  /*-------------------------------------------*/
  /* The structure of the tree is maintained,   */
  /* but the type is drawn at random such that: */
  /* pS = S/(S+L) and pL = L/(S+L)              */
  filename.str("");
  filename << "SL_neutral_" << backup_time << ".txt";
  file.open(filename.str(), std::ios::out | std::ios::trunc);
  file << "AI\n";
  size_t ITERATIONS = 10000;
  Prng* prng = new Prng();
  prng->set_seed((unsigned long int)time(NULL));
  double* scores = new double[ITERATIONS];
  for (size_t rep = 0; rep < ITERATIONS; rep++)
  {
    scores[rep] = 0.0;
  }
  for (it = _node_map.begin(); it != _node_map.end(); ++it)
  {
    if (!it->second->isAlive())
    {
      std::vector<Node*> tagged_nodes;
      tag_offspring(it->second, &tagged_nodes);
      double total = 0.0;
      for (size_t i = 0; i < tagged_nodes.size(); i++)
      {
        if (tagged_nodes[i]->isAlive())
        {
          total += 1.0;
        }
      }
      tagged_nodes.clear();
      if (total > 0.0)
      {
        for (size_t rep = 0; rep < ITERATIONS; rep++)
        {
          double Lcount = (double)prng->binomial(total, total_L/(total_L+total_S));
          double fi = (Lcount/total > 0.5 ? Lcount/total : 1.0-Lcount/total);
          scores[rep] += (1.0-fi)/pow(2.0, total-1.0);
        }
      }
    }
  }
  
  /*-------------------------------------------*/
  /* 4) Write randomized scores in file        */
  /*-------------------------------------------*/
  for (size_t rep = 0; rep < ITERATIONS; rep++)
  {
    file << scores[rep] << "\n";
  }
  file.close();
  
  delete[] scores;
  scores = NULL;
}

/**
 * \brief    Compute the common ancestor SL repartition
 * \details  Compute the repartition of S/L ecotypes in the phylogenetic tree, starting from CA offspring. Save the proportion of S in each offspring in Sp_A and Sp_B
 * \param    double& Sp_1
 * \param    double& Sp_2
 * \return   \e void
 */
void Tree::compute_common_ancestor_SL_repartition( double& Sp_1, double& Sp_2 )
{
  /*---------------------*/
  /* 1) Get CA offspring */
  /*---------------------*/
  Node* CA = get_common_ancestor();
  if (CA == NULL)
  {
    Sp_1 = -1.0;
    Sp_2 = -1.0;
    return;
  }
  if (CA->get_number_of_children() > 2)
  {
    Sp_1 = -1.0;
    Sp_2 = -1.0;
    return;
  }
  Node* offspring1 = CA->get_child(0);
  Node* offspring2 = CA->get_child(1);
  
  /*---------------------*/
  /* 2) Compute Sp_1     */
  /*---------------------*/
  double nb_A = 0.0;
  double nb_B = 0.0;
  Node* node = get_first_node();
  while (node != NULL)
  {
    if (node->isAlive() && node->isAncestor(offspring1->get_id()))
    {
      trophic_level level = node->get_alive_cell()->get_trophic_level();
      if (level == LEVEL_0 || level == LEVEL_1)
      {
        nb_A += 1.0;
      }
      else if (level == LEVEL_2 || level == NO_LEVEL)
      {
        nb_B += 1.0;
      }
    }
    node = get_next_node();
  }
  if (nb_A == 0.0 && nb_B == 0.0)
  {
    Sp_1 = -1.0;
    Sp_2 = -1.0;
    return;
  }
  Sp_1 = nb_B/(nb_A+nb_B);
  
  /*---------------------*/
  /* 3) Compute Sp_2     */
  /*---------------------*/
  nb_A = 0.0;
  nb_B = 0.0;
  node = get_first_node();
  while (node != NULL)
  {
    if (node->isAlive() && node->isAncestor(offspring2->get_id()))
    {
      trophic_level level = node->get_alive_cell()->get_trophic_level();
      if (level == LEVEL_0 || level == LEVEL_1)
      {
        nb_A += 1.0;
      }
      else if (level == LEVEL_2 || level == NO_LEVEL)
      {
        nb_B += 1.0;
      }
    }
    node = get_next_node();
  }
  if (nb_A == 0.0 && nb_B == 0.0)
  {
    Sp_1 = -1.0;
    Sp_2 = -1.0;
    return;
  }
  Sp_2 = nb_B/(nb_A+nb_B);
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Recursive method used to build Newick format file
 * \details  --
 * \param    Node* node
 * \param    size_t parent_time
 * \param    std::stringstream& output
 * \return   \e void
 */
void Tree::inOrderNewick( Node* node, size_t parent_time, std::stringstream& output )
{
  /*--------------------------------------*/
  /* 1) If node is a leaf                 */
  /*--------------------------------------*/
  if (node->isAlive() || node->get_number_of_children() < 2)
  {
    output << node->get_id() << ":" << node->get_replication_report()->get_birth_time()-parent_time;
  }
  
  /*--------------------------------------*/
  /* 2) Else if node has several children */
  /*--------------------------------------*/
  else
  {
    output << "(";
    for (size_t i = 0; i < node->get_number_of_children(); i++)
    {
      inOrderNewick(node->get_child(i), node->get_replication_report()->get_birth_time(), output);
      if (i < node->get_number_of_children()-1)
      {
        output << ", ";
      }
    }
    output << ")";
    output << node->get_replication_report()->get_id() << ":" << node->get_replication_report()->get_birth_time()-parent_time;
  }
}

/**
 * \brief    Tag all the nodes
 * \details  --
 * \param    void
 * \return   \e void
 */
void Tree::tag_tree()
{
  for (_iterator = _node_map.begin(); _iterator != _node_map.end(); ++_iterator)
  {
    _iterator->second->tag();
  }
}

/**
 * \brief    Untag all the nodes
 * \details  --
 * \param    void
 * \return   \e void
 */
void Tree::untag_tree()
{
  for (_iterator = _node_map.begin(); _iterator != _node_map.end(); ++_iterator)
  {
    _iterator->second->untag();
  }
}

/**
 * \brief    Tag all the offspring of this node
 * \details  --
 * \param    Node* node
 * \param    std::vector<Node*>* tagged_nodes
 * \return   \e void
 */
void Tree::tag_offspring( Node* node, std::vector<Node*>* tagged_nodes )
{
  untag_tree();
  tagged_nodes->clear();
  tagged_nodes->push_back(node);
  bool end = false;
  while (!end)
  {
    end = true;
    for (size_t i = 0; i < tagged_nodes->size(); i++)
    {
      for (size_t j = 0; j < tagged_nodes->at(i)->get_number_of_children(); j++)
      {
        if (!tagged_nodes->at(i)->get_child(j)->isTagged())
        {
          end = false;
          tagged_nodes->at(i)->get_child(j)->tag();
          tagged_nodes->push_back(tagged_nodes->at(i)->get_child(j));
        }
      }
    }
  }
}

