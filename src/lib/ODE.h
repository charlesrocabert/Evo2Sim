
/**
 * \file      ODE.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     ODE class declaration
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

#ifndef __EVOEVO__ODE__
#define __EVOEVO__ODE__

#include <iostream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <unordered_map>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <assert.h>

#include "Macros.h"
#include "Enums.h"
#include "Structs.h"
#include "Parameters.h"
#include "ReplicationReport.h"
#include "Genome.h"
#include "InheritedProteins.h"
#include "SpeciesList.h"
#include "Environment.h"


class ODE
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  ODE( void ) = delete;
  ODE( Parameters* parameters, ReplicationReport* replication_report, Genome* genome, InheritedProteins* inherited_proteins, SpeciesList* species_list, Environment* environment, size_t x, size_t y, double* energy, double* metabolic_uptake, double* metabolic_release, size_t* grn_nb_nodes, size_t* grn_nb_edges, size_t* metabolic_nb_nodes, size_t* metabolic_nb_edges, std::vector<int>* inflowing_pumps, std::vector<int>* outflowing_pumps );
  ODE( const ODE& ode ) = delete;
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~ODE( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  inline reaction_list* get_reaction_list( void );
  inline double*        get_X( void );
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  ODE& operator=(const ODE&) = delete;
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void load( bool from_backup, bool new_individual );
  void build_reaction_list( void );
  void build_genome_concentration_vector( bool new_individual );
  void initialize_state( void );
  bool solve( void );
  void update( void );
  
  static int ODE_system_without_energy( double t, const double y[], double dydt[], void* parameters );
  static int ODE_system_with_energy( double t, const double y[], double dydt[], void* parameters );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  void build_system( void );
  void build_X( void );
  
  bool evaluate_promoter_region( size_t promoter_position, std::vector<size_t>* expressed_genes );
  void add_functional_region( size_t promoter_position, std::vector<size_t>* expressed_genes, std::unordered_map<size_t, unsigned char>& nodes_list );
  void add_metabolic_reaction( reaction_type type, reaction_origin origin, int s, int p, double km, double kcat, double delta_g, size_t e );
  void add_inherited_metabolic_reactions( std::unordered_map<size_t, unsigned char>& nodes_list );
  void compute_replication_report_statistics( void );
  void create_reaction_list( void );
  void free_reaction_list( void );
  
#ifdef DEBUG
  void test_reaction_list_structure( void );
#endif
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  
  /*------------------------------------------------------------------ simulation parameters */
  
  Parameters* _parameters; /*!< Parameters */
  
  /*------------------------------------------------------------------ cell variables */
  
  ReplicationReport* _replication_report; /*!< Replication report                       */
  Genome*            _genome;             /*!< Genome                                   */
  InheritedProteins* _inherited_proteins; /*!< Inherited proteins                       */
  SpeciesList*       _species_list;       /*!< List of molecule species                 */
  Environment*       _environment;        /*!< Environment                              */
  size_t             _x;                  /*!< Cell's x coordinate                      */
  size_t             _y;                  /*!< Cell's y coordinate                      */
  double*            _energy;             /*!< Cell's energy amount                     */
  double*            _metabolic_uptake;   /*!< Cell's metabolic uptake                  */
  double*            _metabolic_release;  /*!< Cell's metabolic release                 */
  size_t*            _grn_nb_nodes;       /*!< Number of nodes in the GRN               */
  size_t*            _grn_nb_edges;       /*!< Number of edges in the GRN               */
  size_t*            _metabolic_nb_nodes; /*!< Number of nodes in the metabolic network */
  size_t*            _metabolic_nb_edges; /*!< Number of edges in the metabolic network */
  std::vector<int>*  _inflowing_pumps;    /*!< List of the inflowing pumps              */
  std::vector<int>*  _outflowing_pumps;   /*!< List of the outflowing pumps             */
  
  /*------------------------------------------------------------------ solver variables */
  
  reaction_list*      _reaction_list; /*!< List of reaction rules */
  gsl_odeiv2_system   _system;        /*!< GSL solver system      */
  gsl_odeiv2_control* _control;       /*!< GSL solver control     */
  gsl_odeiv2_step*    _step;          /*!< GSL solver step        */
  gsl_odeiv2_evolve*  _evolve;        /*!< GSL solver evolve      */
  double*             _X;             /*!< State vector X         */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/**
 * \brief    Get the reaction list
 * \details  --
 * \param    void
 * \return   \e reaction_list*
 */
inline reaction_list* ODE::get_reaction_list( void )
{
  return _reaction_list;
}

/**
 * \brief    Get the state vector X
 * \details  --
 * \param    void
 * \return   \e double*
 */
inline double* ODE::get_X( void )
{
  return _X;
}

/*----------------------------
 * SETTERS
 *----------------------------*/


#endif /* defined(__EVOEVO__ODE__) */
