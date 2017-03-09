
/**
 * \file      Structs.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Definition of structures
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

#ifndef __EVOEVO__Structs__
#define __EVOEVO__Structs__

#include <iostream>
#include <vector>
#include <cstring>
#include <unordered_map>
#include <gsl/gsl_odeiv2.h>

#include "Macros.h"
#include "Enums.h"

class Genome;
class InheritedProteins;
class SpeciesList;


/******************************************************************************************/

/**
 * \brief   Genetic unit struct
 * \details Defines the structure of a genetic unit
 */
typedef struct
{
  
  /*------------------------------------------------------------------ Global attributes */
  
  genetic_unit_type      type;              /*!< Type of genetic unit             */
  unsigned long long int identifier;        /*!< Genetic unit identifier          */
  unsigned long long int parent_identifier; /*!< Parental genetic unit identifier */
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  int    s;             /*!< Substrate species tag  */
  int    p;             /*!< Product species tag    */
  double kcat;          /*!< Kcat constant          */
  double kcat_km_ratio; /*!< Kcat/Km ratio constant */
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  int    BS_tag;         /*!< Binding site tag      */
  int    coE_tag;        /*!< Co-enzyme species tag */
  bool   free_activity;  /*!< Free activity         */
  bool   bound_activity; /*!< Bound activity        */
  size_t binding_window; /*!< Window size           */
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  int TF_tag; /*!< Transcription factor tag */
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  double basal_expression_level; /*!< Basal expression level of the promoter */
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  bool functional; /*!< Indicates if the genetic unit is functional or not */
  
} genetic_unit;

/******************************************************************************************/

/**
 * \brief   Genetic sequence struct
 * \details Defines the structure of a genetic sequence
 */
typedef struct
{
  
  size_t        size;        /*!< Size of the array      */
  size_t        buffer_size; /*!< Buffer size            */
  genetic_unit* x;           /*!< Array of genetic units */

} genetic_sequence;

/******************************************************************************************/

/**
 * \brief   Reaction list struct
 * \details Contains the list of reactions encoded in the genome
 */
typedef struct
{
  /*------------------------------------------------------------------ list of genetic regulation network equations */
  
  size_t                      grn_N;                         /*!< Number of GRN equations                             */
  std::vector<size_t>         grn_promoter;                  /*!< List of promoter positions                          */
  std::vector<size_t>         grn_Nenhancer;                 /*!< Number of enhancer bindings per promoter            */
  std::vector<size_t>         grn_Noperator;                 /*!< Number of operator bindings per promoter            */
  std::vector<size_t>         grn_Ngenes;                    /*!< Number of genes regulated per promoter              */
  std::vector<size_t>         grn_enhancer_nb_BS;            /*!< Number of binding site per enhancer                 */
  std::vector<size_t>         grn_enhancer_TF_list;          /*!< List of enhancer binding TF indexes per promoter    */
  std::vector<double>         grn_enhancer_affinity_list;    /*!< List of enhancer binding TF affinities per promoter */
  std::vector<int>            grn_enhancer_coe_list;         /*!< List of TF co-enzymes per enhancer                  */
  std::vector<co_enzyme_type> grn_enhancer_coe_type;         /*!< List of co-enzyme types per enhancer                */
  std::vector<double>         grn_enhancer_coe_km;           /*!< List of co-enzyme Km constant per enhancer          */
  std::vector<size_t>         grn_operator_nb_BS;            /*!< Number of binding site per operator                 */
  std::vector<size_t>         grn_operator_TF_list;          /*!< List of operator binding TF indexes per promoter    */
  std::vector<double>         grn_operator_affinity_list;    /*!< List of operator binding TF affinities per promoter */
  std::vector<int>            grn_operator_coe_list;         /*!< List of TF co-enzymes per operator                  */
  std::vector<co_enzyme_type> grn_operator_coe_type;         /*!< List of co-enzyme types per operator                */
  std::vector<double>         grn_operator_coe_km;           /*!< List of co-enzyme Km constant per operator          */
  std::vector<size_t>         grn_regulated_genes;           /*!< List of regulated genes per promoter                */
  std::vector<double>         grn_beta;                      /*!< List of basal expression levels per promoter        */
  double                      grn_hill_n;                    /*!< Parameter n of the Hill equation                    */
  double                      grn_theta_n;                   /*!< Constant = pow(hill_theta, hill_n)                  */
  double                      grn_protein_degradation_rate;  /*!< Protein degradation rate                            */
  double                      grn_energy_transcription_cost; /*!< Energetical cost of transcription                   */
  double                      grn_energy_degradation_cost;   /*!< Energetical cost of degradation                     */
  double                      grn_energy_dissipation_rate;   /*!< Energy dissipation rate                             */
  
  /*------------------------------------------------------------------ list of metabolic equations */
  
  size_t                       metabolic_N;       /*!< Number of metabolic equations */
  std::vector<reaction_type>   metabolic_type;    /*!< Type of the equation          */
  std::vector<reaction_origin> metabolic_origin;  /*!< Origin of the reaction        */
  std::vector<int>             metabolic_s;       /*!< substrate tag                 */
  std::vector<int>             metabolic_p;       /*!< Product tag                   */
  std::vector<double>          metabolic_km;      /*!< Km constant                   */
  std::vector<double>          metabolic_kcat;    /*!< Kcat constant                 */
  std::vector<double>          metabolic_delta_g; /*!< Energy cost of the reaction   */
  std::vector<size_t>          metabolic_e;       /*!< Enzyme index                  */
  
  /*------------------------------------------------------------------ size of each state vector subspace */
  
  size_t Ngenome;    /*!< Length of the genome state vector       */
  size_t Ninherited; /*!< Length of the inherited proteins vector */
  size_t Ncell;      /*!< Length of the cell state vector         */
  size_t Nenv;       /*!< Length of the environement state vector */
  size_t N;          /*!< Length of the state vector X            */
  
  /*------------------------------------------------------------------ modeling schemes */
  
  bool energy_constraints;    /*!< Indicates if energy constraints are activated   */
  bool membrane_permeability; /*!< Indicates if membrane permeability is activated */
  bool metabolic_inheritance; /*!< Indicates if metabolic inheritance is activated */
  bool enzymatic_inheritance; /*!< Indicates if enzymatic inheritance is activated */
  bool co_enzyme_activity;    /*!< Indicates if co-enzyme activity is activated    */
  
  /*------------------------------------------------------------------ membrane permeability */
  
  double kmp; /*!< Membrane permeability */
  
  /*------------------------------------------------------------------ metabolic exchanges with the environment */
  
  double metabolic_uptake;  /*!< Metabolic uptake  */
  double metabolic_release; /*!< Metabolic release */
  
  /*------------------------------------------------------------------ timestep ratio */
  
  double timestep_ratio; /*!< Timestep ratio between GRN and metabolism timesteps */
  
  /*------------------------------------------------------------------ environment interaction scheme */
  
  bool interacts_with_environment; /*!< Indicates if cells interact with the environment */
  
  /*------------------------------------------------------------------ test variables */
  
#ifdef DEBUG
  Genome*            genome;             /*!< Related cell's genome             */
  InheritedProteins* inherited_proteins; /*!< Related cell's inherited proteins */
  SpeciesList*       cell_species_list;  /*!< Related cell's species list       */
  SpeciesList*       env_species_list;   /*!< Local environment's species list  */
#endif
  
} reaction_list;

/******************************************************************************************/

/**
 * \brief   Variable range
 * \details --
 */
typedef struct
{
  
  double min; /*!< Minimum value of the range */
  double max; /*!< Maximum value of the range */
  
} variable_range;

/******************************************************************************************/

/**
 * \brief   Environment properties
 * \details --
 */
typedef struct
{
  
  size_t                          number_of_init_cycles;   /*!< Number of initialization cycles     */
  variable_range                  species_tag_range;       /*!< Environment species tag range       */
  variable_range                  concentration_range;     /*!< Environment concentration range     */
  variable_range                  number_of_species_range; /*!< Environment number of species range */
  environment_interaction_scheme  interaction_scheme;      /*!< Environment interaction scheme      */
  environment_renewal_scheme      renewal_scheme;          /*!< Environment renewal scheme          */
  environment_variation_scheme    variation_scheme;        /*!< Environment variation scheme        */
  environment_localization_scheme localization_scheme;     /*!< Environment localization scheme     */
  environment_metabolic_scheme    metabolic_scheme;        /*!< Environment metabolic scheme        */
  double                          introduction_rate;       /*!< Environment introduction rate       */
  double                          diffusion_coefficient;   /*!< Environment diffusion coefficient   */
  double                          degradation_rate;        /*!< Environment degradation rate        */
  
} environment_properties;

/******************************************************************************************/


#endif /* defined(__EVOEVO__Structs__) */
