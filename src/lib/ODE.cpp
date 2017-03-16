
/**
 * \file      ODE.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     ODE class definition
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

#include "ODE.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters
 * \param    ReplicationReport* replication_report
 * \param    Genome* genome
 * \param    InheritedProteins* inherited_proteins
 * \param    SpeciesList* species_list
 * \param    Environment* environment
 * \param    size_t x
 * \param    size_t y
 * \param    double* energy
 * \param    double* metabolic_uptake
 * \param    double* metabolic_release
 * \param    size_t* grn_nb_nodes
 * \param    size_t* grn_nb_edges
 * \param    size_t* metabolic_nb_nodes
 * \param    size_t* metabolic_nb_edges
 * \param    std::vector<int>* inflowing_pumps
 * \param    std::vector<int>* outflowing_pumps
 * \return   \e void
 */
ODE::ODE( Parameters* parameters, ReplicationReport* replication_report, Genome* genome, InheritedProteins* inherited_proteins, SpeciesList* species_list, Environment* environment, size_t x, size_t y, double* energy, double* metabolic_uptake, double* metabolic_release, size_t* grn_nb_nodes, size_t* grn_nb_edges, size_t* metabolic_nb_nodes, size_t* metabolic_nb_edges, std::vector<int>* inflowing_pumps, std::vector<int>* outflowing_pumps )
{
  
  /*------------------------------------------------------------------ simulation parameters */
  
  _parameters = parameters;
  
  /*------------------------------------------------------------------ cell variables */
  
  _replication_report = replication_report;
  _genome             = genome;
  _inherited_proteins = inherited_proteins;
  _species_list       = species_list;
  _environment        = environment;
  _x                  = x;
  _y                  = y;
  _energy             = energy;
  _metabolic_uptake   = metabolic_uptake;
  _metabolic_release  = metabolic_release;
  _grn_nb_nodes       = grn_nb_nodes;
  _grn_nb_edges       = grn_nb_edges;
  _metabolic_nb_nodes = metabolic_nb_nodes;
  _metabolic_nb_edges = metabolic_nb_edges;
  _inflowing_pumps    = inflowing_pumps;
  _outflowing_pumps   = outflowing_pumps;
  
  /*------------------------------------------------------------------ solver variables */
  
  _reaction_list = NULL;
  _control       = NULL;
  _step          = NULL;
  _evolve        = NULL;
  _X             = NULL;
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
ODE::~ODE( void )
{
  free_reaction_list();
  gsl_odeiv2_control_free(_control);
  _control = NULL;
  delete[] _X;
  _X = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Load ODE system
 * \details  If the ODE system is loaded just after backup loading, do not modify the genome state vector
 * \param    bool from_backup
 * \param    bool new_individual
 * \return   \e void
 */
void ODE::load( bool from_backup, bool new_individual )
{
  /*--------------------------------------------------------------*/
  /* 1) Build the reaction list, one time for the whole cell life */
  /*--------------------------------------------------------------*/
  build_reaction_list();
  compute_replication_report_statistics();
  
  /*--------------------------------------------------------------*/
  /* 2) Build the genome concentration vector                     */
  /*--------------------------------------------------------------*/
  if (!from_backup)
  {
    build_genome_concentration_vector(new_individual);
  }
  
  /*--------------------------------------------------------------*/
  /* 3) Build the ODE system, one time for the whole cell life    */
  /*--------------------------------------------------------------*/
  build_system();
}

/**
 * \brief    Build the reaction list
 * \details  Parse the genome, and save in the reaction list all available reactions
 * \param    void
 * \return   \e void
 */
void ODE::build_reaction_list( void )
{
  /*-------------------------------------------------*/
  /* 1) Create the reaction list                     */
  /*-------------------------------------------------*/
  create_reaction_list();
  *(_metabolic_nb_nodes) = 0;
  *(_metabolic_nb_edges) = 0;
  std::unordered_map<size_t, unsigned char> nodes_list;
  nodes_list.clear();
  _inflowing_pumps->clear();
  _outflowing_pumps->clear();
  
  /*-------------------------------------------------*/
  /* 2) Load reaction rules list                     */
  /*-------------------------------------------------*/
  /*
   For each promoter, if the promoter region is functional
   (i.e. mandatory presence of E and/or TF types, and optional
   presence of an enhancer - before the P - and an operator
   - after the P -, made of BS types):
   -> Parse the functional region and add the GRN equation
    to the reaction list,
   -> Add metabolic reactions if the functional region
      contains regulated E types.
   */
  size_t* Pi = _genome->get_Pi();
  for (size_t promoter_index = 0; promoter_index < _genome->get_nb_P(); promoter_index++)
  {
    std::vector<size_t> expressed_genes;
    expressed_genes.clear();
    if(evaluate_promoter_region(Pi[promoter_index], &expressed_genes))
    {
      add_functional_region(Pi[promoter_index], &expressed_genes, nodes_list);
    }
  }
  
  /*-------------------------------------------------*/
  /* 3) Load inherited proteins in the reaction list */
  /*-------------------------------------------------*/
  if (_parameters->get_enzymatic_inheritance_scheme())
  {
    add_inherited_metabolic_reactions(nodes_list);
  }
  nodes_list.clear();
}

/**
 * \brief    Build the genome concentration vector
 * \details  The genome concentration vector needs to be initialized. Depending if enzymatic inheritance is activated or not, initial concentrations are zero (inherited concentrations being saved in the InheritedProteins class) or are the basal expression level of each related promoter. In case of new individual creation, enzymatic concentrations are beta.
 * \param    bool new_individual
 * \return   \e void
 */
void ODE::build_genome_concentration_vector( bool new_individual )
{
  double* concentration_vector = _genome->get_concentration_vector();
  for (size_t i = 0; i < _genome->get_size(); i++)
  {
    concentration_vector[i] = 0.0;
  }
  if (!_parameters->get_enzymatic_inheritance_scheme() || new_individual)
  {
    size_t regulated_genes_index = 0;
    for (size_t i = 0; i < _reaction_list->grn_N; i++)
    {
      for (size_t j = 0; j < _reaction_list->grn_Ngenes[i]; j++)
      {
        concentration_vector[_reaction_list->grn_regulated_genes[regulated_genes_index]] = _reaction_list->grn_beta[i]/_reaction_list->grn_protein_degradation_rate;
        regulated_genes_index++;
      }
    }
  }
}

/**
 * \brief    Initialize X and reaction_list states
 * \details  --
 * \param    void
 * \return   \e void
 */
void ODE::initialize_state( void )
{
  /*-------------------------------------------------------------*/
  /* 1) Set the length of the state vector X and each subspaces, */
  /*    set other parameters                                     */
  /*-------------------------------------------------------------*/
  _reaction_list->Ngenome    = _genome->get_size();
  _reaction_list->Ninherited = 0;
  if (_parameters->get_enzymatic_inheritance_scheme())
  {
    _reaction_list->Ninherited = _inherited_proteins->get_size();
  }
  _reaction_list->Ncell = _species_list->get_size();
  _reaction_list->Nenv  = _environment->get_size(_x, _y);
  _reaction_list->N     = _reaction_list->Ngenome+_reaction_list->Ninherited+_reaction_list->Ncell+_reaction_list->Nenv+3;
  
  /*-------------------------------------------------------------*/
  /* 2) Set the dimension of the ODE system parameters space     */
  /*-------------------------------------------------------------*/
  _system.dimension = static_cast<size_t>(_reaction_list->N);
  
  /*-------------------------------------------------------------*/
  /* 3) Build X state vector                                     */
  /*-------------------------------------------------------------*/
  build_X();
}

/**
 * \brief    Solve the ODE system for one timestep
 * \details  --
 * \param    void
 * \return   \e bool
 */
bool ODE::solve( void )
{
  /*-------------------------------------------------------------*/
  /* 1) Set the length of the state vector X and each subspaces, */
  /*    set other parameters                                     */
  /*-------------------------------------------------------------*/
  _reaction_list->Ngenome    = _genome->get_size();
  _reaction_list->Ninherited = 0;
  if (_parameters->get_enzymatic_inheritance_scheme())
  {
    _reaction_list->Ninherited = _inherited_proteins->get_size();
  }
  _reaction_list->Ncell = _species_list->get_size();
  _reaction_list->Nenv  = _environment->get_size(_x, _y);
  _reaction_list->N     = _reaction_list->Ngenome+_reaction_list->Ninherited+_reaction_list->Ncell+_reaction_list->Nenv+3;
  
  /*-------------------------------------------------------------*/
  /* 2) Set the dimension of the ODE system parameters space     */
  /*-------------------------------------------------------------*/
  _system.dimension = static_cast<size_t>(_reaction_list->N);
  
  /*-------------------------------------------------------------*/
  /* 3) Build X state vector                                     */
  /*-------------------------------------------------------------*/
  build_X();
  
  /*-------------------------------------------------------------*/
  /* 4) Solve the ODE system                                     */
  /*-------------------------------------------------------------*/
  size_t N       = _reaction_list->N;
  _step          = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkck, N);
  _evolve        = gsl_odeiv2_evolve_alloc(N);
  double t_start = 0.0;
  double t_end   = _parameters->get_metabolism_timestep();
  double h       = H_INIT;
  
  while (t_start < t_end)
  {
    /*********************************************************/
    /* 4.1) Solve the system for t_start-t_end ODE timesteps */
    /*********************************************************/
    double t_start_old = t_start;
    double h_old       = h;
    
    int status  = gsl_odeiv2_evolve_apply(_evolve, _control, _step, &_system, &t_start, t_end, &h, _X);
    
    if (status != GSL_SUCCESS)
    {
      fprintf(stderr, "Error: an error occured during ODE solving. Exit.\n");
      exit(EXIT_FAILURE);
    }
    
    /*********************************************************/
    /* 4.2) Evaluate the state vector X                      */
    /*********************************************************/
    size_t INDEX_MAX = N-3;
    if (_parameters->get_energy_costs_scheme())
    {
      //INDEX_MAX = N-2;
      INDEX_MAX = N-3;
    }
    bool reset_ODE_system = false;
    bool decrease_h       = false;
    
    /* For each value of X, excepted last values (pos=INDEX_MAX), if the
       value is below -MINIMUM_CONCENTRATION, a significative error is detected,
       the system is reset to previous t_start, and h is reduced.
       If the value is between -MINIMUM_CONCENTRATION and MINIMUM_CONCENTRATION,
       it is simply set to zero.
     */
    for (size_t i = 0; i < INDEX_MAX; i++)
    {
      if (_X[i] <= -MINIMUM_CONCENTRATION)
      {
        decrease_h = true;
      }
      
      else if (_X[i] > -MINIMUM_CONCENTRATION && _X[i] <= MINIMUM_CONCENTRATION)
      {
        _X[i]            = 0.0;
        reset_ODE_system = true;
      }
    }
    
    /*********************************************************/
    /* 4.3) Apply potential modifications and corrections    */
    /*********************************************************/
    if (reset_ODE_system)
    {
      gsl_odeiv2_evolve_reset(_evolve);
    }
    if (decrease_h)
    {
      //std::cout << "DECREASE H TO " << h << " " << t_start << "\n";
      gsl_odeiv2_evolve_reset(_evolve);
      gsl_odeiv2_step_reset(_step);
      memcpy(_X, _evolve->y0, sizeof(double)*N);
      t_start = t_start_old;
      h       = h_old/2.0; /* GSL increases h by 5 at most. */
    }
  }
  
  gsl_odeiv2_step_free(_step);
  _step = NULL;
  gsl_odeiv2_evolve_free(_evolve);
  _evolve = NULL;
  
#ifdef DEBUG
  test_reaction_list_structure();
#endif
  
  return true;
}

/**
 * \brief    Update cell and environment species lists, and energy
 * \details  Return false if a negative value appears, true else
 * \param    void
 * \return   \e void
 */
void ODE::update( void )
{
  size_t Ngenome    = _reaction_list->Ngenome;
  size_t Ninherited = _reaction_list->Ninherited;
  size_t Ncell      = _reaction_list->Ncell;
  
  /*------------------------------------------------------------------*/
  /* 1) Load the new set of enzymatic concentrations                  */
  /*------------------------------------------------------------------*/
  double* genome_vector = _genome->get_concentration_vector();
  for (size_t index = 0; index < Ngenome; index++)
  {
    genome_vector[index] = _X[index];
  }
  
  /*------------------------------------------------------------------*/
  /* 2) Load the new set of inherited enzymatic concentrations        */
  /*------------------------------------------------------------------*/
  if (_parameters->get_enzymatic_inheritance_scheme())
  {
    double* inherited_vector = _inherited_proteins->get_concentration_vector();
    for (size_t index = 0; index < Ninherited; index++)
    {
      inherited_vector[index] = _X[Ngenome+index];
    }
  }
  
  /*------------------------------------------------------------------*/
  /* 3) Evaluate the amount of eaten and released metabolites and     */
  /*    load new concentrations in cell and environment species lists */
  /*------------------------------------------------------------------*/
  for (size_t index = 0; index < Ncell; index++)
  {
    _species_list->set((int)index+1, _X[Ngenome+Ninherited+index]);
    if (_parameters->get_environment_properties()->interaction_scheme == INTERACTION)
    {
      _environment->set(_x, _y, (int)index+1, _X[Ngenome+Ninherited+Ncell+index]);
    }
  }
  
  /*------------------------------------------------------------------*/
  /* 4) Load other cell's variables                                   */
  /*------------------------------------------------------------------*/
  *(_energy)            = _X[_reaction_list->N-3];
  *(_metabolic_uptake)  = _X[_reaction_list->N-2];
  *(_metabolic_release) = _X[_reaction_list->N-1];
  
  /*------------------------------------------------------------------*/
  /* 5) Free the state vector                                         */
  /*------------------------------------------------------------------*/
  delete[] _X;
  _X = NULL;
}

/**
 * \brief    Compute the ODE system without energy
 * \details  --
 * \param    double t
 * \param    const double y[]
 * \param    double dydt[]
 * \param    void* parameters
 * \return   \e int GSL_SUCCESS
 */
int ODE::ODE_system_without_energy( double t, const double y[], double dydt[], void* parameters )
{
  (void)t;
  
  /******************************/
  /* 1) Initialize variables    */
  /******************************/
  
  reaction_list*               list                     = (reaction_list*)parameters;
  size_t                       N                        = list->N;
  size_t                       Ngenome                  = list->Ngenome;
  size_t                       Ninherited               = list->Ninherited;
  size_t                       Ncell                    = list->Ncell;
  double                       kmp                      = list->kmp;
  std::vector<size_t>&         Nenhancer                = list->grn_Nenhancer;
  std::vector<size_t>&         Noperator                = list->grn_Noperator;
  std::vector<size_t>&         Ngenes                   = list->grn_Ngenes;
  std::vector<size_t>&         enhancer_TF_list         = list->grn_enhancer_TF_list;
  std::vector<double>&         enhancer_affinity_list   = list->grn_enhancer_affinity_list;
  std::vector<int>&            enhancer_coe_list        = list->grn_enhancer_coe_list;
  std::vector<co_enzyme_type>& enhancer_coe_type        = list->grn_enhancer_coe_type;
  std::vector<double>&         enhancer_coe_km          = list->grn_enhancer_coe_km;
  std::vector<size_t>&         operator_TF_list         = list->grn_operator_TF_list;
  std::vector<double>&         operator_affinity_list   = list->grn_operator_affinity_list;
  std::vector<int>&            operator_coe_list        = list->grn_operator_coe_list;
  std::vector<co_enzyme_type>& operator_coe_type        = list->grn_operator_coe_type;
  std::vector<double>&         operator_coe_km          = list->grn_operator_coe_km;
  std::vector<size_t>&         regulated_genes          = list->grn_regulated_genes;
  std::vector<double>&         beta                     = list->grn_beta;
  double                       hill_n                   = list->grn_hill_n;
  double                       theta_n                  = list->grn_theta_n;
  double                       protein_degradation_rate = list->grn_protein_degradation_rate;
  double                       timestep_ratio           = list->timestep_ratio;
  bool                         coe_activity             = list->co_enzyme_activity;
  bool                         membrane_permeability    = list->membrane_permeability;
  bool                         enzymatic_inheritance    = list->enzymatic_inheritance;
  std::vector<reaction_type>&  type                     = list->metabolic_type;
  std::vector<int>&            s                        = list->metabolic_s;
  std::vector<int>&            p                        = list->metabolic_p;
  std::vector<double>&         km                       = list->metabolic_km;
  std::vector<double>&         kcat                     = list->metabolic_kcat;
  std::vector<size_t>&         e                        = list->metabolic_e;
  bool                         env_interaction          = list->interacts_with_environment;
  
  /******************************/
  /* 2) Initialize dy/dt vector */
  /******************************/
  
  for (size_t i = 0; i < N; i++)
  {
    dydt[i] = 0.0;
  }
  
  /******************************/
  /* 3) Apply reaction rules    */
  /******************************/
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3.1) APPLY GENETIC REGULATION REACTIONS    */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  size_t enhancer_TF_index     = 0;
  size_t operator_TF_index     = 0;
  size_t regulated_genes_index = 0;
  for (size_t region = 0; region < list->grn_N; region++)
  {
    /*---------------------------------------------------------------------------------------*/
    /* 4.1.1) Compute Ei(t) and Oi(t)                                                        */
    /*---------------------------------------------------------------------------------------*/
    
    /* A) Compute Ei(t) -> enhancer contribution to the transcriptional activity */
    
    double Ei = 0.0;
    
    /* For each TF binding to the enhancer site */
    for (size_t pos = 0; pos < Nenhancer[region]; pos++)
    {
      /* get TF concentration */
      double TF_conc = y[enhancer_TF_list[enhancer_TF_index]];
      
      if (TF_conc > 0.0)
      {
        /* if co-enzyme activity is activated, compute its contribution */
        if (coe_activity)
        {
          double coE_conc = 0.0;
          if (enhancer_coe_list[enhancer_TF_index] <= (int)Ncell)
          {
            coE_conc = y[enhancer_coe_list[enhancer_TF_index]+Ngenome+Ninherited-1];
          }
          if (enhancer_coe_type[enhancer_TF_index] == REPRESSOR)
          {
            TF_conc = TF_conc*(enhancer_coe_km[enhancer_TF_index]/(enhancer_coe_km[enhancer_TF_index]+coE_conc));
          }
          else if (enhancer_coe_type[enhancer_TF_index] == ACTIVATOR)
          {
            TF_conc = TF_conc*(coE_conc/(enhancer_coe_km[enhancer_TF_index]+coE_conc));
          }
          else if (enhancer_coe_type[enhancer_TF_index] == ALWAYS_REPRESSED)
          {
            TF_conc = 0.0;
          }
        }
        Ei += TF_conc*enhancer_affinity_list[enhancer_TF_index];
      }
      
      /* Increment enhancer TF index */
      enhancer_TF_index++;
    }
    
    /* B) Compute Oi(t) -> operator contribution to the transcriptional activity */
    
    double Oi = 0.0;
    
    /* For each TF binding to the enhancer site */
    for (size_t pos = 0; pos < Noperator[region]; pos++)
    {
      /* get TF concentration */
      double TF_conc = y[operator_TF_list[operator_TF_index]];
      
      if (TF_conc > 0.0)
      {
        /* if co-enzyme activity is activated, compute its contribution */
        if (coe_activity)
        {
          double coE_conc = 0.0;
          if (operator_coe_list[operator_TF_index] <= (int)Ncell)
          {
            coE_conc = y[operator_coe_list[operator_TF_index]+Ngenome+Ninherited-1];
          }
          if (operator_coe_type[operator_TF_index] == REPRESSOR)
          {
            TF_conc = TF_conc*(operator_coe_km[operator_TF_index]/(operator_coe_km[operator_TF_index]+coE_conc));
          }
          else if (operator_coe_type[operator_TF_index] == ACTIVATOR)
          {
            TF_conc = TF_conc*(coE_conc/(operator_coe_km[operator_TF_index]+coE_conc));
          }
          else if (operator_coe_type[operator_TF_index] == ALWAYS_REPRESSED)
          {
            TF_conc = 0.0;
          }
        }
        Oi += TF_conc*operator_affinity_list[operator_TF_index];
      }
      
      /* Increment operator TF index */
      operator_TF_index++;
    }
    
    /*---------------------------------------------------------------------------------------*/
    /* 4.1.2) Compute the new transcription rate ei(t)                                       */
    /*---------------------------------------------------------------------------------------*/
    double ei = 0.0;
    if (beta[region] > 0.0)
    {
      if (Ei > 0.0)
      {
        Ei = pow(Ei, hill_n);
      }
      if (Oi > 0.0)
      {
        Oi = pow(Oi, hill_n);
      }
      ei = beta[region] * theta_n / (Oi+theta_n) * (1.0 + (1.0/beta[region] - 1.0) * (Ei/(Ei+theta_n)));
    }
    
    /*---------------------------------------------------------------------------------------*/
    /* 4.1.3) For each regulated gene, compute the next concentration ei(t) - degr_rate*y(t) */
    /*---------------------------------------------------------------------------------------*/
    for (size_t pos = 0; pos < Ngenes[region]; pos++)
    {
      dydt[regulated_genes[regulated_genes_index+pos]] += (ei-y[regulated_genes[regulated_genes_index+pos]]*protein_degradation_rate)*timestep_ratio;
    }
    regulated_genes_index += Ngenes[region];
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4.2) DEGRADATE INHERITED PROTEINS          */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  if (enzymatic_inheritance)
  {
    for (size_t pos = 0; pos < Ninherited; pos++)
    {
      dydt[Ngenome+pos] -= y[Ngenome+pos]*protein_degradation_rate*timestep_ratio;
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4.3) APPLY METABOLIC REACTIONS             */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  for (size_t i = 0; i < list->metabolic_N; i++)
  {
    /*-----------------------------------------------------*/
    /* 4.3.1) If reaction is inflowing pump activity       */
    /*-----------------------------------------------------*/
    if (type[i] == INFLOWING_PUMP_ACTIVITY)
    {
      double dy = (kcat[i] * y[Ngenome+Ninherited+Ncell-1+s[i]] * y[e[i]])/(km[i] + y[Ngenome+Ninherited+Ncell-1+s[i]]);
      if (env_interaction)
      {
        dydt[Ngenome+Ninherited+Ncell-1+s[i]] -= dy;
      }
      dydt[Ngenome+Ninherited-1+p[i]] += dy;
      dydt[N-2]                       += dy;
    }
    
    /*-----------------------------------------------------*/
    /* 4.3.2) If reaction is outflowing pump activity      */
    /*-----------------------------------------------------*/
    else if (type[i] == OUTFLOWING_PUMP_ACTIVITY)
    {
      double dy = (kcat[i] * y[Ngenome+Ninherited-1+s[i]] * y[e[i]])/(km[i] + y[Ngenome+Ninherited-1+s[i]]);
      dydt[Ngenome+Ninherited-1+s[i]] -= dy;
      if (env_interaction)
      {
        dydt[Ngenome+Ninherited+Ncell-1+p[i]] += dy;
      }
      dydt[N-1] += dy;
    }
    
    /*-----------------------------------------------------*/
    /* 4.3.3) If reaction is inner cell catalytic activity */
    /*        consuming energy                             */
    /*-----------------------------------------------------*/
    else if (type[i] == CATALYTIC_CONSUMING_ACTIVITY)
    {
      double dy = (kcat[i] * y[Ngenome+Ninherited-1+s[i]] * y[e[i]])/(km[i] + y[Ngenome+Ninherited-1+s[i]]);
      dydt[Ngenome+Ninherited-1+s[i]] -= dy;
      dydt[Ngenome+Ninherited-1+p[i]] += dy;
    }
    
    /*-----------------------------------------------------*/
    /* 4.3.4) If reaction is inner cell catalytic activity */
    /*        producing energy                             */
    /*-----------------------------------------------------*/
    else if (type[i] == CATALYTIC_REWARDING_ACTIVITY)
    {
      double dy = (kcat[i] * y[Ngenome+Ninherited-1+s[i]] * y[e[i]])/(km[i] + y[Ngenome+Ninherited-1+s[i]]);
      dydt[Ngenome+Ninherited-1+s[i]] -= dy;
      dydt[Ngenome+Ninherited-1+p[i]] += dy;
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4.4) APPLY MEMBRANE PERMEABILITY REACTIONS */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  if (membrane_permeability && kmp > 0.0)
  {
    for (size_t i = 0; i < Ncell; i++)
    {
      dydt[Ngenome+Ninherited+i] += (y[Ngenome+Ninherited+Ncell+i] - y[Ngenome+Ninherited+i])*kmp;
      if (env_interaction)
      {
        dydt[Ngenome+Ninherited+Ncell+i] -= (y[Ngenome+Ninherited+Ncell+i] - y[Ngenome+Ninherited+i])*kmp;
      }
    }
  }
  
  return GSL_SUCCESS;
}

/**
 * \brief    Compute the ODE system with energy
 * \details  --
 * \param    double t
 * \param    const double y[]
 * \param    double dydt[]
 * \param    void* parameters
 * \return   \e int GSL_SUCCESS
 */
int ODE::ODE_system_with_energy( double t, const double y[], double dydt[], void* parameters )
{
  (void)t;
  
  /******************************/
  /* 1) Initialize variables    */
  /******************************/
  
  reaction_list*               list                      = (reaction_list*)parameters;
  size_t                       N                         = list->N;
  size_t                       Ngenome                   = list->Ngenome;
  size_t                       Ninherited                = list->Ninherited;
  size_t                       Ncell                     = list->Ncell;
  double                       kmp                       = list->kmp;
  std::vector<size_t>&         Nenhancer                 = list->grn_Nenhancer;
  std::vector<size_t>&         Noperator                 = list->grn_Noperator;
  std::vector<size_t>&         Ngenes                    = list->grn_Ngenes;
  std::vector<size_t>&         enhancer_TF_list          = list->grn_enhancer_TF_list;
  std::vector<double>&         enhancer_affinity_list    = list->grn_enhancer_affinity_list;
  std::vector<int>&            enhancer_coe_list         = list->grn_enhancer_coe_list;
  std::vector<co_enzyme_type>& enhancer_coe_type         = list->grn_enhancer_coe_type;
  std::vector<double>&         enhancer_coe_km           = list->grn_enhancer_coe_km;
  std::vector<size_t>&         operator_TF_list          = list->grn_operator_TF_list;
  std::vector<double>&         operator_affinity_list    = list->grn_operator_affinity_list;
  std::vector<int>&            operator_coe_list         = list->grn_operator_coe_list;
  std::vector<co_enzyme_type>& operator_coe_type         = list->grn_operator_coe_type;
  std::vector<double>&         operator_coe_km           = list->grn_operator_coe_km;
  std::vector<size_t>&         regulated_genes           = list->grn_regulated_genes;
  std::vector<double>&         beta                      = list->grn_beta;
  double                       hill_n                    = list->grn_hill_n;
  double                       theta_n                   = list->grn_theta_n;
  double                       protein_degradation_rate  = list->grn_protein_degradation_rate;
  double                       energy_transcription_cost = list->grn_energy_transcription_cost;
  double                       energy_degradation_cost   = list->grn_energy_degradation_cost;
  double                       energy_dissipation_rate   = list->grn_energy_dissipation_rate;
  double                       timestep_ratio            = list->timestep_ratio;
  bool                         coe_activity              = list->co_enzyme_activity;
  bool                         membrane_permeability     = list->membrane_permeability;
  bool                         enzymatic_inheritance     = list->enzymatic_inheritance;
  std::vector<reaction_type>&  type                      = list->metabolic_type;
  std::vector<int>&            s                         = list->metabolic_s;
  std::vector<int>&            p                         = list->metabolic_p;
  std::vector<double>&         km                        = list->metabolic_km;
  std::vector<double>&         kcat                      = list->metabolic_kcat;
  std::vector<double>&         delta_g                   = list->metabolic_delta_g;
  std::vector<size_t>&         e                         = list->metabolic_e;
  bool                         env_interaction           = list->interacts_with_environment;
  
  /******************************/
  /* 2) Initialize dy/dt vector */
  /******************************/
  
  for (size_t i = 0; i < N; i++)
  {
    dydt[i] = 0.0;
  }
  
  /******************************/
  /* 3) Apply reaction rules    */
  /******************************/
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3.1) APPLY GENETIC REGULATION REACTIONS    */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  size_t enhancer_TF_index     = 0;
  size_t operator_TF_index     = 0;
  size_t regulated_genes_index = 0;
  for (size_t region = 0; region < list->grn_N; region++)
  {
    /*---------------------------------------------------------------------------------------*/
    /* 4.1.1) Compute Ei(t) and Oi(t)                                                        */
    /*---------------------------------------------------------------------------------------*/
    
    /* A) Compute Ei(t) -> enhancer contribution to the transcriptional activity */
    
    double Ei = 0.0;
    
    /* For each TF binding to the enhancer site */
    for (size_t pos = 0; pos < Nenhancer[region]; pos++)
    {
      /* get TF concentration */
      double TF_conc = y[enhancer_TF_list[enhancer_TF_index]];
      
      if (TF_conc > 0.0)
      {
        /* if co-enzyme activity is activated, compute its contribution */
        if (coe_activity)
        {
          double coE_conc = 0.0;
          if (enhancer_coe_list[enhancer_TF_index] <= (int)Ncell)
          {
            coE_conc = y[enhancer_coe_list[enhancer_TF_index]+Ngenome+Ninherited-1];
          }
          if (enhancer_coe_type[enhancer_TF_index] == REPRESSOR)
          {
            TF_conc = TF_conc*(enhancer_coe_km[enhancer_TF_index]/(enhancer_coe_km[enhancer_TF_index]+coE_conc));
          }
          else if (enhancer_coe_type[enhancer_TF_index] == ACTIVATOR)
          {
            TF_conc = TF_conc*(coE_conc/(enhancer_coe_km[enhancer_TF_index]+coE_conc));
          }
          else if (enhancer_coe_type[enhancer_TF_index] == ALWAYS_REPRESSED)
          {
            TF_conc = 0.0;
          }
        }
        Ei += TF_conc*enhancer_affinity_list[enhancer_TF_index];
      }
      
      /* Increment enhancer TF index */
      enhancer_TF_index++;
    }
    
    /* B) Compute Oi(t) -> operator contribution to the transcriptional activity */
    
    double Oi = 0.0;
    
    /* For each TF binding to the enhancer site */
    for (size_t pos = 0; pos < Noperator[region]; pos++)
    {
      /* get TF concentration */
      double TF_conc = y[operator_TF_list[operator_TF_index]];
      
      if (TF_conc > 0.0)
      {
        /* if co-enzyme activity is activated, compute its contribution */
        if (coe_activity)
        {
          double coE_conc = 0.0;
          if (operator_coe_list[operator_TF_index] <= (int)Ncell)
          {
            coE_conc = y[operator_coe_list[operator_TF_index]+Ngenome+Ninherited-1];
          }
          if (operator_coe_type[operator_TF_index] == REPRESSOR)
          {
            TF_conc = TF_conc*(operator_coe_km[operator_TF_index]/(operator_coe_km[operator_TF_index]+coE_conc));
          }
          else if (operator_coe_type[operator_TF_index] == ACTIVATOR)
          {
            TF_conc = TF_conc*(coE_conc/(operator_coe_km[operator_TF_index]+coE_conc));
          }
          else if (operator_coe_type[operator_TF_index] == ALWAYS_REPRESSED)
          {
            TF_conc = 0.0;
          }
        }
        Oi += TF_conc*operator_affinity_list[operator_TF_index];
      }
      
      /* Increment operator TF index */
      operator_TF_index++;
    }
    
    /*---------------------------------------------------------------------------------------*/
    /* 4.1.2) Compute the new transcription rate ei(t)                                       */
    /*---------------------------------------------------------------------------------------*/
    double ei = 0.0;
    if (beta[region] > 0.0)
    {
      if (Ei > 0.0)
      {
        Ei = pow(Ei, hill_n);
      }
      if (Oi > 0.0)
      {
        Oi = pow(Oi, hill_n);
      }
      ei = beta[region] * theta_n / (Oi+theta_n) * (1.0 + (1.0/beta[region] - 1.0) * (Ei/(Ei+theta_n)));
    }
    
    /*---------------------------------------------------------------------------------------*/
    /* 4.1.3) For each regulated gene, compute the next concentration ei(t) - degr_rate*y(t) */
    /*---------------------------------------------------------------------------------------*/
    for (size_t pos = 0; pos < Ngenes[region]; pos++)
    {
      double prod = ei;
      double degr = y[regulated_genes[regulated_genes_index+pos]]*protein_degradation_rate;
      dydt[regulated_genes[regulated_genes_index+pos]] += (prod-degr)*timestep_ratio;
      dydt[N-3] -= (prod*energy_transcription_cost+degr*energy_degradation_cost)*timestep_ratio;
    }
    regulated_genes_index += Ngenes[region];
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4.2) DEGRADATE INHERITED PROTEINS          */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  if (enzymatic_inheritance)
  {
    for (size_t pos = 0; pos < Ninherited; pos++)
    {
      double degr = y[Ngenome+pos]*protein_degradation_rate;
      dydt[Ngenome+pos] -= degr*timestep_ratio;
      dydt[N-3] -= degr*energy_degradation_cost*timestep_ratio;
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4.3) APPLY METABOLIC REACTIONS             */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  for (size_t i = 0; i < list->metabolic_N; i++)
  {
    /*----------------------------------------------------------------------*/
    /* 4.3.1) if reaction is inflowing pump activity                        */
    /*----------------------------------------------------------------------*/
    if (type[i] == INFLOWING_PUMP_ACTIVITY)
    {
      double dy = (kcat[i] * y[Ngenome+Ninherited+Ncell+s[i]-1] * y[e[i]])/(km[i] + y[Ngenome+Ninherited+Ncell+s[i]-1]);
      if (env_interaction)
      {
        dydt[Ngenome+Ninherited+Ncell+s[i]-1] -= dy;
      }
      dydt[Ngenome+Ninherited+p[i]-1] += dy;
      dydt[N-3] -= dy*delta_g[i];
      dydt[N-2] += dy;
    }
    
    /*----------------------------------------------------------------------*/
    /* 4.3.2) if reaction is outflowing pump activity                       */
    /*----------------------------------------------------------------------*/
    else if (type[i] == OUTFLOWING_PUMP_ACTIVITY)
    {
      double dy = (kcat[i] * y[Ngenome+Ninherited+s[i]-1] * y[e[i]])/(km[i] + y[Ngenome+Ninherited+s[i]-1]);
      dydt[Ngenome+Ninherited+s[i]-1] -= dy;
      if (env_interaction)
      {
        dydt[Ngenome+Ninherited+Ncell+p[i]-1] += dy;
      }
      dydt[N-3] -= dy*delta_g[i];
      dydt[N-1] += dy;
    }
    
    /*----------------------------------------------------------------------*/
    /* 4.3.3) if reaction is inner cell catalytic activity consuming energy */
    /*----------------------------------------------------------------------*/
    else if (type[i] == CATALYTIC_CONSUMING_ACTIVITY)
    {
      double dy = (kcat[i] * y[Ngenome+Ninherited+s[i]-1] * y[e[i]])/(km[i] + y[Ngenome+Ninherited+s[i]-1]);
      dydt[Ngenome+Ninherited+s[i]-1] -= dy;
      dydt[Ngenome+Ninherited+p[i]-1] += dy;
      dydt[N-3] -= dy*delta_g[i];
    }
    
    /*----------------------------------------------------------------------*/
    /* 4.3.4) if reaction is inner cell catalytic activity rewarding energy */
    /*----------------------------------------------------------------------*/
    else if (type[i] == CATALYTIC_REWARDING_ACTIVITY)
    {
      double dy = (kcat[i] * y[Ngenome+Ninherited+s[i]-1] * y[e[i]])/(km[i] + y[Ngenome+Ninherited+s[i]-1]);
      dydt[Ngenome+Ninherited+s[i]-1] -= dy;
      dydt[Ngenome+Ninherited+p[i]-1] += dy;
      dydt[N-3] += dy*delta_g[i];
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4.4) APPLY ENERGY DISSIPATION REACTION     */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  dydt[N-3] -= y[N-3]*energy_dissipation_rate;
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4.5) APPLY MEMBRANE PERMEABILITY REACTIONS */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  if (membrane_permeability && kmp > 0.0)
  {
    for (size_t i = 0; i < Ncell; i++)
    {
      dydt[Ngenome+Ninherited+i] += (y[Ngenome+Ninherited+Ncell+i] - y[Ngenome+Ninherited+i])*kmp;
      if (env_interaction)
      {
        dydt[Ngenome+Ninherited+Ncell+i] -= (y[Ngenome+Ninherited+Ncell+i] - y[Ngenome+Ninherited+i])*kmp;
      }
    }
  }
  
  return GSL_SUCCESS;
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Build gsl system
 * \details  Build GSL ode solver system. State vector size is set at 1 at initialization, and is updated at each timestep
 * \param    void
 * \return   \e void
 */
void ODE::build_system( void )
{
  if (!_parameters->get_energy_costs_scheme())
  {
    _system = {ODE::ODE_system_without_energy, NULL, static_cast<size_t>(1), _reaction_list};
  }
  else
  {
    _system = {ODE::ODE_system_with_energy, NULL, static_cast<size_t>(1), _reaction_list};
  }
  _control = gsl_odeiv2_control_y_new(ERR_ABS, ERR_REL);
}

/**
 * \brief    Build state vector X
 * \details  The index of each species and enzymatic concentrations are saved as is, because a strict identity between ODE state vector indexes and units or species indexes is maintained. Hence, the state vector X is :
 [unit(0), unit(1), ..., unit(Ngenome-1) | unit(0), unit(1), ..., unit(Ninherited-1) | tag(0), tag(1), ..., tag(Ncell-1) | env_tag(0), env_tag(1), ..., env_tag(Nenv-1) | energy | metabolic uptake | metabolic release]
 with Ngenome the length of the genome, Ninherited the length of inherited proteins, Ncell the length of the species list, Ncell the length of the local environment species list.
 Ncell and Nenv should be equal. The length of X is then N = Ngenome+Ninherited+Ncell+Nenv+3.
 The first part of the state vector X contains enzymes and TF concentrations.
 * \param    void
 * \return   \e void
 */
void ODE::build_X( void )
{
  /*------------------------------------------*/
  /* 1) allocate memory for the state vectors */
  /*------------------------------------------*/
  assert(_reaction_list->N > 0);
  _X = new double[_reaction_list->N];
  
  /*------------------------------------------*/
  /* 2) fill vectors with the initial values  */
  /*------------------------------------------*/
  size_t Ngenome    = _reaction_list->Ngenome;
  size_t Ninherited = _reaction_list->Ninherited;
  size_t Ncell      = _reaction_list->Ncell;
  
  /* 2.1) get enzymatic concentrations ------------*/
  for (size_t i = 0; i < Ngenome; i++)
  {
    _X[i] = (_genome->get_concentration_vector()[i] > MINIMUM_CONCENTRATION ? _genome->get_concentration_vector()[i] : 0.0);
  }
  if (_parameters->get_enzymatic_inheritance_scheme())
  {
    for (size_t i = 0; i < Ninherited; i++)
    {
      _X[Ngenome+i] = (_inherited_proteins->get_concentration_vector()[i] > MINIMUM_CONCENTRATION ? _inherited_proteins->get_concentration_vector()[i] : 0.0);
    }
  }
  
  /* 2.2) get metabolic concentrations ------------*/
  for (size_t i = 0; i < Ncell; i++)
  {
    _X[Ngenome+Ninherited+i]       = (_species_list->get_X()[i] > MINIMUM_CONCENTRATION ? _species_list->get_X()[i] : 0.0);
    _X[Ngenome+Ninherited+Ncell+i] = (_environment->get_species_list(_x, _y)->get_X()[i] > MINIMUM_CONCENTRATION ? _environment->get_species_list(_x, _y)->get_X()[i] : 0.0);
  }
  
  /* 2.2) initialize other cell variables ---------*/
  _X[_reaction_list->N-3] = (*_energy);
  _X[_reaction_list->N-2] = 0.0;
  _X[_reaction_list->N-1] = 0.0;
}

/**
 * \brief    Evaluate the promoter region
 * \details  Evaluate if the promoter region is functional (returns true if functional, false else)
 * \param    size_t promoter_position
 * \param    std::vector<size_t>* expressed_genes
 * \return   \e void
 */
bool ODE::evaluate_promoter_region( size_t promoter_position, std::vector<size_t>* expressed_genes )
{
  size_t        genome_size      = _genome->get_size();
  genetic_unit* str              = _genome->get_genetic_sequence()->x;
  bool          genes_found      = false;
  size_t        current_position = (promoter_position+1+genome_size)%genome_size;
  bool          valid_region     = false;
  expressed_genes->clear();
  
  /*--------------------------------------------------------*/
  /* While no stop criteria are encountered (P or NC units) */
  /*--------------------------------------------------------*/
  while (str[current_position].type != PROMOTER && str[current_position].type != NON_CODING)
  {
    /*------------------------------------------------------------------------------*/
    /* A) If a E or TF type is met for the first time, start regulated genes saving */
    /*------------------------------------------------------------------------------*/
    if ((str[current_position].type == ENZYME || str[current_position].type == TRANSCRIPTION_FACTOR) && !genes_found)
    {
      valid_region = true;
      genes_found  = true;
      expressed_genes->push_back(current_position);
    }
    /*------------------------------------------------------------------------------*/
    /* B) Save all (E-TF) types as regulated genes                                  */
    /*------------------------------------------------------------------------------*/
    else if ((str[current_position].type == ENZYME || str[current_position].type == TRANSCRIPTION_FACTOR) && genes_found)
    {
      expressed_genes->push_back(current_position);
    }
    /*------------------------------------------------------------------------------*/
    /* C) At the end of the regulated zone, stop the exploration                    */
    /*------------------------------------------------------------------------------*/
    else if ((str[current_position].type == BINDING_SITE) && genes_found)
    {
      break;
    }
    current_position = (current_position+1+genome_size)%genome_size;
  }
  return valid_region;
}

/**
 * \brief    Add a functional region to the reaction list
 * \details  Add all the reaction rules (regulation and metabolic network) related to one functional region.
 * \param    size_t promoter_position
 * \param    std::vector<size_t>* expressed_genes
 * \param    std::unordered_map<size_t, unsigned char>& nodes_list
 * \return   \e void
 */
void ODE::add_functional_region( size_t promoter_position, std::vector<size_t>* expressed_genes, std::unordered_map<size_t, unsigned char>& nodes_list )
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Initialize some useful variables                     */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  /* 1.1) Get genome attributes ------------*/
  size_t        genome_size = _genome->get_size();
  genetic_unit* genome_str  = _genome->get_genetic_sequence()->x;
  size_t*       genome_TFi  = _genome->get_TFi();
  
  /* 1.2) Get inherited proteins attributes */
  genetic_unit* inherited_str = NULL;
  size_t*       inherited_TFi = NULL;
  if (_parameters->get_enzymatic_inheritance_scheme())
  {
    inherited_str = _inherited_proteins->get_genetic_sequence()->x;
    inherited_TFi = _inherited_proteins->get_TFi();
  }
  
  /* 1.3) Get reaction list length ---------*/
  size_t N = _reaction_list->grn_N;
  
  /* 1.4) Set the promoter functional ------*/
  genome_str[promoter_position].functional = true;
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Evaluate the enhancer site                           */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  _reaction_list->grn_Nenhancer.push_back(0);
  _reaction_list->grn_enhancer_nb_BS.push_back(0);
  size_t current_position = (promoter_position-1+genome_size)%genome_size;
  
  /* FOR EACH BINDING SITE PRECEDING THE PROMOTER */
  while (genome_str[current_position].type == BINDING_SITE)
  {
    genome_str[current_position].functional = true;
    _reaction_list->grn_enhancer_nb_BS[N]++;
    int TF_tag = genome_str[current_position].TF_tag;
    
    /*--------------------------------------------------------------*/
    /* 2.A) For each transcription factor in the genome             */
    /*--------------------------------------------------------------*/
    
    for (size_t index = 0; index < _genome->get_nb_TF(); index++)
    {
      size_t TF_pos = genome_TFi[index];
      
      /* if the transcription factor bind to the binding site, save it */
      if (abs(genome_str[TF_pos].BS_tag-TF_tag) <= (int)genome_str[TF_pos].binding_window)
      {
        double affinity = 1.0;
        if (genome_str[TF_pos].binding_window > 0)
        {
          affinity = 1.0-(double)abs(genome_str[TF_pos].BS_tag-TF_tag)/(double)genome_str[TF_pos].binding_window;
        }
        _reaction_list->grn_enhancer_TF_list.push_back(TF_pos);
        _reaction_list->grn_enhancer_affinity_list.push_back(affinity);
        _reaction_list->grn_enhancer_coe_list.push_back(genome_str[TF_pos].coE_tag);
        _reaction_list->grn_enhancer_coe_km.push_back(fabs(genome_str[TF_pos].kcat/genome_str[TF_pos].kcat_km_ratio));
        if (genome_str[TF_pos].free_activity && !genome_str[TF_pos].bound_activity)
        {
          _reaction_list->grn_enhancer_coe_type.push_back(REPRESSOR);
        }
        else if (!genome_str[TF_pos].free_activity && genome_str[TF_pos].bound_activity)
        {
          _reaction_list->grn_enhancer_coe_type.push_back(ACTIVATOR);
        }
        else if (genome_str[TF_pos].free_activity && genome_str[TF_pos].bound_activity)
        {
          _reaction_list->grn_enhancer_coe_type.push_back(ALWAYS_ACTIVE);
        }
        else if (!genome_str[TF_pos].free_activity && !genome_str[TF_pos].bound_activity)
        {
          _reaction_list->grn_enhancer_coe_type.push_back(ALWAYS_REPRESSED);
        }
        
        /* Increment the number of equations related to this enhancer site */
        _reaction_list->grn_Nenhancer[N]++;
      }
    }
    
    /*--------------------------------------------------------------*/
    /* 2.B) For each transcription factor in the inherited proteins */
    /*--------------------------------------------------------------*/
    
    if (_parameters->get_enzymatic_inheritance_scheme())
    {
      for (size_t index = 0; index < _inherited_proteins->get_nb_TF(); index++)
      {
        size_t TF_pos = inherited_TFi[index];
        
        /* if the transcription factor bind to the binding site, save it */
        if (abs(inherited_str[TF_pos].BS_tag-TF_tag) <= (int)inherited_str[TF_pos].binding_window)
        {
          double affinity = 1.0;
          if (inherited_str[TF_pos].binding_window > 0)
          {
            affinity = 1.0-(double)abs(inherited_str[TF_pos].BS_tag-TF_tag)/(double)inherited_str[TF_pos].binding_window;
          }
          _reaction_list->grn_enhancer_TF_list.push_back(TF_pos+genome_size); /* [genome | inherited | ...] */
          _reaction_list->grn_enhancer_affinity_list.push_back(affinity);
          _reaction_list->grn_enhancer_coe_list.push_back(inherited_str[TF_pos].coE_tag);
          _reaction_list->grn_enhancer_coe_km.push_back(fabs(inherited_str[TF_pos].kcat/inherited_str[TF_pos].kcat_km_ratio));
          if (inherited_str[TF_pos].free_activity && !inherited_str[TF_pos].bound_activity)
          {
            _reaction_list->grn_enhancer_coe_type.push_back(REPRESSOR);
          }
          else if (!inherited_str[TF_pos].free_activity && inherited_str[TF_pos].bound_activity)
          {
            _reaction_list->grn_enhancer_coe_type.push_back(ACTIVATOR);
          }
          else if (inherited_str[TF_pos].free_activity && inherited_str[TF_pos].bound_activity)
          {
            _reaction_list->grn_enhancer_coe_type.push_back(ALWAYS_ACTIVE);
          }
          else if (!inherited_str[TF_pos].free_activity && !inherited_str[TF_pos].bound_activity)
          {
            _reaction_list->grn_enhancer_coe_type.push_back(ALWAYS_REPRESSED);
          }
          
          /* Increment the number of equations related to this enhancer site */
          _reaction_list->grn_Nenhancer[N]++;
        }
      }
    }
    
    /*--------------------------------------------------------------*/
    /* 2.C) Go to the previous position in the genome               */
    /*--------------------------------------------------------------*/
    
    current_position = (current_position-1+genome_size)%genome_size;
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) Evaluate the operator site                           */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  _reaction_list->grn_Noperator.push_back(0);
  _reaction_list->grn_operator_nb_BS.push_back(0);
  current_position = (promoter_position+1+genome_size)%genome_size;
  
  /* FOR EACH BINDING SITE FOLLOWING THE PROMOTER */
  while (genome_str[current_position].type == BINDING_SITE)
  {
    genome_str[current_position].functional = true;
    _reaction_list->grn_operator_nb_BS[N]++;
    int TF_tag = genome_str[current_position].TF_tag;
    
    /*--------------------------------------------------------------*/
    /* 3.A) For each transcription factor in the genome             */
    /*--------------------------------------------------------------*/
    
    for (size_t index = 0; index < _genome->get_nb_TF(); index++)
    {
      size_t TF_pos = genome_TFi[index];
      
      /* If the transcription factor bind to the binding site, save it */
      if (abs(genome_str[TF_pos].BS_tag-TF_tag) <= (int)genome_str[TF_pos].binding_window)
      {
        double affinity = 1.0;
        if (genome_str[TF_pos].binding_window > 0)
        {
          affinity = 1.0-(double)abs(genome_str[TF_pos].BS_tag-TF_tag)/(double)genome_str[TF_pos].binding_window;
        }
        _reaction_list->grn_operator_TF_list.push_back(TF_pos);
        _reaction_list->grn_operator_affinity_list.push_back(affinity);
        _reaction_list->grn_operator_coe_list.push_back(genome_str[TF_pos].coE_tag);
        _reaction_list->grn_operator_coe_km.push_back(fabs(genome_str[TF_pos].kcat/genome_str[TF_pos].kcat_km_ratio));
        if (genome_str[TF_pos].free_activity && !genome_str[TF_pos].bound_activity)
        {
          _reaction_list->grn_operator_coe_type.push_back(REPRESSOR);
        }
        else if (!genome_str[TF_pos].free_activity && genome_str[TF_pos].bound_activity)
        {
          _reaction_list->grn_operator_coe_type.push_back(ACTIVATOR);
        }
        else if (genome_str[TF_pos].free_activity && genome_str[TF_pos].bound_activity)
        {
          _reaction_list->grn_operator_coe_type.push_back(ALWAYS_ACTIVE);
        }
        else if (!genome_str[TF_pos].free_activity && !genome_str[TF_pos].bound_activity)
        {
          _reaction_list->grn_operator_coe_type.push_back(ALWAYS_REPRESSED);
        }
        
        /* Increment the number of equations related to this enhancer site */
        _reaction_list->grn_Noperator[N]++;
      }
    }
    
    /*--------------------------------------------------------------*/
    /* 3.B) For each transcription factor in the inherited proteins */
    /*--------------------------------------------------------------*/
    
    if (_parameters->get_enzymatic_inheritance_scheme())
    {
      for (size_t index = 0; index < _inherited_proteins->get_nb_TF(); index++)
      {
        size_t TF_pos = inherited_TFi[index];
        
        /* If the transcription factor bind to the binding site, save it */
        if (abs(inherited_str[TF_pos].BS_tag-TF_tag) <= (int)inherited_str[TF_pos].binding_window)
        {
          double affinity = 1.0;
          if (inherited_str[TF_pos].binding_window > 0)
          {
            affinity = 1.0-(double)abs(inherited_str[TF_pos].BS_tag-TF_tag)/(double)inherited_str[TF_pos].binding_window;
          }
          _reaction_list->grn_operator_TF_list.push_back(TF_pos+genome_size); /* [genome | inherited | ...] */
          _reaction_list->grn_operator_affinity_list.push_back(affinity);
          _reaction_list->grn_operator_coe_list.push_back(inherited_str[TF_pos].coE_tag);
          _reaction_list->grn_operator_coe_km.push_back(fabs(inherited_str[TF_pos].kcat/inherited_str[TF_pos].kcat_km_ratio));
          if (inherited_str[TF_pos].free_activity && !inherited_str[TF_pos].bound_activity)
          {
            _reaction_list->grn_operator_coe_type.push_back(REPRESSOR);
          }
          else if (!inherited_str[TF_pos].free_activity && inherited_str[TF_pos].bound_activity)
          {
            _reaction_list->grn_operator_coe_type.push_back(ACTIVATOR);
          }
          else if (inherited_str[TF_pos].free_activity && inherited_str[TF_pos].bound_activity)
          {
            _reaction_list->grn_operator_coe_type.push_back(ALWAYS_ACTIVE);
          }
          else if (!inherited_str[TF_pos].free_activity && !inherited_str[TF_pos].bound_activity)
          {
            _reaction_list->grn_operator_coe_type.push_back(ALWAYS_REPRESSED);
          }
          
          /* Increment the number of equations related to this enhancer site */
          _reaction_list->grn_Noperator[N]++;
        }
      }
    }
    
    /*--------------------------------------------------------------*/
    /* 3.B) Go to the next position in the genome                   */
    /*--------------------------------------------------------------*/
    
    current_position = (current_position+1+genome_size)%genome_size;
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4) Evaluate regulated genes and add metabolic reactions */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  _reaction_list->grn_Ngenes.push_back(expressed_genes->size());
  for (size_t gene_index = 0; gene_index < _reaction_list->grn_Ngenes[N]; gene_index++)
  {
    size_t pos = expressed_genes->at(gene_index);
    genome_str[pos].functional = true;
    _reaction_list->grn_regulated_genes.push_back(pos);
    
    genetic_unit myunit = genome_str[pos];
    
    if (myunit.type == ENZYME)
    {
      /**************************/
      /* Case 4.A) Sext -> Sint */
      /**************************/
      if (myunit.s == myunit.p && myunit.kcat > 0.0)
      {
        add_metabolic_reaction(INFLOWING_PUMP_ACTIVITY, FROM_GENOME, myunit.s, myunit.p, myunit.kcat/myunit.kcat_km_ratio, myunit.kcat, _parameters->get_energy_pumping_cost(), pos);
        if (nodes_list.find(myunit.s) == nodes_list.end())
        {
          _inflowing_pumps->push_back(myunit.s);
          nodes_list[myunit.s] = 1;
          (*_metabolic_nb_nodes)++;
        }
        (*_metabolic_nb_edges)++;
      }
      /**************************/
      /* Case 4.B) Sint -> Sext */
      /**************************/
      else if (myunit.s == myunit.p && myunit.kcat < 0.0)
      {
        add_metabolic_reaction(OUTFLOWING_PUMP_ACTIVITY, FROM_GENOME, myunit.s, myunit.p, -myunit.kcat/myunit.kcat_km_ratio, -myunit.kcat, _parameters->get_energy_pumping_cost(), pos);
        if (nodes_list.find(myunit.s) == nodes_list.end())
        {
          _outflowing_pumps->push_back(myunit.s);
          nodes_list[myunit.s] = 1;
          (*_metabolic_nb_nodes)++;
        }
        (*_metabolic_nb_edges)++;
      }
      /**************************/
      /* Case 4.C) Sint -> Pint */
      /**************************/
      else if (myunit.s != myunit.p && fabs(myunit.s-myunit.p) <= _parameters->get_maximum_reaction_size() && myunit.kcat > 0.0)
      {
        if (myunit.s < myunit.p)
        {
          double delta_g = (myunit.p-myunit.s)*_parameters->get_energy_enzymatic_cost();
          assert(delta_g >= 0.0);
          add_metabolic_reaction(CATALYTIC_CONSUMING_ACTIVITY, FROM_GENOME, myunit.s, myunit.p, myunit.kcat/myunit.kcat_km_ratio, myunit.kcat, delta_g, pos);
          if (nodes_list.find(myunit.s) == nodes_list.end())
          {
            nodes_list[myunit.s] = 1;
            (*_metabolic_nb_nodes)++;
          }
          if (nodes_list.find(myunit.p) == nodes_list.end())
          {
            nodes_list[myunit.p] = 1;
            (*_metabolic_nb_nodes)++;
          }
          (*_metabolic_nb_edges)++;
        }
        else if (myunit.s > myunit.p)
        {
          double delta_g = (myunit.s-myunit.p)*_parameters->get_energy_enzymatic_cost();
          assert(delta_g >= 0.0);
          add_metabolic_reaction(CATALYTIC_REWARDING_ACTIVITY, FROM_GENOME, myunit.s, myunit.p, myunit.kcat/myunit.kcat_km_ratio, myunit.kcat, delta_g, pos);
          if (nodes_list.find(myunit.s) == nodes_list.end())
          {
            nodes_list[myunit.s] = 1;
            (*_metabolic_nb_nodes)++;
          }
          if (nodes_list.find(myunit.p) == nodes_list.end())
          {
            nodes_list[myunit.p] = 1;
            (*_metabolic_nb_nodes)++;
          }
          (*_metabolic_nb_edges)++;
        }
      }
      /**************************/
      /* Case 4.D) Pint -> Sint */
      /**************************/
      else if (myunit.s != myunit.p && fabs(myunit.s-myunit.p) <= _parameters->get_maximum_reaction_size() && myunit.kcat < 0.0)
      {
        if (myunit.p < myunit.s)
        {
          double delta_g = (myunit.s-myunit.p)*_parameters->get_energy_enzymatic_cost();
          assert(delta_g >= 0.0);
          add_metabolic_reaction(CATALYTIC_CONSUMING_ACTIVITY, FROM_GENOME, myunit.p, myunit.s, -myunit.kcat/myunit.kcat_km_ratio, -myunit.kcat, delta_g, pos);
          if (nodes_list.find(myunit.s) == nodes_list.end())
          {
            nodes_list[myunit.s] = 1;
            (*_metabolic_nb_nodes)++;
          }
          if (nodes_list.find(myunit.p) == nodes_list.end())
          {
            nodes_list[myunit.p] = 1;
            (*_metabolic_nb_nodes)++;
          }
          (*_metabolic_nb_edges)++;
        }
        else if (myunit.p > myunit.s)
        {
          double delta_g = (myunit.p-myunit.s)*_parameters->get_energy_enzymatic_cost();
          assert(delta_g >= 0.0);
          add_metabolic_reaction(CATALYTIC_REWARDING_ACTIVITY, FROM_GENOME, myunit.p, myunit.s, -myunit.kcat/myunit.kcat_km_ratio, -myunit.kcat, delta_g, pos);
          if (nodes_list.find(myunit.s) == nodes_list.end())
          {
            nodes_list[myunit.s] = 1;
            (*_metabolic_nb_nodes)++;
          }
          if (nodes_list.find(myunit.p) == nodes_list.end())
          {
            nodes_list[myunit.p] = 1;
            (*_metabolic_nb_nodes)++;
          }
          (*_metabolic_nb_edges)++;
        }
      }
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 5) Set promoter constants                               */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  _reaction_list->grn_promoter.push_back(promoter_position);
  _reaction_list->grn_beta.push_back(genome_str[promoter_position].basal_expression_level);
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 6) Increment the number of GRN reactions                */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  _reaction_list->grn_N++;
}

/**
 * \brief    Add a metabolic reaction in the reaction list
 * \details  --
 * \param    reaction_type type
 * \param    reaction_origin origin
 * \param    int s
 * \param    int p
 * \param    double km
 * \param    double kcat
 * \param    double delta_g
 * \param    size_t e
 * \return   \e void
 */
void ODE::add_metabolic_reaction( reaction_type type, reaction_origin origin, int s, int p, double km, double kcat, double delta_g, size_t e )
{
  _reaction_list->metabolic_type.push_back(type);
  _reaction_list->metabolic_origin.push_back(origin);
  _reaction_list->metabolic_s.push_back(s);
  _reaction_list->metabolic_p.push_back(p);
  _reaction_list->metabolic_km.push_back(km);
  _reaction_list->metabolic_kcat.push_back(kcat);
  _reaction_list->metabolic_delta_g.push_back(delta_g);
  _reaction_list->metabolic_e.push_back(e);
  _reaction_list->metabolic_N++;
}

/**
 * \brief    Add inherited proteins metabolic reactions
 * \details  --
 * \param    std::unordered_map<size_t, unsigned char>& nodes_list
 * \return   \e void
 */
void ODE::add_inherited_metabolic_reactions( std::unordered_map<size_t, unsigned char>& nodes_list )
{
  size_t        Ngenome = _genome->get_size();
  genetic_unit* str     = _inherited_proteins->get_genetic_sequence()->x;
  size_t*       Ei      = _inherited_proteins->get_Ei();
  /*---------------------------*/
  /* For each inherited enzyme */
  /*---------------------------*/
  for (size_t i = 0; i < _inherited_proteins->get_nb_E(); i++)
  {
    size_t pos = Ei[i];
    
    genetic_unit myunit = str[pos];
    
    /**************************/
    /* Case 4.A) Sext -> Sint */
    /**************************/
    if (myunit.s == myunit.p && myunit.kcat > 0.0)
    {
      add_metabolic_reaction(INFLOWING_PUMP_ACTIVITY, FROM_INHERITED_PROTEINS, myunit.s, myunit.p, myunit.kcat/myunit.kcat_km_ratio, myunit.kcat, _parameters->get_energy_pumping_cost(), pos+Ngenome);
      if (nodes_list.find(myunit.s) == nodes_list.end())
      {
        _inflowing_pumps->push_back(myunit.s);
        nodes_list[myunit.s] = 1;
        (*_metabolic_nb_nodes)++;
      }
      (*_metabolic_nb_edges)++;
    }
    /**************************/
    /* Case 4.B) Sint -> Sext */
    /**************************/
    else if (myunit.s == myunit.p && myunit.kcat < 0.0)
    {
      add_metabolic_reaction(OUTFLOWING_PUMP_ACTIVITY, FROM_INHERITED_PROTEINS, myunit.s, myunit.p, -myunit.kcat/myunit.kcat_km_ratio, -myunit.kcat, _parameters->get_energy_pumping_cost(), pos+Ngenome);
      if (nodes_list.find(myunit.s) == nodes_list.end())
      {
        _outflowing_pumps->push_back(myunit.s);
        nodes_list[myunit.s] = 1;
        (*_metabolic_nb_nodes)++;
      }
      (*_metabolic_nb_edges)++;
    }
    /**************************/
    /* Case 4.C) Sint -> Pint */
    /**************************/
    else if (myunit.s != myunit.p && fabs(myunit.s-myunit.p) <= _parameters->get_maximum_reaction_size() && myunit.kcat > 0.0)
    {
      if (myunit.s < myunit.p)
      {
        double delta_g = (myunit.p-myunit.s)*_parameters->get_energy_enzymatic_cost();
        assert(delta_g >= 0.0);
        add_metabolic_reaction(CATALYTIC_CONSUMING_ACTIVITY, FROM_INHERITED_PROTEINS, myunit.s, myunit.p, myunit.kcat/myunit.kcat_km_ratio, myunit.kcat, delta_g, pos+Ngenome);
        if (nodes_list.find(myunit.s) == nodes_list.end())
        {
          nodes_list[myunit.s] = 1;
          (*_metabolic_nb_nodes)++;
        }
        if (nodes_list.find(myunit.p) == nodes_list.end())
        {
          nodes_list[myunit.p] = 1;
          (*_metabolic_nb_nodes)++;
        }
        (*_metabolic_nb_edges)++;
      }
      else if (myunit.s > myunit.p)
      {
        double delta_g = (myunit.s-myunit.p)*_parameters->get_energy_enzymatic_cost();
        assert(delta_g >= 0.0);
        add_metabolic_reaction(CATALYTIC_REWARDING_ACTIVITY, FROM_INHERITED_PROTEINS, myunit.s, myunit.p, myunit.kcat/myunit.kcat_km_ratio, myunit.kcat, delta_g, pos+Ngenome);
        if (nodes_list.find(myunit.s) == nodes_list.end())
        {
          nodes_list[myunit.s] = 1;
          (*_metabolic_nb_nodes)++;
        }
        if (nodes_list.find(myunit.p) == nodes_list.end())
        {
          nodes_list[myunit.p] = 1;
          (*_metabolic_nb_nodes)++;
        }
        (*_metabolic_nb_edges)++;
      }
    }
    /**************************/
    /* Case 4.D) Pint -> Sint */
    /**************************/
    else if (myunit.s != myunit.p && fabs(myunit.s-myunit.p) <= _parameters->get_maximum_reaction_size() && myunit.kcat < 0.0)
    {
      if (myunit.p < myunit.s)
      {
        double delta_g = (myunit.s-myunit.p)*_parameters->get_energy_enzymatic_cost();
        assert(delta_g >= 0.0);
        add_metabolic_reaction(CATALYTIC_CONSUMING_ACTIVITY, FROM_INHERITED_PROTEINS, myunit.p, myunit.s, -myunit.kcat/myunit.kcat_km_ratio, -myunit.kcat, delta_g, pos+Ngenome);
        if (nodes_list.find(myunit.s) == nodes_list.end())
        {
          nodes_list[myunit.s] = 1;
          (*_metabolic_nb_nodes)++;
        }
        if (nodes_list.find(myunit.p) == nodes_list.end())
        {
          nodes_list[myunit.p] = 1;
          (*_metabolic_nb_nodes)++;
        }
        (*_metabolic_nb_edges)++;
      }
      else if (myunit.p > myunit.s)
      {
        double delta_g = (myunit.p-myunit.s)*_parameters->get_energy_enzymatic_cost();
        assert(delta_g >= 0.0);
        add_metabolic_reaction(CATALYTIC_REWARDING_ACTIVITY, FROM_INHERITED_PROTEINS, myunit.p, myunit.s, -myunit.kcat/myunit.kcat_km_ratio, -myunit.kcat, delta_g, pos+Ngenome);
        if (nodes_list.find(myunit.s) == nodes_list.end())
        {
          nodes_list[myunit.s] = 1;
          (*_metabolic_nb_nodes)++;
        }
        if (nodes_list.find(myunit.p) == nodes_list.end())
        {
          nodes_list[myunit.p] = 1;
          (*_metabolic_nb_nodes)++;
        }
        (*_metabolic_nb_edges)++;
      }
    }
  }
}

/**
 * \brief    Compute replication report statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void ODE::compute_replication_report_statistics( void )
{
  *(_grn_nb_nodes) = 0;
  *(_grn_nb_edges) = 0;
  std::unordered_map<size_t, char> grn_nodes;
  grn_nodes.clear();
  
  genetic_unit* genome_str    = _genome->get_genetic_sequence()->x;
  genetic_unit* inherited_str = NULL;
  if (_parameters->get_enzymatic_inheritance_scheme())
  {
    inherited_str = _inherited_proteins->get_genetic_sequence()->x;
  }
  
  /****************************************/
  /* 1) Initialize statistics             */
  /****************************************/
  _replication_report->set_nb_functional_regions(_reaction_list->grn_N);
  size_t nb_enhancers                = 0;
  size_t nb_operators                = 0;
  size_t nb_E_regions                = 0;
  size_t nb_TF_regions               = 0;
  size_t nb_mixed_regions            = 0;
  double mean_functional_region_size = 0.0;
  double mean_E_region_size          = 0.0;
  double mean_TF_region_size         = 0.0;
  double mean_mixed_region_size      = 0.0;
  double mean_enhancer_size          = 0.0;
  double mean_operator_size          = 0.0;
  double mean_operon_size            = 0.0;
  double mean_E_operon_size          = 0.0;
  double mean_TF_operon_size         = 0.0;
  double mean_mixed_operon_size      = 0.0;
  
  /****************************************/
  /* 2) For each functional region        */
  /****************************************/
  size_t enhancer_TF_index     = 0;
  size_t operator_TF_index     = 0;
  size_t regulated_genes_index = 0;
  for (size_t region = 0; region < _reaction_list->grn_N; region++)
  {
    /*-----------------------------*/
    /* 2.1) Evaluate enhancer data */
    /*-----------------------------*/
    if (_reaction_list->grn_Nenhancer[region] > 0)
    {
      /* Evaluate number of TFs */
      *(_grn_nb_edges) += _reaction_list->grn_Nenhancer[region]*_reaction_list->grn_Ngenes[region];
      nb_enhancers++;
      mean_enhancer_size += _reaction_list->grn_enhancer_nb_BS[region];
      
      /* Evaluate the nodes */
      for (size_t TFindex = 0; TFindex < _reaction_list->grn_Nenhancer[region]; TFindex++)
      {
        size_t TFtag = _reaction_list->grn_enhancer_TF_list[enhancer_TF_index];
        if (grn_nodes.find(TFtag) == grn_nodes.end())
        {
          grn_nodes[TFtag] = 1;
        }
        if (_parameters->get_co_enzyme_activity_scheme())
        {
          *(_grn_nb_edges) += 1;
          size_t COEtag = _reaction_list->grn_enhancer_coe_list[enhancer_TF_index]+_reaction_list->Ngenome+_reaction_list->Ninherited;
          if (grn_nodes.find(COEtag) == grn_nodes.end())
          {
            grn_nodes[COEtag] = 1;
          }
        }
        enhancer_TF_index++;
      }
    }
    
    /*-----------------------------*/
    /* 2.2) Evaluate operator data */
    /*-----------------------------*/
    if (_reaction_list->grn_Noperator[region] > 0)
    {
      /* Evaluate number of TFs */
      *(_grn_nb_edges) += _reaction_list->grn_Noperator[region]*_reaction_list->grn_Ngenes[region];
      nb_operators++;
      mean_operator_size += _reaction_list->grn_operator_nb_BS[region];
      
      /* Evaluate the nodes */
      for (size_t TFindex = 0; TFindex < _reaction_list->grn_Noperator[region]; TFindex++)
      {
        size_t TFtag = _reaction_list->grn_operator_TF_list[operator_TF_index];
        if (grn_nodes.find(TFtag) == grn_nodes.end())
        {
          grn_nodes[TFtag] = 1;
        }
        if (_parameters->get_co_enzyme_activity_scheme())
        {
          *(_grn_nb_edges) += 1;
          size_t COEtag = _reaction_list->grn_operator_coe_list[operator_TF_index]+_reaction_list->Ngenome+_reaction_list->Ninherited;
          if (grn_nodes.find(COEtag) == grn_nodes.end())
          {
            grn_nodes[COEtag] = 1;
          }
        }
        operator_TF_index++;
      }
    }
    
    /*-----------------------------*/
    /* 2.3) Evaluate operon data   */
    /*-----------------------------*/
    mean_operon_size += _reaction_list->grn_Ngenes[region];
    
    bool E_region  = false;
    bool TF_region = false;
    for (size_t i = 0; i < _reaction_list->grn_Ngenes[region]; i++)
    {
      if (genome_str[_reaction_list->grn_regulated_genes[regulated_genes_index]].type == ENZYME)
      {
        E_region = true;
      }
      else if (genome_str[_reaction_list->grn_regulated_genes[regulated_genes_index]].type == TRANSCRIPTION_FACTOR)
      {
        TF_region = true;
      }
      regulated_genes_index++;
    }
    mean_functional_region_size += _reaction_list->grn_enhancer_nb_BS[region]+_reaction_list->grn_operator_nb_BS[region]+_reaction_list->grn_Ngenes[region]+1;
    if (E_region && !TF_region)
    {
      nb_E_regions++;
      mean_E_operon_size += _reaction_list->grn_Ngenes[region];
      mean_E_region_size += _reaction_list->grn_enhancer_nb_BS[region]+_reaction_list->grn_operator_nb_BS[region]+_reaction_list->grn_Ngenes[region]+1;
    }
    else if (!E_region && TF_region)
    {
      nb_TF_regions++;
      mean_TF_operon_size += _reaction_list->grn_Ngenes[region];
      mean_TF_region_size += _reaction_list->grn_enhancer_nb_BS[region]+_reaction_list->grn_operator_nb_BS[region]+_reaction_list->grn_Ngenes[region]+1;
    }
    else if (E_region && TF_region)
    {
      nb_mixed_regions++;
      mean_mixed_operon_size += _reaction_list->grn_Ngenes[region];
      mean_mixed_region_size += _reaction_list->grn_enhancer_nb_BS[region]+_reaction_list->grn_operator_nb_BS[region]+_reaction_list->grn_Ngenes[region]+1;
    }
  }
  
  *(_grn_nb_nodes) = grn_nodes.size();
  grn_nodes.clear();
  
  /****************************************/
  /* 3) Compute mean values               */
  /****************************************/
  if (_reaction_list->grn_N > 0)
  {
    mean_functional_region_size /= _reaction_list->grn_N;
    mean_operon_size            /= _reaction_list->grn_N;
  }
  if (nb_E_regions > 0)
  {
    mean_E_region_size /= nb_E_regions;
    mean_E_operon_size /= nb_E_regions;
  }
  if (nb_TF_regions > 0)
  {
    mean_TF_region_size /= nb_TF_regions;
    mean_TF_operon_size /= nb_TF_regions;
  }
  if (nb_mixed_regions > 0)
  {
    mean_mixed_region_size /= nb_mixed_regions;
    mean_mixed_operon_size /= nb_mixed_regions;
  }
  if (nb_enhancers > 0)
  {
    mean_enhancer_size /= nb_enhancers;
  }
  if (nb_operators > 0)
  {
    mean_operator_size /= nb_operators;
  }
  
  /****************************************/
  /* 4) Save GRN data                     */
  /****************************************/
  _replication_report->set_nb_enhancers(nb_enhancers);
  _replication_report->set_nb_operators(nb_operators);
  _replication_report->set_nb_E_regions(nb_E_regions);
  _replication_report->set_nb_TF_regions(nb_TF_regions);
  _replication_report->set_nb_mixed_regions(nb_mixed_regions);
  _replication_report->set_mean_functional_region_size(mean_functional_region_size);
  _replication_report->set_mean_E_region_size(mean_E_region_size);
  _replication_report->set_mean_TF_region_size(mean_TF_region_size);
  _replication_report->set_mean_mixed_region_size(mean_mixed_region_size);
  _replication_report->set_mean_enhancer_size(mean_enhancer_size);
  _replication_report->set_mean_operator_size(mean_operator_size);
  _replication_report->set_mean_operon_size(mean_operon_size);
  _replication_report->set_mean_E_operon_size(mean_E_operon_size);
  _replication_report->set_mean_TF_operon_size(mean_TF_operon_size);
  _replication_report->set_mean_mixed_operon_size(mean_mixed_operon_size);
  
  /****************************************/
  /* 5) Compute the regulation redundancy */
  /****************************************/
  double mean_regulation_redundancy = 0.0;
  size_t number_of_reactions        = 0;
  
  /*### Explore the genome ###*/
  bool*   visited = new bool[_genome->get_nb_TF()];
  size_t* TFi     = _genome->get_TFi();
  for (size_t i = 0; i < _genome->get_nb_TF(); i++)
  {
    visited[i] = false;
  }
  for (size_t i = 0; i < _genome->get_nb_TF(); i++)
  {
    if (!visited[i])
    {
      genetic_unit* first_unit     = _genome->get_genetic_unit(TFi[i]);
      int           BS_tag         = first_unit->BS_tag;
      int           coE_tag        = first_unit->coE_tag;
      bool          free_activity  = first_unit->free_activity;
      bool          bound_activity = first_unit->bound_activity;
      
      size_t count = 1;
      number_of_reactions++;
      for (size_t j = i+1; j < _genome->get_nb_TF(); j++)
      {
        genetic_unit* second_unit = _genome->get_genetic_unit(TFi[j]);
        if (!visited[j] && BS_tag == second_unit->BS_tag && coE_tag == second_unit->coE_tag && free_activity == second_unit->free_activity && bound_activity == second_unit->bound_activity)
        {
          count++;
          visited[j] = true;
        }
      }
      mean_regulation_redundancy += count;
      visited[i] = true;
    }
  }
  delete[] visited;
  visited = NULL;
  
  /*### Explore the inherited proteins ###*/
  /*
  if (_parameters->get_enzymatic_inheritance_scheme())
  {
    visited = new bool[_inherited_proteins->get_nb_TF()];
    TFi     = _inherited_proteins->get_TFi();
    for (size_t i = 0; i < _inherited_proteins->get_nb_TF(); i++)
    {
      visited[i] = false;
    }
    for (size_t i = 0; i < _inherited_proteins->get_nb_TF(); i++)
    {
      if (!visited[i])
      {
        genetic_unit* first_unit     = _inherited_proteins->get_genetic_unit(TFi[i]);
        int           BS_tag         = first_unit->BS_tag;
        int           coE_tag        = first_unit->coE_tag;
        bool          free_activity  = first_unit->free_activity;
        bool          bound_activity = first_unit->bound_activity;
        
        size_t count = 1;
        number_of_reactions++;
        for (size_t j = i+1; j < _inherited_proteins->get_nb_TF(); j++)
        {
          genetic_unit* second_unit = _inherited_proteins->get_genetic_unit(TFi[j]);
          if (!visited[j] && BS_tag == second_unit->BS_tag && coE_tag == second_unit->coE_tag && free_activity == second_unit->free_activity && bound_activity == second_unit->bound_activity)
          {
            count++;
            visited[j] = true;
          }
        }
        mean_regulation_redundancy += count;
        visited[i] = true;
      }
    }
    delete[] visited;
    visited = NULL;
  }
  */
  
  if (number_of_reactions > 0)
  {
    mean_regulation_redundancy /= number_of_reactions;
  }
  _replication_report->set_mean_regulation_redundancy(mean_regulation_redundancy);
  
  /****************************************/
  /* 6) Compute the metabolic redundancy  */
  /****************************************/
  double mean_metabolic_redundancy = 0.0;
  number_of_reactions              = 0;
  visited                          = new bool[_reaction_list->metabolic_N];
  for (size_t i = 0; i < _reaction_list->metabolic_N; i++)
  {
    visited[i] = false;
  }
  for (size_t i = 0; i < _reaction_list->metabolic_N; i++)
  {
    if (_reaction_list->metabolic_origin[i] == FROM_GENOME && !visited[i])
    {
      int s        = _reaction_list->metabolic_s[i];
      int p        = _reaction_list->metabolic_p[i];
      double km    = _reaction_list->metabolic_km[i];
      double kcat  = _reaction_list->metabolic_kcat[i];
      size_t count = 1;
      number_of_reactions++;
      for (size_t j = i+1; j < _reaction_list->metabolic_N; j++)
      {
        if (!visited[j] && s == _reaction_list->metabolic_s[j] && p == _reaction_list->metabolic_p[j] && km == _reaction_list->metabolic_km[j] && kcat == _reaction_list->metabolic_kcat[j])
        {
          count++;
          visited[j] = true;
        }
      }
      mean_metabolic_redundancy += count;
      visited[i] = true;
    }
  }
  delete[] visited;
  visited = NULL;
  if (number_of_reactions > 0)
  {
    mean_metabolic_redundancy /= number_of_reactions;
  }
  _replication_report->set_mean_metabolic_redundancy(mean_metabolic_redundancy);
  
  /****************************************/
  /* 7) Compute functional size           */
  /****************************************/
  size_t functional_size = 0;
  for (size_t pos = 0; pos < _genome->get_size(); pos++)
  {
    if (genome_str[pos].functional)
    {
      functional_size++;
    }
  }
  _replication_report->set_genome_functional_size(functional_size);
  (void)inherited_str;
}

/**
 * \brief    Create the reaction list
 * \details  Depending if enzymatic inheritance is activated or not, buffer sizes vary
 * \param    void
 * \return   \e void
 */
void ODE::create_reaction_list( void )
{
  /*--------------------------------------------------------*/
  /* 1) Reserve reaction list memory                        */
  /*--------------------------------------------------------*/
  _reaction_list = new reaction_list;
  
  _reaction_list->grn_promoter.reserve(GRN_REACTIONS_BUFFER*sizeof(size_t));
  _reaction_list->grn_Nenhancer.reserve(GRN_REACTIONS_BUFFER*sizeof(size_t));
  _reaction_list->grn_Noperator.reserve(GRN_REACTIONS_BUFFER*sizeof(size_t));
  _reaction_list->grn_Ngenes.reserve(GRN_REACTIONS_BUFFER*sizeof(size_t));
  _reaction_list->grn_enhancer_nb_BS.reserve(GRN_REACTIONS_BUFFER*sizeof(size_t));
  _reaction_list->grn_enhancer_TF_list.reserve(GRN_REACTIONS_BUFFER*sizeof(size_t));
  _reaction_list->grn_enhancer_affinity_list.reserve(GRN_REACTIONS_BUFFER*sizeof(double));
  _reaction_list->grn_enhancer_coe_list.reserve(GRN_REACTIONS_BUFFER*sizeof(int));
  _reaction_list->grn_enhancer_coe_type.reserve(GRN_REACTIONS_BUFFER*sizeof(co_enzyme_type));
  _reaction_list->grn_enhancer_coe_km.reserve(GRN_REACTIONS_BUFFER*sizeof(double));
  _reaction_list->grn_operator_nb_BS.reserve(GRN_REACTIONS_BUFFER*sizeof(size_t));
  _reaction_list->grn_operator_TF_list.reserve(GRN_REACTIONS_BUFFER*sizeof(size_t));
  _reaction_list->grn_operator_affinity_list.reserve(GRN_REACTIONS_BUFFER*sizeof(double));
  _reaction_list->grn_operator_coe_list.reserve(GRN_REACTIONS_BUFFER*sizeof(int));
  _reaction_list->grn_operator_coe_type.reserve(GRN_REACTIONS_BUFFER*sizeof(co_enzyme_type));
  _reaction_list->grn_operator_coe_km.reserve(GRN_REACTIONS_BUFFER*sizeof(double));
  _reaction_list->grn_regulated_genes.reserve(GRN_REACTIONS_BUFFER*sizeof(size_t));
  _reaction_list->grn_beta.reserve(GRN_REACTIONS_BUFFER*sizeof(double));
  
  _reaction_list->metabolic_type.reserve(METABOLIC_REACTIONS_BUFFER*sizeof(reaction_type));
  _reaction_list->metabolic_origin.reserve(METABOLIC_REACTIONS_BUFFER*sizeof(reaction_origin));
  _reaction_list->metabolic_s.reserve(METABOLIC_REACTIONS_BUFFER*sizeof(int));
  _reaction_list->metabolic_p.reserve(METABOLIC_REACTIONS_BUFFER*sizeof(int));
  _reaction_list->metabolic_km.reserve(METABOLIC_REACTIONS_BUFFER*sizeof(double));
  _reaction_list->metabolic_kcat.reserve(METABOLIC_REACTIONS_BUFFER*sizeof(double));
  _reaction_list->metabolic_delta_g.reserve(METABOLIC_REACTIONS_BUFFER*sizeof(double));
  _reaction_list->metabolic_e.reserve(METABOLIC_REACTIONS_BUFFER*sizeof(size_t));
  
  /*--------------------------------------------------------*/
  /* 2) Initialize reaction counters and GRN parameters     */
  /*--------------------------------------------------------*/
  _reaction_list->grn_N       = 0;
  _reaction_list->metabolic_N = 0;
  
  _reaction_list->grn_hill_n                    = _parameters->get_hill_function_n();
  _reaction_list->grn_theta_n                   = pow(_parameters->get_hill_function_theta(), _parameters->get_hill_function_n());
  _reaction_list->grn_protein_degradation_rate  = _parameters->get_protein_degradation_rate();
  _reaction_list->grn_energy_transcription_cost = _parameters->get_energy_transcription_cost();
  _reaction_list->grn_energy_degradation_cost   = _parameters->get_energy_degradation_cost();
  _reaction_list->grn_energy_dissipation_rate   = _parameters->get_energy_dissipation_rate();
  
  /*--------------------------------------------------------*/
  /* 3) Initialize size of each state vector subspace       */
  /*--------------------------------------------------------*/
  _reaction_list->Ngenome    = 0;
  _reaction_list->Ninherited = 0;
  _reaction_list->Ncell      = 0;
  _reaction_list->Nenv       = 0;
  _reaction_list->N          = 0;
  
  /*--------------------------------------------------------*/
  /* 4) Load simulation schemes                             */
  /*--------------------------------------------------------*/
  _reaction_list->energy_constraints    = _parameters->get_energy_costs_scheme();
  _reaction_list->membrane_permeability = _parameters->get_membrane_permeability_scheme();
  _reaction_list->metabolic_inheritance = _parameters->get_metabolic_inheritance_scheme();
  _reaction_list->enzymatic_inheritance = _parameters->get_enzymatic_inheritance_scheme();
  _reaction_list->co_enzyme_activity    = _parameters->get_co_enzyme_activity_scheme();
  
  /*--------------------------------------------------------*/
  /* 5) Load membrane permeability constant                 */
  /*--------------------------------------------------------*/
  _reaction_list->kmp = _parameters->get_membrane_permeability();
  
  /*--------------------------------------------------------*/
  /* 6) Initialize metabolic exchanges with the environment */
  /*--------------------------------------------------------*/
  _reaction_list->metabolic_uptake  = 0.0;
  _reaction_list->metabolic_release = 0.0;
  
  /*--------------------------------------------------------*/
  /* 7) Load the timestep ratio                             */
  /*--------------------------------------------------------*/
  _reaction_list->timestep_ratio = _parameters->get_genetic_regulation_network_timestep()/_parameters->get_metabolism_timestep();
  
  /*--------------------------------------------------------*/
  /* 8) Load the environment interaction scheme             */
  /*--------------------------------------------------------*/
  _reaction_list->interacts_with_environment = false;
  if (_parameters->get_environment_properties()->interaction_scheme == INTERACTION)
  {
    _reaction_list->interacts_with_environment = true;
  }
  
  /*--------------------------------------------------------*/
  /* 9) Load test variables                                 */
  /*--------------------------------------------------------*/
#ifdef DEBUG
  _reaction_list->genome             = _genome;
  _reaction_list->inherited_proteins = _inherited_proteins;
  _reaction_list->cell_species_list  = _species_list;
  _reaction_list->env_species_list   = _environment->get_species_list(_x, _y);
#endif
}

/**
 * \brief    Free reaction list memory
 * \details  --
 * \param    void
 * \return   \e void
 */
void ODE::free_reaction_list( void )
{
  /*-------------------------------------------*/
  /* 1) Delete GRN reaction list               */
  /*-------------------------------------------*/
  _reaction_list->grn_promoter.clear();
  _reaction_list->grn_Nenhancer.clear();
  _reaction_list->grn_Noperator.clear();
  _reaction_list->grn_Ngenes.clear();
  _reaction_list->grn_enhancer_nb_BS.clear();
  _reaction_list->grn_enhancer_TF_list.clear();
  _reaction_list->grn_enhancer_affinity_list.clear();
  _reaction_list->grn_enhancer_coe_list.clear();
  _reaction_list->grn_enhancer_coe_type.clear();
  _reaction_list->grn_enhancer_coe_km.clear();
  _reaction_list->grn_operator_nb_BS.clear();
  _reaction_list->grn_operator_TF_list.clear();
  _reaction_list->grn_operator_affinity_list.clear();
  _reaction_list->grn_operator_coe_list.clear();
  _reaction_list->grn_operator_coe_type.clear();
  _reaction_list->grn_operator_coe_km.clear();
  _reaction_list->grn_regulated_genes.clear();
  _reaction_list->grn_beta.clear();
  
  /*-------------------------------------------*/
  /* 2) Delete metabolic network reaction list */
  /*-------------------------------------------*/
  _reaction_list->metabolic_type.clear();
  _reaction_list->metabolic_origin.clear();
  _reaction_list->metabolic_s.clear();
  _reaction_list->metabolic_p.clear();
  _reaction_list->metabolic_km.clear();
  _reaction_list->metabolic_kcat.clear();
  _reaction_list->metabolic_delta_g.clear();
  _reaction_list->metabolic_e.clear();
  
  /*-------------------------------------------*/
  /* 3) Delete reaction list and state vector  */
  /*-------------------------------------------*/
  delete _reaction_list;
  _reaction_list = NULL;
}

/**
 * \brief    Test the structure of the reaction list
 * \details  --
 * \param    void
 * \return   \e void
 */
#ifdef DEBUG
void ODE::test_reaction_list_structure( void )
{
  size_t N          = _reaction_list->N;
  size_t Ngenome    = _reaction_list->Ngenome;
  size_t Ninherited = _reaction_list->Ninherited;
  size_t Ncell      = _reaction_list->Ncell;

  Genome*            genome   = _reaction_list->genome;
  InheritedProteins* inhprot  = _reaction_list->inherited_proteins;
  SpeciesList*       cell_vec = _reaction_list->cell_species_list;
  SpeciesList*       env_vec  = _reaction_list->env_species_list;
  
  assert(cell_vec->get_size() == Ncell);
  assert(env_vec->get_size() == _reaction_list->Nenv);
  assert(cell_vec->get_size() == env_vec->get_size());
  
  // A) TEST GRN STRUCTURE :
  // -----------------------
  std::vector<size_t>& Nenhancer              = _reaction_list->grn_Nenhancer;
  std::vector<size_t>& Noperator              = _reaction_list->grn_Noperator;
  std::vector<size_t>& Ngenes                 = _reaction_list->grn_Ngenes;
  std::vector<size_t>& enhancer_TF_list       = _reaction_list->grn_enhancer_TF_list;
  std::vector<double>& enhancer_affinity_list = _reaction_list->grn_enhancer_affinity_list;
  std::vector<int>&    enhancer_coe_list      = _reaction_list->grn_enhancer_coe_list;
  std::vector<size_t>& operator_TF_list       = _reaction_list->grn_operator_TF_list;
  std::vector<double>& operator_affinity_list = _reaction_list->grn_operator_affinity_list;
  std::vector<int>&    operator_coe_list      = _reaction_list->grn_operator_coe_list;
  std::vector<size_t>& regulated_genes        = _reaction_list->grn_regulated_genes;
  std::vector<double>& beta                   = _reaction_list->grn_beta;
  bool                 coe_activity           = _reaction_list->co_enzyme_activity;
  bool                 enzymatic_inheritance  = _reaction_list->enzymatic_inheritance;
  
  size_t enhancer_TF_index     = 0;
  size_t operator_TF_index     = 0;
  size_t regulated_genes_index = 0;
  for (size_t region = 0; region < _reaction_list->grn_N; region++)
  {
    // Test enhancer region
    for (size_t TF_index = 0; TF_index < Nenhancer[region]; TF_index++)
    {
      // Test if TF index corresponds to a TF unit
      size_t TF_pos         = enhancer_TF_list[enhancer_TF_index];
      genetic_unit* TF_unit = NULL;
      if (TF_pos < Ngenome)
      {
        TF_unit = genome->get_genetic_unit(TF_pos);
      }
      else
      {
        TF_unit = inhprot->get_genetic_unit(TF_pos-Ngenome);
      }
      assert(TF_unit->type == TRANSCRIPTION_FACTOR);
      // Test if co-enzyme index is correct
      if (coe_activity)
      {
        assert(enhancer_coe_list[enhancer_TF_index] > 0);
        assert(enhancer_coe_list[enhancer_TF_index]+Ngenome+Ninherited-1 >= Ngenome+Ninherited);
        assert(enhancer_coe_list[enhancer_TF_index] == TF_unit->coE_tag);
      }
      // Test TF affinity values
      assert(enhancer_affinity_list[enhancer_TF_index] >= 0.0);
      assert(enhancer_affinity_list[enhancer_TF_index] <= 1.0);
      
      // Increment enhancer TF index
      enhancer_TF_index++;
    }
    
    // Test operator region
    for (size_t TF_index = 0; TF_index < Noperator[region]; TF_index++)
    {
      // Test if TF index corresponds to a TF unit
      size_t TF_pos         = operator_TF_list[operator_TF_index];
      genetic_unit* TF_unit = NULL;
      if (TF_pos < Ngenome)
      {
        TF_unit = genome->get_genetic_unit(TF_pos);
      }
      else
      {
        TF_unit = inhprot->get_genetic_unit(TF_pos-Ngenome);
      }
      assert(TF_unit->type == TRANSCRIPTION_FACTOR);
      // Test if co-enzyme index is correct
      if (coe_activity)
      {
        assert(operator_coe_list[operator_TF_index] > 0);
        assert(operator_coe_list[operator_TF_index]+Ngenome+Ninherited-1 >= Ngenome+Ninherited);
        assert(operator_coe_list[operator_TF_index] == TF_unit->coE_tag);
      }
      // Test TF affinity values
      assert(operator_affinity_list[operator_TF_index] >= 0.0);
      assert(operator_affinity_list[operator_TF_index] <= 1.0);
      
      // Increment operator TF index
      operator_TF_index++;
    }
    
    // Test promoter values
    genetic_unit* prom_unit = genome->get_genetic_unit(_reaction_list->grn_promoter[region]);
    assert(prom_unit->type == PROMOTER);
    assert(beta[region] == prom_unit->basal_expression_level);
    
    //Test regulated genes
    for (size_t gene_index = 0; gene_index < Ngenes[region]; gene_index++)
    {
      assert(regulated_genes[regulated_genes_index] < genome->get_size());
      genetic_unit* reg_unit = genome->get_genetic_unit(regulated_genes[regulated_genes_index]);
      assert(reg_unit->type == ENZYME || reg_unit->type == TRANSCRIPTION_FACTOR);
      regulated_genes_index++;
    }
  }
  
  // Test inherited proteins structure
  if (enzymatic_inheritance)
  {
    for (size_t index = 0; index < Ninherited; index++)
    {
      assert(index < inhprot->get_size());
      genetic_unit* inh_unit = inhprot->get_genetic_unit(index);
      assert(inh_unit->type == ENZYME || inh_unit->type == TRANSCRIPTION_FACTOR);
    }
  }
  
  // A) TEST METABOLIC NETWORK STRUCTURE :
  // -------------------------------------
  std::vector<reaction_type>& type = _reaction_list->metabolic_type;
  std::vector<int>&           s    = _reaction_list->metabolic_s;
  std::vector<int>&           p    = _reaction_list->metabolic_p;
  std::vector<double>&        km   = _reaction_list->metabolic_km;
  std::vector<double>&        kcat = _reaction_list->metabolic_kcat;
  std::vector<size_t>&        e    = _reaction_list->metabolic_e;
  
  for (size_t i = 0; i < _reaction_list->metabolic_N; i++)
  {
    genetic_unit* E_unit = NULL;
    if (e[i] < Ngenome)
    {
      E_unit = genome->get_genetic_unit(e[i]);
    }
    else
    {
      E_unit = inhprot->get_genetic_unit(e[i]-Ngenome);
    }
    if (type[i] == INFLOWING_PUMP_ACTIVITY)
    {
      assert(E_unit->type == ENZYME);
      assert(E_unit->s == s[i]);
      assert(E_unit->p == p[i]);
      assert(E_unit->kcat/E_unit->kcat_km_ratio == km[i]);
      assert(E_unit->kcat == kcat[i]);
      assert(s[i] > 0);
      assert(s[i] <= (int)cell_vec->get_size());
      assert(p[i] > 0);
      assert(p[i] <= (int)cell_vec->get_size());
      assert(km[i] > 0.0);
      assert(kcat[i] > 0.0);
      assert(Ngenome+Ninherited+Ncell+s[i]-1 >= Ngenome+Ninherited+Ncell);
      assert(Ngenome+Ninherited+Ncell+s[i]-1 < N-3);
      assert(Ngenome+Ninherited+p[i]-1 >= Ngenome+Ninherited);
      assert(Ngenome+Ninherited+p[i]-1 < Ngenome+Ninherited+Ncell);
    }
    else if (type[i] == OUTFLOWING_PUMP_ACTIVITY)
    {
      assert(E_unit->type == ENZYME);
      assert(E_unit->s == s[i]);
      assert(E_unit->p == p[i]);
      assert(E_unit->kcat/E_unit->kcat_km_ratio == -km[i]);
      assert(E_unit->kcat == -kcat[i]);
      assert(s[i] > 0);
      assert(s[i] <= (int)cell_vec->get_size());
      assert(p[i] > 0);
      assert(p[i] <= (int)cell_vec->get_size());
      assert(km[i] > 0.0);
      assert(kcat[i] > 0.0);
      assert(Ngenome+Ninherited+s[i]-1 >= Ngenome+Ninherited);
      assert(Ngenome+Ninherited+s[i]-1 < Ngenome+Ninherited+Ncell);
      assert(Ngenome+Ninherited+Ncell+p[i]-1 >= Ngenome+Ninherited+Ncell);
      assert(Ngenome+Ninherited+Ncell+p[i]-1 < N-3);
    }
    else if (type[i] == CATALYTIC_CONSUMING_ACTIVITY || type[i] == CATALYTIC_REWARDING_ACTIVITY)
    {
      assert(E_unit->type == ENZYME);
      assert(E_unit->s == s[i] || E_unit->s == p[i]);
      assert(E_unit->p == p[i] || E_unit->p == s[i]);
      assert(E_unit->kcat/E_unit->kcat_km_ratio == km[i] || E_unit->kcat/E_unit->kcat_km_ratio == -km[i]);
      assert(E_unit->kcat == kcat[i] || E_unit->kcat == -kcat[i]);
      assert(s[i] > 0);
      assert(s[i] <= (int)cell_vec->get_size());
      assert(p[i] > 0);
      assert(p[i] <= (int)cell_vec->get_size());
      assert(km[i] > 0.0);
      assert(kcat[i] > 0.0);
      assert(Ngenome+Ninherited+s[i]-1 >= Ngenome+Ninherited);
      assert(Ngenome+Ninherited+s[i]-1 < Ngenome+Ninherited+Ncell);
      assert(Ngenome+Ninherited+p[i]-1 >= Ngenome+Ninherited);
      assert(Ngenome+Ninherited+p[i]-1 < Ngenome+Ninherited+Ncell);
    }
  }
}
#endif
