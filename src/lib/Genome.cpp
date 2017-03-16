
/**
 * \file      Genome.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Genome class definition
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

#include "Genome.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters
 * \param    Prng* prng
 * \param    ReplicationReport* replication_report
 * \return   \e void
 */
Genome::Genome( Parameters* parameters, Prng* prng, ReplicationReport* replication_report )
{
  /*------------------------------------------------------------------ pseudorandom numbers generator */
  
  _prng = prng;
  
  /*------------------------------------------------------------------ simulation parameters */
  
  _parameters = parameters;
  
  /*------------------------------------------------------------------ genetic sequence */
  
  create_genetic_sequence();
  
  /*------------------------------------------------------------------ concentration vector */
  
  _concentration_vector = NULL;
  
  /*------------------------------------------------------------------ statistical data */
  
  _nb_NC            = 0;
  _nb_E             = 0;
  _nb_TF            = 0;
  _nb_BS            = 0;
  _nb_P             = 0;
  _nb_inner_enzymes = 0;
  _nb_inflow_pumps  = 0;
  _nb_outflow_pumps = 0;
  
  /*------------------------------------------------------------------ replication report */
  
  _replication_report = replication_report;
  
  /*------------------------------------------------------------------ genetic unit positions */
  
  _TFi = NULL;
  _Pi  = NULL;
}

/**
 * \brief    Constructor from backup file
 * \details  Load Genome class from backup file
 * \param    Parameters* parameters
 * \param    Prng* prng
 * \param    ReplicationReport* replication_report
 * \param    gzFile backup_file
 * \return   \e void
 */
Genome::Genome( Parameters* parameters, Prng* prng, ReplicationReport* replication_report, gzFile backup_file )
{
  /*------------------------------------------------------------------ pseudorandom numbers generator */
  
  _prng = prng;
  
  /*------------------------------------------------------------------ simulation parameters */
  
  _parameters = parameters;
  
  /*------------------------------------------------------------------ genetic sequence */
  
  load_genetic_sequence(backup_file);
  
  /*------------------------------------------------------------------ concentration vector */
  
  _concentration_vector = NULL;
  if (_genetic_sequence->size > 0)
  {
    _concentration_vector = new double[_genetic_sequence->size];
    for (size_t i = 0; i < _genetic_sequence->size; i++)
    {
      gzread( backup_file, &_concentration_vector[i], sizeof(_concentration_vector[i]) );
    }
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  gzread( backup_file, &_nb_NC,            sizeof(_nb_NC) );
  gzread( backup_file, &_nb_E,             sizeof(_nb_E) );
  gzread( backup_file, &_nb_TF,            sizeof(_nb_TF) );
  gzread( backup_file, &_nb_BS,            sizeof(_nb_BS) );
  gzread( backup_file, &_nb_P,             sizeof(_nb_P) );
  gzread( backup_file, &_nb_inner_enzymes, sizeof(_nb_inner_enzymes) );
  gzread( backup_file, &_nb_inflow_pumps,  sizeof(_nb_inflow_pumps) );
  gzread( backup_file, &_nb_outflow_pumps, sizeof(_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ replication report */
  
  _replication_report = replication_report;
  
  /*------------------------------------------------------------------ genetic unit positions */
  
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
    for (size_t i = 0; i < _nb_TF; i++)
    {
      gzread( backup_file, &_TFi[i], sizeof(_TFi[i]) );
    }
  }
  _Pi = NULL;
  if (_nb_P > 0)
  {
    _Pi = new size_t[_nb_P];
    for (size_t i = 0; i < _nb_P; i++)
    {
      gzread( backup_file, &_Pi[i], sizeof(_Pi[i]) );
    }
  }
}

/**
 * \brief    Copy constructor
 * \details  --
 * \param    const Genome& genome
 * \param    Prng* prng
 * \param    ReplicationReport* replication_report
 * \return   \e void
 */
Genome::Genome( const Genome& genome, Prng* prng, ReplicationReport* replication_report )
{
  /*------------------------------------------------------------------ pseudorandom numbers generator */
  
  _prng = prng;
  
  /*------------------------------------------------------------------ simulation parameters */
  
  _parameters = genome._parameters;
  
  /*------------------------------------------------------------------ genetic sequence */
  
  copy_genetic_sequence(genome._genetic_sequence);
  
  /*------------------------------------------------------------------ concentration vector */
  
  _concentration_vector = NULL;
  if (_genetic_sequence->size)
  {
    _concentration_vector = new double[_genetic_sequence->size];
    memcpy(_concentration_vector, genome._concentration_vector, sizeof(double)*_genetic_sequence->size);
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  _nb_NC            = genome._nb_NC;
  _nb_E             = genome._nb_E;
  _nb_TF            = genome._nb_TF;
  _nb_BS            = genome._nb_BS;
  _nb_P             = genome._nb_P;
  _nb_inner_enzymes = genome._nb_inner_enzymes;
  _nb_inflow_pumps  = genome._nb_inflow_pumps;
  _nb_outflow_pumps = genome._nb_outflow_pumps;
  
  /*------------------------------------------------------------------ replication report */
  
  _replication_report = replication_report;
  
  /*------------------------------------------------------------------ genetic unit positions */
  
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
    memcpy(_TFi, genome._TFi, sizeof(size_t)*_nb_TF);
  }
  _Pi = NULL;
  if (_nb_P > 0)
  {
    _Pi = new size_t[_nb_P];
    memcpy(_Pi, genome._Pi, sizeof(size_t)*_nb_P);
  }
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
Genome::~Genome( void )
{
  delete_genetic_sequence();
  delete[] _concentration_vector;
  _concentration_vector = NULL;
  delete[] _TFi;
  _TFi = NULL;
  delete[] _Pi;
  _Pi = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Initialize the concentration vector at zero (after mutation)
 * \details  --
 * \param    void
 * \return   \e void
 */
void Genome::initialize_concentration_vector( void )
{
  delete[] _concentration_vector;
  _concentration_vector = new double[_genetic_sequence->size];
  for (size_t i = 0; i < _genetic_sequence->size; i++)
  {
    _concentration_vector[i] = 0.0;
  }
}

/**
 * \brief    Mutate the genome
 * \details  Apply rearrangements, then apply point mutations
 * \param    const double* mutation_rates
 * \return   \e void
 */
void Genome::mutate( const double* mutation_rates )
{
  /*---------------------------------------------------------*/
  /* 1) Clear the replication report                         */
  /*---------------------------------------------------------*/
  //_replication_report->clear();
  
  /*---------------------------------------------------------*/
  /* 2) Do rearrangements                                    */
  /*---------------------------------------------------------*/
  do_rearrangements(mutation_rates);
  
  /*---------------------------------------------------------*/
  /* 3) Do horizontal gene transfer                          */
  /*---------------------------------------------------------*/
  if (_parameters->get_hgt_rate() > 0.0)
  {
    do_hgt(mutation_rates);
  }
  
  /*---------------------------------------------------------*/
  /* 4) Do point mutations                                   */
  /*---------------------------------------------------------*/
  do_point_mutations(mutation_rates);
  
  /*---------------------------------------------------------*/
  /* 5) Compute the indexes list for each genetic unit types */
  /*---------------------------------------------------------*/
  genetic_unit* str = _genetic_sequence->x;
  delete[] _TFi;
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
  }
  delete[] _Pi;
  _Pi = NULL;
  if (_nb_P > 0)
  {
    _Pi  = new size_t[_nb_P];
  }
  size_t TFcount = 0;
  size_t Pcount  = 0;
  for (size_t pos = 0; pos < _genetic_sequence->size; pos++)
  {
    if (str[pos].type == TRANSCRIPTION_FACTOR)
    {
      _TFi[TFcount] = pos;
      TFcount++;
    }
    else if (str[pos].type == PROMOTER)
    {
      _Pi[Pcount] = pos;
      Pcount++;
    }
    str[pos].functional = false;
  }
  
  /*---------------------------------------------------------*/
  /* 6) Initialize concentration vector                      */
  /*---------------------------------------------------------*/
  initialize_concentration_vector();
}

/**
 * \brief    Shuffle the genome
 * \details  Shuffle genome order
 * \param    void
 * \return   \e void
 */
void Genome::shuffle( void )
{
  for (size_t i = 0; i < _genetic_sequence->size; i++)
  {
    size_t       pos1 = (size_t)_prng->uniform(0, (int)_genetic_sequence->size-1);
    size_t       pos2 = (size_t)_prng->uniform(0, (int)_genetic_sequence->size-1);
    genetic_unit tmp  = _genetic_sequence->x[pos1];
    _genetic_sequence->x[pos1] = _genetic_sequence->x[pos2];
    _genetic_sequence->x[pos2] = tmp;
  }
}

/**
 * \brief    Save in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void Genome::save( gzFile backup_file )
{
  /*------------------------------------------------------------------ genetic sequence */
  
  save_genetic_sequence(backup_file);
  
  /*------------------------------------------------------------------ concentration vector */
  
  if (_genetic_sequence->size > 0)
  {
    for (size_t i = 0; i < _genetic_sequence->size; i++)
    {
      gzwrite( backup_file, &_concentration_vector[i], sizeof(_concentration_vector[i]) );
    }
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  gzwrite( backup_file, &_nb_NC,            sizeof(_nb_NC) );
  gzwrite( backup_file, &_nb_E,             sizeof(_nb_E) );
  gzwrite( backup_file, &_nb_TF,            sizeof(_nb_TF) );
  gzwrite( backup_file, &_nb_BS,            sizeof(_nb_BS) );
  gzwrite( backup_file, &_nb_P,             sizeof(_nb_P) );
  gzwrite( backup_file, &_nb_inner_enzymes, sizeof(_nb_inner_enzymes) );
  gzwrite( backup_file, &_nb_inflow_pumps,  sizeof(_nb_inflow_pumps) );
  gzwrite( backup_file, &_nb_outflow_pumps, sizeof(_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ genetic unit positions */
  
  if (_nb_TF > 0)
  {
    for (size_t i = 0; i < _nb_TF; i++)
    {
      gzwrite( backup_file, &_TFi[i], sizeof(_TFi[i]) );
    }
  }
  if (_nb_P > 0)
  {
    for (size_t i = 0; i < _nb_P; i++)
    {
      gzwrite( backup_file, &_Pi[i], sizeof(_Pi[i]) );
    }
  }
}

/**
 * \brief    Replace genome data
 * \details  --
 * \param    Genome* genome
 * \return   \e void
 */
void Genome::replace_data( Genome* genome )
{
  /*------------------------------------------------------------------ genetic sequence */
  
  delete_genetic_sequence();
  copy_genetic_sequence(genome->get_genetic_sequence());
  
  /*------------------------------------------------------------------ concentration vector */
  
  delete[] _concentration_vector;
  _concentration_vector = NULL;
  if (_genetic_sequence->size > 0)
  {
    _concentration_vector = new double[_genetic_sequence->size];
    memcpy(_concentration_vector, genome->get_concentration_vector(), sizeof(double)*_genetic_sequence->size);
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  _nb_NC            = genome->get_nb_NC();
  _nb_E             = genome->get_nb_E();
  _nb_TF            = genome->get_nb_TF();
  _nb_BS            = genome->get_nb_BS();
  _nb_P             = genome->get_nb_P();
  _nb_inner_enzymes = genome->get_nb_inner_enzymes();
  _nb_inflow_pumps  = genome->get_nb_inflow_pumps();
  _nb_outflow_pumps = genome->get_nb_outflow_pumps();
  
  /*------------------------------------------------------------------ genetic unit positions */
  
  delete[] _TFi;
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
    memcpy(_TFi, genome->get_TFi(), sizeof(size_t)*_nb_TF);
  }
  delete[] _Pi;
  _Pi = NULL;
  if (_nb_P > 0)
  {
    _Pi = new size_t[_nb_P];
    memcpy(_Pi, genome->get_Pi(), sizeof(size_t)*_nb_P);
  }
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Do point mutations
 * \details  Browse the genome and apply Gene::mutate() on each gene
 * \param    const double* mutation_rates
 * \return   \e void
 */
void Genome::do_point_mutations( const double* mutation_rates )
{
  _nb_NC            = 0;
  _nb_E             = 0;
  _nb_TF            = 0;
  _nb_BS            = 0;
  _nb_P             = 0;
  _nb_inner_enzymes = 0;
  _nb_inflow_pumps  = 0;
  _nb_outflow_pumps = 0;
  genetic_unit* str = _genetic_sequence->x;
  for (size_t pos = 0; pos < _genetic_sequence->size; pos++)
  {
    /*-------------------------*/
    /* 1.1) mutate the gene    */
    /*-------------------------*/
    do_genetic_unit_mutation(mutation_rates, pos);
    /*-------------------------*/
    /* 1.2) compute statistics */
    /*-------------------------*/
    if (str[pos].type == NON_CODING)
    {
      _nb_NC++;
    }
    else if (str[pos].type == ENZYME)
    {
      _nb_E++;
      if (str[pos].kcat > 0.0 && str[pos].s == str[pos].p)
      {
        _nb_inflow_pumps++;
      }
      else if (str[pos].kcat < 0.0 && str[pos].s == str[pos].p)
      {
        _nb_outflow_pumps++;
      }
      else if (str[pos].kcat != 0.0 && str[pos].s != str[pos].p)
      {
        _nb_inner_enzymes++;
      }
    }
    else if (str[pos].type == TRANSCRIPTION_FACTOR)
    {
      _nb_TF++;
    }
    else if (str[pos].type == BINDING_SITE)
    {
      _nb_BS++;
    }
    else if (str[pos].type == PROMOTER)
    {
      _nb_P++;
    }
  }
}

/*
 * \brief    Do rearrangements (duplications + deletions + translocations + inversions)
 * \details  Each breakpoint can mutate one of the two adjacent genes
 * \param    const double* mutation_rates
 * \return   \e void
 */
void Genome::do_rearrangements( const double* mutation_rates )
{
  if (_genetic_sequence->size == 0)
  {
    return;
  }
  size_t ndup = _prng->binomial(_genetic_sequence->size, mutation_rates[DUPLICATION_RATE]);
  size_t ndel = _prng->binomial(_genetic_sequence->size, mutation_rates[DELETION_RATE]);
  size_t ntra = _prng->binomial(_genetic_sequence->size, mutation_rates[TRANSLOCATION_RATE]);
  size_t ninv = _prng->binomial(_genetic_sequence->size, mutation_rates[INVERSION_RATE]);
  size_t nrea = ndup + ndel + ntra + ninv;
  for (size_t i = nrea; i >= 1; i--)
  {
    size_t draw = (size_t)_prng->uniform(0, (int)i);
    /*------------------*/
    /* Do duplication   */
    /*------------------*/
    if ( draw < ndup )
    {
      do_duplication(mutation_rates);
      ndup--;
    }
    /*------------------*/
    /* Do deletion      */
    /*------------------*/
    else if ( draw < ndup + ndel )
    {
      do_deletion(mutation_rates);
      ndel--;
    }
    /*------------------*/
    /* Do translocation */
    /*------------------*/
    else if ( draw < ndup + ndel + ntra )
    {
      do_translocation(mutation_rates);
      ntra--;
    }
    /*------------------*/
    /* Do inversion     */
    /*------------------*/
    else
    {
      do_inversion(mutation_rates);
      ninv--;
    }
    check_genome_size();
  }
}

/**
 * \brief    Do a duplication
 * \details  Draw uniformly two breakpoints in the genome (start and end), and copy the sequence at a third random breakpoint.
 * \param    const double* mutation_rates
 * \return   \e void
 */
void Genome::do_duplication( const double* mutation_rates )
{
  if (_genetic_sequence->size == 0)
  {
    return;
  }
  
  /*********************************************/
  /* Draw duplication parameters               */
  /*********************************************/
  
  size_t N      = _genetic_sequence->size;
  size_t start  = (size_t)_prng->uniform(0, (int)N-1);
  size_t end    = (size_t)_prng->uniform(0, (int)N-1);
  size_t insert = (size_t)_prng->uniform(0, (int)N-1);
  size_t size   = 0;
  
  /******************************************/
  /* Count the number of units per type     */
  /******************************************/
  
  size_t nb_NC = 0;
  size_t nb_E  = 0;
  size_t nb_TF = 0;
  size_t nb_BS = 0;
  size_t nb_P  = 0;
  size_t index = start;
  while (index != end)
  {
    if (_genetic_sequence->x[index].type == NON_CODING)
    {
      nb_NC++;
    }
    else if (_genetic_sequence->x[index].type == ENZYME)
    {
      nb_E++;
    }
    else if (_genetic_sequence->x[index].type == TRANSCRIPTION_FACTOR)
    {
      nb_TF++;
    }
    else if (_genetic_sequence->x[index].type == BINDING_SITE)
    {
      nb_BS++;
    }
    else if (_genetic_sequence->x[index].type == PROMOTER)
    {
      nb_P++;
    }
    index = (index+1+_genetic_sequence->size)%_genetic_sequence->size;
  }
  if (_genetic_sequence->x[index].type == NON_CODING)
  {
    nb_NC++;
  }
  else if (_genetic_sequence->x[index].type == ENZYME)
  {
    nb_E++;
  }
  else if (_genetic_sequence->x[index].type == TRANSCRIPTION_FACTOR)
  {
    nb_TF++;
  }
  else if (_genetic_sequence->x[index].type == BINDING_SITE)
  {
    nb_BS++;
  }
  else if (_genetic_sequence->x[index].type == PROMOTER)
  {
    nb_P++;
  }
  
  /*********************************************/
  /* Compute the duplication when start <= end */
  /*********************************************/
  
  if (start <= end)
  {
    size_t duplicata_size = end-start+1;
    size_t shift_size     = N-insert-1;
    size_t new_size       = N+duplicata_size;
    if (_genetic_sequence->size+new_size > _genetic_sequence->buffer_size)
    {
      increase_buffer_size(new_size);
    }
    genetic_unit* duplicata = NULL;
    genetic_unit* shift     = NULL;
    if (duplicata_size > 0)
    {
      duplicata = new genetic_unit[duplicata_size];
      for (size_t i = 0; i < duplicata_size; i++)
      {
        assert(start+i < N);
        duplicata[i] = _genetic_sequence->x[start+i];
      }
    }
    if (shift_size > 0)
    {
      shift = new genetic_unit[shift_size];
      for (size_t i = 0; i < shift_size; i++)
      {
        assert(insert+1+i < N);
        shift[i] = _genetic_sequence->x[insert+1+i];
      }
    }
    for (size_t i = 0; i < duplicata_size; i++)
    {
      assert(insert+1+i < new_size);
      _genetic_sequence->x[insert+1+i] = duplicata[i];
    }
    for (size_t i = 0; i < shift_size; i++)
    {
      assert(insert+1+duplicata_size+i < new_size);
      _genetic_sequence->x[insert+1+duplicata_size+i] = shift[i];
    }
    _genetic_sequence->size = new_size;
    size_t pos1 = (insert+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos2 = (insert+1+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos1, pos2);
    size_t pos3 = (insert+duplicata_size+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos4 = (insert+duplicata_size+1+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos3, pos4);
    delete[] duplicata;
    duplicata = NULL;
    delete[] shift;
    shift = NULL;
    size = duplicata_size;
  }
  
  /*********************************************/
  /* Compute the duplication when start > end  */
  /*********************************************/
  
  if (start > end)
  {
    size_t duplicata_size = N-start+end+1;
    size_t shift_size     = N-insert-1;
    size_t new_size       = N+duplicata_size;
    if (_genetic_sequence->size+new_size > _genetic_sequence->buffer_size)
    {
      increase_buffer_size(new_size);
    }
    genetic_unit* duplicata = NULL;
    genetic_unit* shift     = NULL;
    if (duplicata_size > 0)
    {
      duplicata = new genetic_unit[duplicata_size];
      for (size_t i = 0; i < N-start; i++)
      {
        assert(start+i < N);
        duplicata[i] = _genetic_sequence->x[start+i];
      }
      for (size_t i = 0; i < end+1; i++)
      {
        assert(i < N);
        duplicata[N-start+i] = _genetic_sequence->x[i];
      }
    }
    if (shift_size > 0)
    {
      shift = new genetic_unit[shift_size];
      for (size_t i = 0; i < shift_size; i++)
      {
        assert(insert+1+i < N);
        shift[i] = _genetic_sequence->x[insert+1+i];
      }
    }
    for (size_t i = 0; i < duplicata_size; i++)
    {
      assert(insert+1+i < new_size);
      _genetic_sequence->x[insert+1+i] = duplicata[i];
    }
    for (size_t i = 0; i < shift_size; i++)
    {
      assert(insert+1+duplicata_size+i < new_size);
      _genetic_sequence->x[insert+1+duplicata_size+i] = shift[i];
    }
    _genetic_sequence->size = new_size;
    size_t pos1 = (insert+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos2 = (insert+1+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos1, pos2);
    size_t pos3 = (insert+duplicata_size+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos4 = (insert+duplicata_size+1+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos3, pos4);
    delete[] duplicata;
    duplicata = NULL;
    delete[] shift;
    shift = NULL;
    size = duplicata_size;
  }
  
  /*********************************************/
  /* Save the mutation event                   */
  /*********************************************/
  
  MutationEvent* new_event = new MutationEvent(DUPLICATION, start, end, insert, size, nb_NC, nb_E, nb_TF, nb_BS, nb_P);
  _replication_report->add_mutation_event(new_event);
}

/**
 * \brief    Do a deletion
 * \details  Draw uniformly two breakpoints in the genome (start and end), and delete the sequence.
 * \param    const double* mutation_rates
 * \return   \e void
 */
void Genome::do_deletion( const double* mutation_rates )
{
  if (_genetic_sequence->size == 0)
  {
    return;
  }
  
  /******************************************/
  /* Draw deletion parameters               */
  /******************************************/
  
  size_t N     = _genetic_sequence->size;
  size_t start = (size_t)_prng->uniform(0, (int)N-1);
  size_t end   = (size_t)_prng->uniform(0, (int)N-1);
  size_t size  = 0;
  
  /******************************************/
  /* Count the number of units per type     */
  /******************************************/
  
  size_t nb_NC = 0;
  size_t nb_E  = 0;
  size_t nb_TF = 0;
  size_t nb_BS = 0;
  size_t nb_P  = 0;
  size_t index = start;
  while (index != end)
  {
    if (_genetic_sequence->x[index].type == NON_CODING)
    {
      nb_NC++;
    }
    else if (_genetic_sequence->x[index].type == ENZYME)
    {
      nb_E++;
    }
    else if (_genetic_sequence->x[index].type == TRANSCRIPTION_FACTOR)
    {
      nb_TF++;
    }
    else if (_genetic_sequence->x[index].type == BINDING_SITE)
    {
      nb_BS++;
    }
    else if (_genetic_sequence->x[index].type == PROMOTER)
    {
      nb_P++;
    }
    index = (index+1+_genetic_sequence->size)%_genetic_sequence->size;
  }
  if (_genetic_sequence->x[index].type == NON_CODING)
  {
    nb_NC++;
  }
  else if (_genetic_sequence->x[index].type == ENZYME)
  {
    nb_E++;
  }
  else if (_genetic_sequence->x[index].type == TRANSCRIPTION_FACTOR)
  {
    nb_TF++;
  }
  else if (_genetic_sequence->x[index].type == BINDING_SITE)
  {
    nb_BS++;
  }
  else if (_genetic_sequence->x[index].type == PROMOTER)
  {
    nb_P++;
  }
  
  /******************************************/
  /* Compute the deletion when start <= end */
  /******************************************/
  
  if (start <= end)
  {
    size_t pos1 = (start-1+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos2 = (end+1+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos1, pos2);
    size_t shift_size   = N-end-1;
    size_t new_size     = N-end+start-1;
    genetic_unit* shift = NULL;
    if (shift_size > 0)
    {
      shift = new genetic_unit[shift_size];
      for (size_t i = 0; i < shift_size; i++)
      {
        assert(end+1+i < N);
        shift[i] = _genetic_sequence->x[end+1+i];
      }
    }
    for (size_t i = 0; i < shift_size; i++)
    {
      assert(start+i < new_size);
      _genetic_sequence->x[start+i] = shift[i];
    }
    _genetic_sequence->size = new_size;
    if (_genetic_sequence->size <= _genetic_sequence->buffer_size-3*GENOME_BUFFER/2)
    {
      decrease_buffer_size();
    }
    delete[] shift;
    shift = NULL;
    size = N-_genetic_sequence->size;
  }
  
  /******************************************/
  /* Compute the deletion when start > end  */
  /******************************************/
  
  if (start > end)
  {
    size_t pos1 = (end+1+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos2 = (start-1+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos1, pos2);
    size_t shift_size   = start-end-1;
    size_t new_size     = start-end-1;
    genetic_unit* shift = NULL;
    if (shift_size > 0)
    {
      shift = new genetic_unit[shift_size];
      for (size_t i = 0; i < shift_size; i++)
      {
        assert(end+1+i < N);
        shift[i] = _genetic_sequence->x[end+1+i];
      }
    }
    for (size_t i = 0; i < shift_size; i++)
    {
      assert(i < new_size);
      _genetic_sequence->x[i] = shift[i];
    }
    _genetic_sequence->size = new_size;
    if (_genetic_sequence->size <= _genetic_sequence->buffer_size-3*GENOME_BUFFER/2)
    {
      decrease_buffer_size();
    }
    delete[] shift;
    shift = NULL;
    size = N-_genetic_sequence->size;
  }
  
  /******************************************/
  /* Save the mutation event                */
  /******************************************/
  
  MutationEvent* new_event = new MutationEvent(DELETION, start, end, 0, size, nb_NC, nb_E, nb_TF, nb_BS, nb_P);
  _replication_report->add_mutation_event(new_event);
}

/**
 * \brief    Do a translocation
 * \details  Draw uniformly two breakpoints in the genome (start and end), and translocate the sequence at a third random breakpoint.
 * \param    const double* mutation_rates
 * \return   \e void
 */
void Genome::do_translocation( const double* mutation_rates )
{
  if (_genetic_sequence->size == 0)
  {
    return;
  }
  
  /***********************************************/
  /* Draw first translocation parameters         */
  /***********************************************/
  
  size_t N      = _genetic_sequence->size;
  size_t start  = (size_t)_prng->uniform(0, (int)N-1);
  size_t end    = (size_t)_prng->uniform(0, (int)N-1);
  size_t insert = 0;
  size_t size   = 0;
  
  /***********************************************/
  /* Count the number of units per type          */
  /***********************************************/
  
  size_t nb_NC = 0;
  size_t nb_E  = 0;
  size_t nb_TF = 0;
  size_t nb_BS = 0;
  size_t nb_P  = 0;
  size_t index = start;
  while (index != end)
  {
    if (_genetic_sequence->x[index].type == NON_CODING)
    {
      nb_NC++;
    }
    else if (_genetic_sequence->x[index].type == ENZYME)
    {
      nb_E++;
    }
    else if (_genetic_sequence->x[index].type == TRANSCRIPTION_FACTOR)
    {
      nb_TF++;
    }
    else if (_genetic_sequence->x[index].type == BINDING_SITE)
    {
      nb_BS++;
    }
    else if (_genetic_sequence->x[index].type == PROMOTER)
    {
      nb_P++;
    }
    index = (index+1+_genetic_sequence->size)%_genetic_sequence->size;
  }
  if (_genetic_sequence->x[index].type == NON_CODING)
  {
    nb_NC++;
  }
  else if (_genetic_sequence->x[index].type == ENZYME)
  {
    nb_E++;
  }
  else if (_genetic_sequence->x[index].type == TRANSCRIPTION_FACTOR)
  {
    nb_TF++;
  }
  else if (_genetic_sequence->x[index].type == BINDING_SITE)
  {
    nb_BS++;
  }
  else if (_genetic_sequence->x[index].type == PROMOTER)
  {
    nb_P++;
  }
  
  /***********************************************/
  /* Compute the translocation when start <= end */
  /***********************************************/
  
  if (start <= end)
  {
    /*----------------------------------------------*/
    /* 1) Step 1 : copy the duplicata and delete it */
    /*----------------------------------------------*/
    size_t duplicata_size   = end-start+1;
    size_t shift_size       = N-end-1;
    size_t new_size         = N-end+start-1;
    genetic_unit* duplicata = NULL;
    genetic_unit* shift     = NULL;
    if (duplicata_size > 0)
    {
      duplicata = new genetic_unit[duplicata_size];
      for (size_t i = 0; i < duplicata_size; i++)
      {
        assert(start+i < N);
        duplicata[i] = _genetic_sequence->x[start+i];
      }
    }
    if (shift_size > 0)
    {
      shift = new genetic_unit[shift_size];
      for (size_t i = 0; i < shift_size; i++)
      {
        assert(end+1+i < N);
        shift[i] = _genetic_sequence->x[end+1+i];
      }
    }
    for (size_t i = 0; i < shift_size; i++)
    {
      assert(start+i < new_size);
      _genetic_sequence->x[start+i] = shift[i];
    }
    _genetic_sequence->size = new_size;
    delete[] shift;
    shift = NULL;
    
    /*----------------------------------------------*/
    /* 2) Step 2 : insert the copied sequence       */
    /*----------------------------------------------*/
    N        = _genetic_sequence->size;
    new_size = _genetic_sequence->size+duplicata_size;
    if (N > 0)
    {
      insert     = (size_t)_prng->uniform(0, (int)N-1);
      shift_size = N-insert-1;
      if (shift_size > 0)
      {
        shift = new genetic_unit[shift_size];
        for (size_t i = 0; i < shift_size; i++)
        {
          assert(insert+1+i < N);
          shift[i] = _genetic_sequence->x[insert+1+i];
        }
      }
      for (size_t i = 0; i < duplicata_size; i++)
      {
        assert(insert+1+i < new_size);
        _genetic_sequence->x[insert+1+i] = duplicata[i];
      }
      for (size_t i = 0; i < shift_size; i++)
      {
        assert(insert+1+duplicata_size+i < new_size);
        _genetic_sequence->x[insert+1+duplicata_size+i] = shift[i];
      }
    }
    else
    {
      for (size_t i = 0; i < duplicata_size; i++)
      {
        assert(i < new_size);
        _genetic_sequence->x[i] = duplicata[i];
      }
    }
    _genetic_sequence->size = new_size;
    size_t pos1 = (start-1+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos2 = (end+1+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos1, pos2);
    size_t pos3 = (insert+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos4 = (insert+1+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos3, pos4);
    size_t pos5 = (insert+duplicata_size+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos6 = (insert+duplicata_size+1+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos5, pos6);
    delete[] duplicata;
    duplicata = NULL;
    delete[] shift;
    shift = NULL;
    size = duplicata_size;
  }
  
  /***********************************************/
  /* Compute the translocation when start > end  */
  /***********************************************/
  
  if (start > end)
  {
    /*----------------------------------------------*/
    /* 1) Step 1 : copy the duplicata and delete it */
    /*----------------------------------------------*/
    size_t duplicata_size   = N-start+end+1;
    size_t shift_size       = start-end-1;
    size_t new_size         = start-end-1;
    genetic_unit* duplicata = NULL;
    genetic_unit* shift     = NULL;
    if (duplicata_size > 0)
    {
      duplicata = new genetic_unit[duplicata_size];
      for (size_t i = 0; i < N-start; i++)
      {
        assert(start+i < N);
        duplicata[i] = _genetic_sequence->x[start+i];
      }
      for (size_t i = 0; i < end+1; i++)
      {
        assert(i < N);
        duplicata[N-start+i] = _genetic_sequence->x[i];
      }
    }
    if (shift_size > 0)
    {
      shift = new genetic_unit[shift_size];
      for (size_t i = 0; i < shift_size; i++)
      {
        assert(end+1+i < N);
        shift[i] = _genetic_sequence->x[end+1+i];
      }
    }
    for (size_t i = 0; i < shift_size; i++)
    {
      assert(i < new_size);
      _genetic_sequence->x[i] = shift[i];
    }
    _genetic_sequence->size = new_size;
    delete[] shift;
    shift = NULL;
    
    /*----------------------------------------------*/
    /* 2) Step 2 : insert the copied sequence       */
    /*----------------------------------------------*/
    N        = _genetic_sequence->size;
    new_size = _genetic_sequence->size+duplicata_size;
    if (N > 0)
    {
      insert     = (size_t)_prng->uniform(0, (int)N-1);
      shift_size = N-insert-1;
      if (shift_size > 0)
      {
        shift = new genetic_unit[shift_size];
        for (size_t i = 0; i < shift_size; i++)
        {
          assert(insert+1+i < N);
          shift[i] = _genetic_sequence->x[insert+1+i];
        }
      }
      for (size_t i = 0; i < duplicata_size; i++)
      {
        assert(insert+1+i < new_size);
        _genetic_sequence->x[insert+1+i] = duplicata[i];
      }
      for (size_t i = 0; i < shift_size; i++)
      {
        assert(insert+1+duplicata_size+i < new_size);
        _genetic_sequence->x[insert+1+duplicata_size+i] = shift[i];
      }
    }
    else
    {
      new_size = _genetic_sequence->size+duplicata_size;
      for (size_t i = 0; i < duplicata_size; i++)
      {
        assert(i < new_size);
        _genetic_sequence->x[i] = duplicata[i];
      }
    }
    _genetic_sequence->size = new_size;
    size_t pos1 = (start-1+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos2 = (end+1+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos1, pos2);
    size_t pos3 = (insert+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos4 = (insert+1+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos3, pos4);
    size_t pos5 = (insert+duplicata_size+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos6 = (insert+duplicata_size+1+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos5, pos6);
    delete[] duplicata;
    duplicata = NULL;
    delete[] shift;
    shift = NULL;
    size = duplicata_size;
  }
  
  /***********************************************/
  /* Save the mutation event                     */
  /***********************************************/
  
  MutationEvent* new_event = new MutationEvent(TRANSLOCATION, start, end, insert, size, nb_NC, nb_E, nb_TF, nb_BS, nb_P);
  _replication_report->add_mutation_event(new_event);
}

/**
 * \brief    Do an inversion
 * \details  Draw uniformly two breakpoints in the genome (start end end), and revert the sequence.
 * \param    const double* mutation_rates
 * \return   \e void
 */
void Genome::do_inversion( const double* mutation_rates )
{
  if (_genetic_sequence->size == 0)
  {
    return;
  }
  
  /*******************************************/
  /* Draw inversion parameters               */
  /*******************************************/
  
  size_t start         = (size_t)_prng->uniform(0, (int)_genetic_sequence->size-1);
  size_t end           = (size_t)_prng->uniform(0, (int)_genetic_sequence->size-1);
  size_t current_start = start;
  size_t current_end   = end;
  size_t size          = 0;
  
  /*******************************************/
  /* Count the number of units per type      */
  /*******************************************/
  
  size_t nb_NC = 0;
  size_t nb_E  = 0;
  size_t nb_TF = 0;
  size_t nb_BS = 0;
  size_t nb_P  = 0;
  size_t index = start;
  while (index != end)
  {
    if (_genetic_sequence->x[index].type == NON_CODING)
    {
      nb_NC++;
    }
    else if (_genetic_sequence->x[index].type == ENZYME)
    {
      nb_E++;
    }
    else if (_genetic_sequence->x[index].type == TRANSCRIPTION_FACTOR)
    {
      nb_TF++;
    }
    else if (_genetic_sequence->x[index].type == BINDING_SITE)
    {
      nb_BS++;
    }
    else if (_genetic_sequence->x[index].type == PROMOTER)
    {
      nb_P++;
    }
    index = (index+1+_genetic_sequence->size)%_genetic_sequence->size;
  }
  if (_genetic_sequence->x[index].type == NON_CODING)
  {
    nb_NC++;
  }
  else if (_genetic_sequence->x[index].type == ENZYME)
  {
    nb_E++;
  }
  else if (_genetic_sequence->x[index].type == TRANSCRIPTION_FACTOR)
  {
    nb_TF++;
  }
  else if (_genetic_sequence->x[index].type == BINDING_SITE)
  {
    nb_BS++;
  }
  else if (_genetic_sequence->x[index].type == PROMOTER)
  {
    nb_P++;
  }
  
  /*******************************************/
  /* Compute the inversion when start <= end */
  /*******************************************/
  
  if (start <= end)
  {
    size = end-start+1;
    while (current_start < current_end)
    {
      genetic_unit tmp                    = _genetic_sequence->x[current_start];
      _genetic_sequence->x[current_start] = _genetic_sequence->x[current_end];
      _genetic_sequence->x[current_end]   = tmp;
      current_start = (current_start+1+_genetic_sequence->size)%_genetic_sequence->size;
      current_end   = (current_end-1+_genetic_sequence->size)%_genetic_sequence->size;
    }
    size_t pos1 = (start-1+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos2 = (start+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos1, pos2);
    size_t pos3 = (end+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos4 = (end+1+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos3, pos4);
  }
  
  /*******************************************/
  /* Compute the inversion when start > end  */
  /*******************************************/
  
  else if (start > end)
  {
    size = _genetic_sequence->size-start+end+1;
    while (current_start > current_end)
    {
      genetic_unit tmp                    = _genetic_sequence->x[current_start];
      _genetic_sequence->x[current_start] = _genetic_sequence->x[current_end];
      _genetic_sequence->x[current_end]   = tmp;
      current_start = (current_start+1+_genetic_sequence->size)%_genetic_sequence->size;
      current_end   = (current_end-1+_genetic_sequence->size)%_genetic_sequence->size;
    }
    size_t pos1 = (end-1+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos2 = (end+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos1, pos2);
    size_t pos3 = (start+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos4 = (start+1+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos3, pos4);
  }
  
  /*******************************************/
  /* Save the mutation event                 */
  /*******************************************/
  
  MutationEvent* new_event = new MutationEvent(INVERSION, start, end, 0, size, nb_NC, nb_E, nb_TF, nb_BS, nb_P);
  _replication_report->add_mutation_event(new_event);
}

/**
 * \brief    Perform horizontal gene transfer (HGT)
 * \details  The number of new random sequences introduced in the genome depends on a Poisson law of parameter 'HGT_RATE'. Pasted genetic units at insertion point undergo BLX-a crossover
 * \param    const double* mutation_rates
 * \return   \e void
 */
void Genome::do_hgt( const double* mutation_rates )
{
  /*--------------------------------------*/
  /* 1) Draw the number of HGT to perform */
  /*--------------------------------------*/
  int nb_draws = _prng->poisson(_parameters->get_hgt_rate());
  
  /*--------------------------------------*/
  /* 2) For each new HGT                  */
  /*--------------------------------------*/
  for (int i = 0; i < nb_draws; i++)
  {
    /* 2.1) Draw the HGT sequence */
    size_t insert       = 0;
    genetic_unit* shift = NULL;
    size_t nb_NC        = 0;
    size_t nb_E         = 0;
    size_t nb_TF        = 0;
    size_t nb_BS        = 0;
    size_t nb_P         = 0;
    size_t size         = (size_t)_prng->uniform(HGT_MIN_SIZE, HGT_MAX_SIZE);
    genetic_unit* hgt   = new genetic_unit[size];
    for (size_t j = 0; j < size; j++)
    {
      draw_random_genetic_unit(hgt[j]);
      if (hgt[j].type == NON_CODING)
      {
        nb_NC++;
      }
      else if (hgt[j].type == ENZYME)
      {
        nb_E++;
      }
      else if (hgt[j].type == TRANSCRIPTION_FACTOR)
      {
        nb_TF++;
      }
      else if (hgt[j].type == BINDING_SITE)
      {
        nb_BS++;
      }
      else if (hgt[j].type == PROMOTER)
      {
        nb_P++;
      }
    }
    
    /* 2.2) Insert the HGT sequence in the genome */
    if (_genetic_sequence->size == 0)
    {
      insert = 0;
    }
    else
    {
      insert = _prng->uniform(0, (int)_genetic_sequence->size-1);
    }
    
    shift = new genetic_unit[_genetic_sequence->size-insert];
    memcpy(shift, &_genetic_sequence->x[insert], sizeof(genetic_unit)*(_genetic_sequence->size-insert));
    memcpy(&_genetic_sequence->x[insert], hgt, sizeof(genetic_unit)*size);
    memcpy(&_genetic_sequence->x[insert+size], shift, sizeof(genetic_unit)*(_genetic_sequence->size-insert));
    _genetic_sequence->size += size;
    
    /* 2.3) Delete sequences */
    delete[] hgt;
    hgt = NULL;
    delete[] shift;
    shift = NULL;
    
    /* 2.4) perform crossover */
    size_t pos1 = (insert-1+_genetic_sequence->size)%_genetic_sequence->size;
    size_t pos2 = (insert+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos1, pos2);
    pos1 = (insert+size-1+_genetic_sequence->size)%_genetic_sequence->size;
    pos2 = (insert+size+_genetic_sequence->size)%_genetic_sequence->size;
    do_permutation(mutation_rates, pos1, pos2);
    
    /* 2.5) Save the mutation event */
    MutationEvent* new_event = new MutationEvent(HGT, insert, size, nb_NC, nb_E, nb_TF, nb_BS, nb_P);
    _replication_report->add_mutation_event(new_event);
  }
}

/**
 * \brief    Mutate a genetic unit
 * \details  Apply point mutations on a genetic unit
 * \param    const double* mutation_rates
 * \param    size_t pos
 * \return   \e bool
 */
void Genome::do_genetic_unit_mutation( const double* mutation_rates, size_t pos )
{
  MutationVector* dX           = new MutationVector();
  genetic_unit& unit           = _genetic_sequence->x[pos];
  dX->get_dX()->type           = unit.type;
  dX->get_dX()->free_activity  = unit.free_activity;
  dX->get_dX()->bound_activity = unit.bound_activity;
  dX->get_dX()->binding_window = unit.binding_window;
  bool mutate                  = dX->draw(_prng, mutation_rates);
  
  /*----------------------------------------*/
  /* if at least one point mutation occured */
  /*----------------------------------------*/
  if (mutate)
  {
    genetic_unit* dx = dX->get_dX();
    
    /*------------------------------------------------------------------ Global attributes */
    
    if (dx->type == NON_CODING ||
        dx->type == E_TO_NC_TRANSITION ||
        dx->type == TF_TO_NC_TRANSITION ||
        dx->type == BS_TO_NC_TRANSITION ||
        dx->type == P_TO_NC_TRANSITION)
    {
      unit.type = NON_CODING;
    }
    else if (dx->type == ENZYME ||
             dx->type == NC_TO_E_TRANSITION ||
             dx->type == TF_TO_E_TRANSITION ||
             dx->type == BS_TO_E_TRANSITION ||
             dx->type == P_TO_E_TRANSITION)
    {
      unit.type = ENZYME;
    }
    else if (dx->type == TRANSCRIPTION_FACTOR ||
             dx->type == NC_TO_TF_TRANSITION ||
             dx->type == E_TO_TF_TRANSITION ||
             dx->type == BS_TO_TF_TRANSITION ||
             dx->type == P_TO_TF_TRANSITION)
    {
      unit.type = TRANSCRIPTION_FACTOR;
    }
    else if (dx->type == BINDING_SITE ||
             dx->type == NC_TO_BS_TRANSITION ||
             dx->type == E_TO_BS_TRANSITION ||
             dx->type == TF_TO_BS_TRANSITION ||
             dx->type == P_TO_BS_TRANSITION)
    {
      unit.type = BINDING_SITE;
    }
    else if (dx->type == PROMOTER ||
             dx->type == NC_TO_P_TRANSITION ||
             dx->type == E_TO_P_TRANSITION ||
             dx->type == TF_TO_P_TRANSITION ||
             dx->type == BS_TO_P_TRANSITION)
    {
      unit.type = PROMOTER;
    }
    
    /*------------------------------------------------------------------ Enzyme type (E) attributes */
    
    unit.s += dx->s;
    unit.s  = (unit.s > 0 ? unit.s : 1);
    unit.p += dx->p;
    unit.p  = (unit.p > 0 ? unit.p : 1);
    
    double kcat = fabs(unit.kcat);
    kcat = pow(10.0, log10(kcat)+dx->kcat);
    if (unit.kcat < 0.0)
    {
      if (kcat < pow(10.0, KCAT_MIN_LOG))
      {
        unit.kcat = pow(10.0, KCAT_MIN_LOG);
      }
      else if (kcat > pow(10.0, KCAT_MAX_LOG))
      {
        unit.kcat = -pow(10.0, KCAT_MAX_LOG);
      }
      else
      {
        unit.kcat = -kcat;
      }
    }
    else if (unit.kcat > 0.0)
    {
      if (kcat < pow(10.0, KCAT_MIN_LOG))
      {
        unit.kcat = -pow(10.0, KCAT_MIN_LOG);
      }
      else if (kcat > pow(10.0, KCAT_MAX_LOG))
      {
        unit.kcat = pow(10.0, KCAT_MAX_LOG);
      }
      else
      {
        unit.kcat = kcat;
      }
    }
    
    unit.kcat_km_ratio = pow(10.0, log10(unit.kcat_km_ratio)+dx->kcat_km_ratio);
    if (unit.kcat_km_ratio < pow(10.0, KCAT_KM_RATIO_MIN_LOG))
    {
      unit.kcat_km_ratio = pow(10.0, KCAT_KM_RATIO_MIN_LOG);
    }
    else if (unit.kcat_km_ratio > pow(10.0, KCAT_KM_RATIO_MAX_LOG))
    {
      unit.kcat_km_ratio = pow(10.0, KCAT_KM_RATIO_MAX_LOG);
    }
    
    /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
    
    unit.BS_tag         += dx->BS_tag;
    unit.coE_tag        += dx->coE_tag;
    unit.coE_tag         = (unit.coE_tag > 0 ? unit.coE_tag : 1);
    unit.free_activity   = dx->free_activity;
    unit.bound_activity  = dx->bound_activity;
    unit.binding_window  = dx->binding_window;
    
    /*------------------------------------------------------------------ Binding site type (BS) attributes */
    
    unit.TF_tag += dx->TF_tag;
    
    /*------------------------------------------------------------------ Promoter type (P) attributes */
    
    unit.basal_expression_level += dx->basal_expression_level;
    unit.basal_expression_level  = (unit.basal_expression_level > 0.0 ? unit.basal_expression_level : 0.0);
    unit.basal_expression_level  = (unit.basal_expression_level < 1.0 ? unit.basal_expression_level : 1.0);
    
    /*------------------------------------------------------------------ Save point mutation event */
    
    MutationEvent* new_event = new MutationEvent(POINT_MUTATION, pos, dX);
    _replication_report->add_mutation_event(new_event);
  }
  else
  {
    delete dX;
    dX = NULL;
  }
}

/**
 * \brief    Mutate a genetic unit at breakpoints
 * \details  Apply point mutations on a genetic unit
 * \param    const double* mutation_rates
 * \param    size_t pos
 * \return   \e void
 */
void Genome::do_genetic_unit_mutation_at_breakpoints( const double* mutation_rates, size_t pos )
{
  MutationVector* dX           = new MutationVector();
  genetic_unit& p              = _genetic_sequence->x[pos];
  dX->get_dX()->type           = p.type;
  dX->get_dX()->free_activity  = p.free_activity;
  dX->get_dX()->bound_activity = p.bound_activity;
  dX->get_dX()->binding_window = p.binding_window;
  bool mutate                  = dX->breakpoint_draw(_prng, mutation_rates);
  
  /*----------------------------------------*/
  /* if at least one point mutation occured */
  /*----------------------------------------*/
  if (mutate)
  {
    genetic_unit* dx = dX->get_dX();
    
    /*------------------------------------------------------------------ Global attributes */
    
    if (dx->type == NON_CODING ||
        dx->type == E_TO_NC_TRANSITION ||
        dx->type == TF_TO_NC_TRANSITION ||
        dx->type == BS_TO_NC_TRANSITION ||
        dx->type == P_TO_NC_TRANSITION)
    {
      p.type = NON_CODING;
    }
    else if (dx->type == ENZYME ||
             dx->type == NC_TO_E_TRANSITION ||
             dx->type == TF_TO_E_TRANSITION ||
             dx->type == BS_TO_E_TRANSITION ||
             dx->type == P_TO_E_TRANSITION)
    {
      p.type = ENZYME;
    }
    else if (dx->type == TRANSCRIPTION_FACTOR ||
             dx->type == NC_TO_TF_TRANSITION ||
             dx->type == E_TO_TF_TRANSITION ||
             dx->type == BS_TO_TF_TRANSITION ||
             dx->type == P_TO_TF_TRANSITION)
    {
      p.type = TRANSCRIPTION_FACTOR;
    }
    else if (dx->type == BINDING_SITE ||
             dx->type == NC_TO_BS_TRANSITION ||
             dx->type == E_TO_BS_TRANSITION ||
             dx->type == TF_TO_BS_TRANSITION ||
             dx->type == P_TO_BS_TRANSITION)
    {
      p.type = BINDING_SITE;
    }
    else if (dx->type == PROMOTER ||
             dx->type == NC_TO_P_TRANSITION ||
             dx->type == E_TO_P_TRANSITION ||
             dx->type == TF_TO_P_TRANSITION ||
             dx->type == BS_TO_P_TRANSITION)
    {
      p.type = PROMOTER;
    }
    
    /*------------------------------------------------------------------ Enzyme type (E) attributes */
    
    p.s += dx->s;
    p.s  = (p.s > 0 ? p.s : 1);
    p.p += dx->p;
    p.p  = (p.p > 0 ? p.p : 1);
    
    double kcat = fabs(p.kcat);
    kcat = pow(10.0, log10(kcat)+dx->kcat);
    if (p.kcat < 0.0)
    {
      if (kcat < pow(10.0, KCAT_MIN_LOG))
      {
        p.kcat = pow(10.0, KCAT_MIN_LOG);
      }
      else if (kcat > pow(10.0, KCAT_MAX_LOG))
      {
        p.kcat = -pow(10.0, KCAT_MAX_LOG);
      }
      else
      {
        p.kcat = -kcat;
      }
    }
    else if (p.kcat > 0.0)
    {
      if (kcat < pow(10.0, KCAT_MIN_LOG))
      {
        p.kcat = -pow(10.0, KCAT_MIN_LOG);
      }
      else if (kcat > pow(10.0, KCAT_MAX_LOG))
      {
        p.kcat = pow(10.0, KCAT_MAX_LOG);
      }
      else
      {
        p.kcat = kcat;
      }
    }
    
    p.kcat_km_ratio = pow(10.0, log10(p.kcat_km_ratio)+dx->kcat_km_ratio);
    if (p.kcat_km_ratio < pow(10.0, KCAT_KM_RATIO_MIN_LOG))
    {
      p.kcat_km_ratio = pow(10.0, KCAT_KM_RATIO_MIN_LOG);
    }
    else if (p.kcat_km_ratio > pow(10.0, KCAT_KM_RATIO_MAX_LOG))
    {
      p.kcat_km_ratio = pow(10.0, KCAT_KM_RATIO_MAX_LOG);
    }
    
    /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
    
    p.BS_tag         += dx->BS_tag;
    p.coE_tag        += dx->coE_tag;
    p.coE_tag         = (p.coE_tag > 0 ? p.coE_tag : 1);
    p.free_activity   = dx->free_activity;
    p.bound_activity  = dx->bound_activity;
    p.binding_window  = dx->binding_window;
    
    /*------------------------------------------------------------------ Binding site type (BS) attributes */
    
    p.TF_tag += dx->TF_tag;
    
    /*------------------------------------------------------------------ Promoter type (P) attributes */
    
    p.basal_expression_level += dx->basal_expression_level;
    p.basal_expression_level  = (p.basal_expression_level > 0.0 ? p.basal_expression_level : 0.0);
    p.basal_expression_level  = (p.basal_expression_level < 1.0 ? p.basal_expression_level : 1.0);
    
    /*------------------------------------------------------------------ Save point mutation event */
    
    MutationEvent* new_event = new MutationEvent(POINT_MUTATION, pos, dX);
    _replication_report->add_mutation_event(new_event);
  }
  else
  {
    delete dX;
    dX = NULL;
  }
}

/**
 * \brief    Do a permutation
 * \details  Perform a permutation between pos1 and pos2 units
 * \param    const double* mutation_rates
 * \param    size_t pos1
 * \param    size_t pos2
 * \return   \e void
 */
void Genome::do_permutation( const double* mutation_rates, size_t pos1, size_t pos2 )
{
  genetic_unit& p1 = _genetic_sequence->x[pos1];
  genetic_unit& p2 = _genetic_sequence->x[pos2];
  
  /*------------------------------------------------------------------ Global attributes */
  
  if (_prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    genetic_unit_type tmp = p1.type;
    p1.type = p2.type;
    p2.type = tmp;
  }
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  if (_prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    int tmp = p1.s;
    p1.s = p2.s;
    p2.s = tmp;
  }
  if (_prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    int tmp = p1.p;
    p1.p = p2.p;
    p2.p = tmp;
  }
  if (_prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    double tmp = p1.kcat;
    p1.kcat = p2.kcat;
    p2.kcat = tmp;
  }
  if (_prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    double tmp = p1.kcat_km_ratio;
    p1.kcat_km_ratio = p2.kcat_km_ratio;
    p2.kcat_km_ratio = tmp;
  }
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  if (_prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    int tmp = p1.BS_tag;
    p1.BS_tag = p2.BS_tag;
    p2.BS_tag = tmp;
  }
  if (_prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    int tmp = p1.coE_tag;
    p1.coE_tag = p2.coE_tag;
    p2.coE_tag = tmp;
  }
  if (_prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    bool tmp = p1.free_activity;
    p1.free_activity = p2.free_activity;
    p2.free_activity = tmp;
  }
  if (_prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    bool tmp = p1.bound_activity;
    p1.bound_activity = p2.bound_activity;
    p2.bound_activity = tmp;
  }
  if (_prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    size_t tmp = p1.binding_window;
    p1.binding_window = p2.binding_window;
    p2.binding_window = tmp;
  }
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  if (_prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    int tmp = p1.TF_tag;
    p1.TF_tag = p2.TF_tag;
    p2.TF_tag = tmp;
  }
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  if (_prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    double tmp = p1.basal_expression_level;
    p1.basal_expression_level = p2.basal_expression_level;
    p2.basal_expression_level = tmp;
  }
}

/**
 * \brief    Create a default genetic sequence
 * \details  --
 * \param    void
 * \return   \e void
 */
void Genome::create_genetic_sequence( void )
{
  _genetic_sequence              = new genetic_sequence;
  _genetic_sequence->x           = new genetic_unit[GENOME_BUFFER];
  _genetic_sequence->size        = 0;
  _genetic_sequence->buffer_size = GENOME_BUFFER;
}

/**
 * \brief    Copy a genetic sequence
 * \details  --
 * \param    const genetic_sequence* model
 * \return   \e void
 */
void Genome::copy_genetic_sequence( const genetic_sequence* model )
{
  _genetic_sequence              = new genetic_sequence;
  _genetic_sequence->x           = new genetic_unit[model->buffer_size];
  memcpy(_genetic_sequence->x, model->x, sizeof(genetic_unit)*model->buffer_size);
  _genetic_sequence->size        = model->size;
  _genetic_sequence->buffer_size = model->buffer_size;
}

/**
 * \brief    Delete the genetic sequence
 * \details  --
 * \param    void
 * \return   \e void
 */
void Genome::delete_genetic_sequence( void )
{
  delete[] _genetic_sequence->x;
  _genetic_sequence->x = NULL;
  delete _genetic_sequence;
  _genetic_sequence = NULL;
}

/**
 * \brief    Load the genetic sequence from backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void Genome::load_genetic_sequence( gzFile backup_file )
{
  _genetic_sequence = new genetic_sequence;
  gzread( backup_file, &_genetic_sequence->size,        sizeof(_genetic_sequence->size) );
  gzread( backup_file, &_genetic_sequence->buffer_size, sizeof(_genetic_sequence->buffer_size) );
  _genetic_sequence->x = new genetic_unit[_genetic_sequence->buffer_size];
  for (size_t i = 0; i < _genetic_sequence->size; i++)
  {
    load_genetic_unit(backup_file, _genetic_sequence->x[i]);
  }
}

/**
 * \brief    Save the genetic sequence in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void Genome::save_genetic_sequence( gzFile backup_file )
{
  gzwrite( backup_file, &_genetic_sequence->size,        sizeof(_genetic_sequence->size) );
  gzwrite( backup_file, &_genetic_sequence->buffer_size, sizeof(_genetic_sequence->buffer_size) );
  for (size_t i = 0; i < _genetic_sequence->size; i++)
  {
    save_genetic_unit(backup_file, _genetic_sequence->x[i]);
  }
}

/**
 * \brief    Load a genetic unit in backup file
 * \details  --
 * \param    gzFile backup_file
 * \param    genetic_unit& unit
 * \return   \e void
 */
void Genome::load_genetic_unit( gzFile backup_file, genetic_unit& unit )
{
  /*------------------------------------------------------------------ global attributes */
  
  gzread( backup_file, &unit.type,              sizeof(unit.type) );
  gzread( backup_file, &unit.identifier,        sizeof(unit.identifier) );
  gzread( backup_file, &unit.parent_identifier, sizeof(unit.parent_identifier) );
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  gzread( backup_file, &unit.s,             sizeof(unit.s) );
  gzread( backup_file, &unit.p,             sizeof(unit.p) );
  gzread( backup_file, &unit.kcat,          sizeof(unit.kcat) );
  gzread( backup_file, &unit.kcat_km_ratio, sizeof(unit.kcat_km_ratio) );
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  gzread( backup_file, &unit.BS_tag,         sizeof(unit.BS_tag) );
  gzread( backup_file, &unit.coE_tag,        sizeof(unit.coE_tag) );
  gzread( backup_file, &unit.free_activity,  sizeof(unit.free_activity) );
  gzread( backup_file, &unit.bound_activity, sizeof(unit.bound_activity) );
  gzread( backup_file, &unit.binding_window, sizeof(unit.binding_window) );
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  gzread( backup_file, &unit.TF_tag, sizeof(unit.TF_tag) );
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  gzread( backup_file, &unit.basal_expression_level, sizeof(unit.basal_expression_level) );
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  gzread( backup_file, &unit.functional, sizeof(unit.functional) );
}

/**
 * \brief    Save a genetic unit in backup file
 * \details  --
 * \param    gzFile backup_file
 * \param    genetic_unit& unit
 * \return   \e void
 */
void Genome::save_genetic_unit( gzFile backup_file, genetic_unit& unit )
{
  /*------------------------------------------------------------------ global attributes */
  
  gzwrite( backup_file, &unit.type,              sizeof(unit.type) );
  gzwrite( backup_file, &unit.identifier,        sizeof(unit.identifier) );
  gzwrite( backup_file, &unit.parent_identifier, sizeof(unit.parent_identifier) );
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  gzwrite( backup_file, &unit.s,             sizeof(unit.s) );
  gzwrite( backup_file, &unit.p,             sizeof(unit.p) );
  gzwrite( backup_file, &unit.kcat,          sizeof(unit.kcat) );
  gzwrite( backup_file, &unit.kcat_km_ratio, sizeof(unit.kcat_km_ratio) );
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  gzwrite( backup_file, &unit.BS_tag,         sizeof(unit.BS_tag) );
  gzwrite( backup_file, &unit.coE_tag,        sizeof(unit.coE_tag) );
  gzwrite( backup_file, &unit.free_activity,  sizeof(unit.free_activity) );
  gzwrite( backup_file, &unit.bound_activity, sizeof(unit.bound_activity) );
  gzwrite( backup_file, &unit.binding_window, sizeof(unit.binding_window) );
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  gzwrite( backup_file, &unit.TF_tag, sizeof(unit.TF_tag) );
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  gzwrite( backup_file, &unit.basal_expression_level, sizeof(unit.basal_expression_level) );
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  gzwrite( backup_file, &unit.functional, sizeof(unit.functional) );
}

/**
 * \brief    Draw a random genetic unit
 * \details  --
 * \param    genetic_unit& unit
 * \return   \e void
 */
void Genome::draw_random_genetic_unit( genetic_unit& unit )
{
  /*------------------------------------------------------------------ Global attributes */
  
  double probas[5]       = {1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0};
  unit.type              = (genetic_unit_type)_prng->roulette_wheel(probas, 1.0, 5);
  unit.identifier        = 0;
  unit.parent_identifier = 0;
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  unit.s             = _prng->uniform(1, 100);
  unit.p             = _prng->uniform(1, 100);
  unit.kcat          = pow(10.0, _prng->uniform()*(KCAT_MAX_LOG-KCAT_MIN_LOG)+KCAT_MIN_LOG);
  unit.kcat_km_ratio = pow(10.0, _prng->uniform()*(KCAT_KM_RATIO_MAX_LOG-KCAT_KM_RATIO_MIN_LOG)+KCAT_KM_RATIO_MIN_LOG);
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  unit.BS_tag         = _prng->uniform(-500, 500);
  unit.coE_tag        = _prng->uniform(1, 100);
  unit.free_activity  = (_prng->uniform() < 0.5 ? true : false);
  unit.bound_activity = (_prng->uniform() < 0.5 ? true : false);
  unit.binding_window = _parameters->get_transcription_factor_binding_window();
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  unit.TF_tag = _prng->uniform(-500, 500);
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  unit.basal_expression_level = _prng->uniform();
}

/**
 * \brief    Increase buffer size relatively to the new genome size
 * \details  --
 * \param    size_t new_size
 * \return   \e void
 */
void Genome::increase_buffer_size( size_t new_size )
{
  assert(new_size <= _genetic_sequence->size*2);
  _genetic_sequence->buffer_size = (new_size/GENOME_BUFFER+1)*GENOME_BUFFER;
  genetic_unit* new_x = new genetic_unit[_genetic_sequence->buffer_size];
  memcpy(new_x, _genetic_sequence->x, sizeof(genetic_unit)*_genetic_sequence->size);
  delete[] _genetic_sequence->x;
  _genetic_sequence->x = new_x;
}

/**
 * \brief    Decrease buffer size relatively to the genome size
 * \details  --
 * \param    void
 * \return   \e void
 */
void Genome::decrease_buffer_size( void )
{
  if (_genetic_sequence->buffer_size > GENOME_BUFFER)
  {
    _genetic_sequence->buffer_size = (_genetic_sequence->size/GENOME_BUFFER+1)*GENOME_BUFFER;
    assert(_genetic_sequence->buffer_size >= _genetic_sequence->size);
    genetic_unit* new_x = new genetic_unit[_genetic_sequence->buffer_size];
    memcpy(new_x, _genetic_sequence->x, sizeof(genetic_unit)*_genetic_sequence->size);
    delete[] _genetic_sequence->x;
    _genetic_sequence->x = new_x;
  }
}

/**
 * \brief    Check if the genome size reaches the maximum genome size
 * \details  --
 * \param    void
 * \return   \e void
 */
void Genome::check_genome_size( void )
{
  if (_genetic_sequence->size > MAXIMUM_GENOME_SIZE)
  {
    delete[] _genetic_sequence->x;
    _genetic_sequence->x = new genetic_unit[GENOME_BUFFER];
    _genetic_sequence->size = 0;
    _genetic_sequence->buffer_size = GENOME_BUFFER;
  }
}
