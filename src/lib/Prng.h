
/**
 * \file      Prng.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Prng class declaration
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

#ifndef __EVOEVO__Prng__
#define __EVOEVO__Prng__

#include <iostream>
#include <sstream>
#include <cstring>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <zlib.h>
#include <assert.h>

#include "Macros.h"
#include "Enums.h"


class Prng
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  Prng( void );
  Prng( FILE* backup_file );
  Prng( const Prng& prng );
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~Prng( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  inline void set_seed( unsigned long int seed );
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  double uniform( void );
  int    uniform( int min, int max );
  int    bernouilli( double p );
  size_t binomial( size_t n, double p );
  void   multinomial( unsigned int* draws, double* probas, int N, int K );
  double gaussian( double mu, double sigma );
  int    exponential( double mu );
  int    poisson( double mu );
  int    roulette_wheel( double* probas, double sum, int N );
  void   shuffle( void* base, size_t n, size_t size );
  void   save( FILE* backup_file );
  
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
  gsl_rng* _prng; /*!< Pseudorandom numbers generator */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/*----------------------------
 * SETTERS
 *----------------------------*/

/**
 * \brief    Set prng seed
 * \details  --
 * \param    unsigned long int seed
 * \return   \e double
 */
inline void Prng::set_seed( unsigned long int seed )
{
  gsl_rng_set(_prng, seed);
}


#endif /* defined(__EVOEVO__Prng__) */
