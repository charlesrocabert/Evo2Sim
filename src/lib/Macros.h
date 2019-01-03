
/**
 * \file      Macros.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Definition of macros
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

#ifndef __EVOEVO__Macros__
#define __EVOEVO__Macros__


/*-------------------------------------*/
/* 1) SEED GENERATION MACROS           */
/*-------------------------------------*/

#define MAXIMUM_SEED 100000000 /*!< Maximum drawable seed */

/*-------------------------------------*/
/* 2) MEMORY MANAGEMENT MACROS         */
/*-------------------------------------*/

#define GENOME_BUFFER              500 /*!< Genome buffer size                        */
#define INHERITED_PROTEINS_BUFFER  500 /*!< Inherited proteins buffer size            */
#define SPECIES_LIST_BUFFER        5   /*!< Species list buffer size                  */
#define GRN_REACTIONS_BUFFER       20  /*!< Number of GRN reactions buffer size       */
#define METABOLIC_REACTIONS_BUFFER 10  /*!< Number of metabolic reactions buffer size */

/*-------------------------------------*/
/* 3) GENOME SIZE MACROS               */
/*-------------------------------------*/

#define MAXIMUM_GENOME_SIZE 10000 /* Maximum genome size */

/*-------------------------------------*/
/* 4) MUTATION RATES MACROS            */
/*-------------------------------------*/

#define NUMBER_OF_MUTATION_RATES 15 /*!< Number of mutation rates */

/*-------------------------------------*/
/* 5) HGT MACROS                       */
/*-------------------------------------*/

#define HGT_MIN_SIZE 1  /*!< HGT minimum size */
#define HGT_MAX_SIZE 10 /*!< HGT maximum size */

/*-------------------------------------*/
/* 6) ODE SOLVER MACROS                */
/*-------------------------------------*/

#define ERR_ABS 1e-30 /*!< ODE solver absolute precision */
#define ERR_REL 1e-02 /*!< ODE solver relative precision */
#define H_INIT  1e-01 /*!< ODE solver initial timestep   */
#define H_MAX   5     /*!< ODE solver maximum timestep   */

/*-------------------------------------*/
/* 7) MICHAELIS MENTEN EQUATION MACROS */
/*-------------------------------------*/

#define KCAT_MIN_LOG          -3.0 /*!< Minimum log10 boundarie of Kcat          */
#define KCAT_MAX_LOG          -1.0 /*!< Maximum log10 boundarie of Kcat          */
#define KCAT_KM_RATIO_MIN_LOG -4.0 /*!< Minimum log10 boundarie of Kcat/Km ratio */
#define KCAT_KM_RATIO_MAX_LOG -4.0 /*!< Maximum log10 boundarie of Kcat/Km ratio */

/*-------------------------------------*/
/* 8) CONCENTRATIONS EVOLUTION MACROS  */
/*-------------------------------------*/

#define MINIMUM_CONCENTRATION                  1e-18 /*!< Minimum concentration threshold                    */
#define MINIMUM_SCORE                          1e-03 /*!< Minimum sustainable score to survive               */
#define MINIMUM_HERITABLE_ENZYME_CONCENTRATION 1e-06 /*!< All enzymes below this threshold are not heritable */

/*-------------------------------------*/
/* 9) PRIME NUMBERS LIST               */
/*-------------------------------------*/

#define PRIME_NUMBERS_LIST_LENGTH  5000

/*-------------------------------------*/
/* 10) GRAPHICS MACROS                 */
/*-------------------------------------*/

#define FRAMERATE      0  /*!< Framerate of the graphic display     */
#define CELL_SCALE     5  /*!< Cell's scale in the graphic display  */
#define CELL_SPACE     1  /*!< Space displayed between each cell    */
#define GRADIENT_SCALE 30 /*!< Gradient's scale                     */
#define GRADIENT_SIZE  10 /*!< Gradient's size (nb of color points) */
#define TEXT_SCALE     50 /*!< Text's scale                         */
#define SPAN           5  /*!< Span of the render window            */


#endif /* defined(__EVOEVO__Macros__) */
