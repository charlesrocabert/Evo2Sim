
/**
 * \file      GraphicDisplay.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     GraphicDisplay class declaration
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

#ifndef __Evo2Sim__GraphicDisplay__
#define __Evo2Sim__GraphicDisplay__

#include <iostream>
#include <SFML/Graphics.hpp>
#include <cmath>
#include <assert.h>

#include "Macros.h"
#include "Structs.h"
#include "Enums.h"
#include "Parameters.h"
#include "Cell.h"
#include "Population.h"
#include "Simulation.h"


class GraphicDisplay
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  GraphicDisplay( void ) = delete;
  GraphicDisplay( Parameters* parameters, Simulation* simulation, sf::RenderWindow* pop_window, sf::RenderWindow* env_window, std::string font_path );
  GraphicDisplay( const GraphicDisplay& graphicDisplay ) = delete;
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~GraphicDisplay( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  GraphicDisplay& operator=(const GraphicDisplay&) = delete;
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  bool display( void );
  bool display_population_and_environment( void );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  void      build_gradients( size_t step );
  sf::Color get_color( double value );
  sf::Color get_color_by_trophic_level( trophic_level level );
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  
  /*------------------------------------------------------------------ main parameters */
  
  Parameters*         _parameters;   /*!< Simulation parameters      */
  Simulation*         _simulation;   /*!< Simulation                 */
  sf::RenderWindow*   _pop_window;   /*!< Population render window   */
  sf::RenderWindow*   _env_window;   /*!< Environment render window  */
  sf::RectangleShape* _lattice;      /*!< Lattice shape              */
  sf::RectangleShape* _cell;         /*!< Cell shape                 */
  sf::Vertex*         _pop_gradient; /*!< Population color gradient  */
  sf::Vertex*         _env_gradient; /*!< Environment color gradient */
  std::string         _font_path;    /*!< Font src path              */
  sf::Font            _font;         /*!< Font                       */
  sf::Text            _text_up;      /*!< Up text                    */
  sf::Text            _text_down;    /*!< Down text                  */
  
  /*------------------------------------------------------------------ events management */
  
  int _color_scheme; /*!< Color scheme (0:score, 1:trophic level, 2:trophic group) */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/*----------------------------
 * SETTERS
 *----------------------------*/


#endif /* defined(__Evo2Sim__GraphicDisplay__) */
