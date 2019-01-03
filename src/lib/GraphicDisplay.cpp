
/**
 * \file      GraphicDisplay.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     GraphicDisplay class definition
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

#include "GraphicDisplay.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters
 * \param    Simulation* simulation
 * \param    sf::RenderWindow* pop_window
 * \param    sf::RenderWindow* env_window
 * \param    std::string font_path
 * \return   \e void
 */
GraphicDisplay::GraphicDisplay( Parameters* parameters, Simulation* simulation, sf::RenderWindow* pop_window, sf::RenderWindow* env_window, std::string font_path )
{
  /*------------------------------------------------------------------ main parameters */
  
  _parameters = parameters;
  _simulation = simulation;
  _pop_window = pop_window;
  _env_window = env_window;
  
  _lattice = new sf::RectangleShape(sf::Vector2f(_parameters->get_width()*CELL_SCALE+(_parameters->get_width()-1)*CELL_SPACE, _parameters->get_height()*CELL_SCALE+(_parameters->get_height()-1)*CELL_SPACE));
  _lattice->setFillColor(sf::Color::Black);
  _lattice->setPosition(SPAN, SPAN);
  _cell = new sf::RectangleShape(sf::Vector2f(CELL_SCALE, CELL_SCALE));
  
  build_gradients(GRADIENT_SIZE);
  
  _font_path = font_path;
  if (!_font.loadFromFile(_font_path))
  {
    printf("ERROR: Font not found in GraphicDisplay constructor.\n");
    exit(EXIT_FAILURE);
  }
  _text_up.setFont(_font);
  _text_up.setCharacterSize(10);
  _text_up.setColor(sf::Color::White);
  _text_up.setPosition(_lattice->getSize().x+3*SPAN+GRADIENT_SCALE, SPAN-SPAN/2);
  
  _text_down.setFont(_font);
  _text_down.setCharacterSize(10);
  _text_down.setColor(sf::Color::White);
  _text_down.setString(std::to_string(0.0));
  _text_down.setPosition(_lattice->getSize().x+3*SPAN+GRADIENT_SCALE, _lattice->getSize().y-SPAN/2);
  
  /*------------------------------------------------------------------ events management */
  
  _color_scheme = 2;
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
GraphicDisplay::~GraphicDisplay( void )
{
  delete _lattice;
  _lattice = NULL;
  delete _cell;
  _cell = NULL;
  delete _pop_gradient;
  _pop_gradient = NULL;
  delete _env_gradient;
  _env_gradient = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Display windows
 * \details  --
 * \param    void
 * \return   \e bool
 */
bool GraphicDisplay::display( void )
{
  if (_pop_window->isOpen() && _env_window->isOpen())
  {
    /* A) Display graphics ------------------------*/
    _pop_window->clear(sf::Color(50,50,50));
    _env_window->clear(sf::Color(50,50,50));
    display_population_and_environment();
    _pop_window->display();
    _env_window->display();
    
    _pop_window->setTitle("Population ("+std::to_string(_simulation->get_population()->get_time())+" timesteps)");
    _env_window->setTitle("Environment ("+std::to_string(_simulation->get_population()->get_time())+" timesteps)");
    
    /* B) Events management -----------------------*/
    sf::Event event;
    while (_pop_window->pollEvent(event))
    {
      /* Window's close event */
      if (event.type == sf::Event::Closed)
      {
        _pop_window->close();
        return false;
      }
      if (event.type == sf::Event::KeyPressed)
      {
        if (event.key.code == sf::Keyboard::Left)
        {
          _color_scheme = (_color_scheme+1+3)%3;
          if (_color_scheme == 0)
          {
            printf("> Switch to trophic group colouring\n");
          }
          else if (_color_scheme == 1)
          {
            printf("> Switch to cell's score colouring\n");
          }
          else if (_color_scheme == 2)
          {
            printf("> Switch to trophic level colouring\n");
          }
        }
        else if (event.key.code == sf::Keyboard::Right)
        {
          _color_scheme = (_color_scheme-1+3)%3;
          if (_color_scheme == 0)
          {
            printf("> Switch to trophic group colouring\n");
          }
          else if (_color_scheme == 1)
          {
            printf("> Switch to cell's score colouring\n");
          }
          else if (_color_scheme == 2)
          {
            printf("> Switch to trophic level colouring\n");
          }
        }
      }
    }
    while (_env_window->pollEvent(event))
    {
      /* Window's close event */
      if (event.type == sf::Event::Closed)
      {
        _env_window->close();
        return false;
      }
    }
    return true;
  }
  return false;
}

/**
 * \brief    Display population and environment grids
 * \details  --
 * \param    void
 * \return   \e bool
 */
bool GraphicDisplay::display_population_and_environment( void )
{
  _pop_window->draw(*_lattice);
  _env_window->draw(*_lattice);
  double min_score  = _simulation->get_min_score();
  double max_score  = _simulation->get_max_score();
  double min_amount = _simulation->get_environment()->get_min_amount();
  double max_amount = _simulation->get_environment()->get_max_amount();
  size_t xcount     = 0;
  size_t ycount     = 0;
  for (size_t x = 0; x < _parameters->get_width(); x++)
  {
    ycount = 0;
    for (size_t y = 0; y < _parameters->get_height(); y++)
    {
      Cell& cell        = *_simulation->get_population()->get_cell(x, y);
      double env_amount = _simulation->get_environment()->get_species_list(x, y)->get_amount();
      double score      = cell.get_score();
      double value      = (score-min_score)/(max_score-min_score);
      _cell->setPosition(sf::Vector2f(SPAN+x*CELL_SCALE+xcount, SPAN+y*CELL_SCALE+ycount));
      if (cell.isAlive())
      {
        if (_color_scheme == 0)
        {
          _cell->setFillColor(sf::Color(cell.get_red_color(), cell.get_green_color(), cell.get_blue_color()));
        }
        else if (_color_scheme == 1)
        {
          _cell->setFillColor(get_color(value < 1.0 ? value : 1.0));
        }
        else if (_color_scheme == 2)
        {
          _cell->setFillColor(get_color_by_trophic_level(cell.get_trophic_level()));
        }
        _pop_window->draw(*_cell);
      }
      value = (env_amount-min_amount)/(max_amount-min_amount);
      _cell->setFillColor(get_color(value < 1.0 ? value : 1.0));
      _env_window->draw(*_cell);
      ycount += CELL_SPACE;
    }
    xcount += CELL_SPACE;
  }
  _pop_window->draw(_pop_gradient, GRADIENT_SIZE*2, sf::TrianglesStrip);
  _env_window->draw(_env_gradient, GRADIENT_SIZE*2, sf::TrianglesStrip);
  
  _text_up.setString(std::to_string(max_score));
  _pop_window->draw(_text_up);
  
  _text_down.setString(std::to_string(min_score));
  _pop_window->draw(_text_down);
  
  _text_up.setString(std::to_string(max_amount));
  _env_window->draw(_text_up);
  
  _text_down.setString(std::to_string(min_amount));
  _env_window->draw(_text_down);
  
  return true;
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Build color gradients
 * \details  --
 * \param    int size
 * \return   \e void
 */
void GraphicDisplay::build_gradients( size_t size )
{
  /*-----------------------------*/
  /* 1) Create vertices array    */
  /*-----------------------------*/
  _pop_gradient = new sf::Vertex[2*size];
  _env_gradient = new sf::Vertex[2*size];
  
  /*-----------------------------*/
  /* 2) Get lattice size         */
  /*-----------------------------*/
  sf::Vector2f vec = _lattice->getSize();
  
  /*-----------------------------*/
  /* 3) Initialize coordinates   */
  /*-----------------------------*/
  double x1 = vec.x+2*SPAN;
  double x2 = x1+GRADIENT_SCALE;
  
  /*-----------------------------*/
  /* 4) Initialize gradient step */
  /*-----------------------------*/
  double step = vec.y/(size-1);
  
  /*-----------------------------*/
  /* 5) Fill color points        */
  /*-----------------------------*/
  double y = SPAN;
  double value = 1.0;
  for (size_t i = 0; i < 2*size; i += 2)
  {
    _pop_gradient[i]   = sf::Vertex(sf::Vector2f(x1, y), get_color(value));
    _pop_gradient[i+1] = sf::Vertex(sf::Vector2f(x2, y), get_color(value));
    
    _env_gradient[i]   = sf::Vertex(sf::Vector2f(x1, y), get_color(value));
    _env_gradient[i+1] = sf::Vertex(sf::Vector2f(x2, y), get_color(value));
    
    y += step;
    value -= step/vec.y;
  }
}

/**
 * \brief    Get color from value between 0 and 1
 * \details  --
 * \param    double value
 * \return   \e sf::Color
 */
sf::Color GraphicDisplay::get_color( double value )
{
  value = pow(value, 2.0);
  if (value < 0.2)
  {
    return sf::Color((0.1) * 255.0, (0.2) * 255.0, (0.3) * 255.0);
  }
  else if (value >= 0.2 && value < 0.7)
  {
    return sf::Color((0.1 + ((value-0.2) / 0.5)*0.9) * 255.0, (0.2 - ((value-0.2) / 0.5)*0.2) * 255.0, (0.3 + ((value-0.2) / 0.5)*0.7) * 255.0);
  }
  else if (value >= 0.7 && value < 0.9)
  {
    return sf::Color((1.0) * 255.0, ((value-0.7) / 0.2) * 255.0, (1.0 - (value-0.7) / 0.2) * 255.0);
  }
  else if (value >= 0.9)
  {
    return sf::Color((1.0) * 255.0, (1.0) * 255.0, ((value-0.9) / 0.1) * 255.0);
  }
  return sf::Color(0.0, 255.0, 0.0);
  /*
  value = pow(value, 3.0);
  double rfactor = 2.0;
  if (value < 0.2)
  {
    double r = (83-34+rfactor)*value+34+rfactor;
    double g = (42-22)*value+22;
    double b = (90-48)*value+48;
    return sf::Color(r, g, b);
  }
  else if (value >= 0.2 && value < 0.4)
  {
    double r = (142-83+rfactor)*value+83+rfactor;
    double g = (73-42)*value+42;
    double b = (121-90)*value+90;
    return sf::Color(r, g, b);
  }
  else if (value >= 0.4 && value < 0.6)
  {
    double r = (190-142+rfactor)*value+142+rfactor;
    double g = (117-73)*value+73;
    double b = (143-121)*value+121;
    return sf::Color(r, g, b);
  }
  else if (value >= 0.6 && value < 0.8)
  {
    double r = (223-190+rfactor)*value+190+rfactor;
    double g = (174-117)*value+117;
    double b = (175-143)*value+143;
    return sf::Color(r, g, b);
  }
  else if (value >= 0.8)
  {
    double r = (247-223+rfactor)*value+223+rfactor;
    double g = (234-174)*value+174;
    double b = (227-175)*value+175;
    return sf::Color(r, g, b);
  }
   */
}

/**
 * \brief    Get color from the trophic level
 * \details  --
 * \param    trophic_level level
 * \return   \e sf::Color
 */
sf::Color GraphicDisplay::get_color_by_trophic_level( trophic_level level )
{
  if (level == NO_LEVEL)
  {
    return sf::Color(109, 109, 109);
  }
  else if (level == LEVEL_0)
  {
    return sf::Color(108, 0, 110);
  }
  else if (level == LEVEL_1)
  {
    return sf::Color(50, 117, 193);
  }
  else if (level == LEVEL_2)
  {
    return sf::Color(116, 169, 44);
  }
  return sf::Color::White;
}
