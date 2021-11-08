#!/bin/bash

#****************************************************************************
# Evo2Sim (Evolution of Evolution Simulator)
# -------------------------------------------
# Digital evolution model dedicated to
# bacterial in silico experimental evolution.
#
# Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon
# Web: https://github.com/charlesrocabert/Evo2Sim
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#****************************************************************************

make clean
rm -rf ../build/bin/*
rm -rf CMakeFiles
rm cmake_install.cmake
rm CMakeCache.txt
rm Config.h
rm libEvo2Sim.a
rm Makefile
