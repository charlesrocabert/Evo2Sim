# <img src="logo/logo_evo2sim_small.png" width="120"> Evo2Sim 1.0.1

Evo<sup>2</sup>Sim is a new multi-scale model of <em>in silico</em> experimental evolution, the virtual pendant of experimental evolution in laboratory. The software is equipped with the whole tool case of experimental setups, competition assays, phylogenetic analysis, and, most importantly, allowing for evolvable ecological interactions.

Digital organisms with an evolvable genome structure, encoding evolvable genetic regulation and metabolic networks are evolved for tens of thousands of generations in environments mimicking the dynamics of real controlled environments, including chemostat or batch culture.
Evo<sup>2</sup>Sim was first developed during the EvoEvo project (http://www.evoevo.eu/), a FP7-ICT project funded by the European Commission (FP7-ICT-610427).

# License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

# Community

Evo<sup>2</sup>Sim is developed by Charles Rocabert, Carole Knibbe and Guillaume Beslon, under the EvoEvo project. The list of contributors is displayed in [AUTHORS.md](AUTHORS.md). You shall find more details on EvoEvo community on http://www.evoevo.eu/community/.

# Installation instructions

Installation instructions are also available in the [User Manual](doc/user_manual/user_manual.pdf).

Download the latest release of Evo<sup>2</sup>Sim, and save it to a directory of your choice. Open a terminal and use the <code>cd</code> command to navigate to this directory. Then follow the steps below to compile and build the executables.

### 1. Supported platforms
Evo<sup>2</sup>Sim software has been successfully tested on Ubuntu 12.04 LTS, Ubuntu 14.04 LTS, OSX 10.9.5 (Maverick) and OSX 10.10.1 (Yosemite).

### 2. Required dependencies
* A C++ compiler (GCC, LLVM, ...) ∙ CMake (command line version)
* zlib
* GSL
* CBLAS
* TBB
* R (packages ape and RColorBrewer are needed)

### 3. Optional dependencies (for graphical outputs)
* X11 (or XQuartz on latest OSX version)
* SFML 2
* matplotlib (this python library is needed for the script track_cell.py (see below)

