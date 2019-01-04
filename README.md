<h1 align="center">Evo<sup>2</sup>Sim</h1>
<p align="center">
    <img src="logo/logo_evo2sim_small.png">
    <br/>
    <em>Evolution of Evolution Simulator</em>
    <br/><br/>
    <a href="https://github.com/charlesrocabert/Evo2Sim/releases/latest"><img src="https://img.shields.io/badge/version- 1.0.2-green.svg" /></a>&nbsp;<a href="https://github.com/charlesrocabert/Evo2Sim/releases/latest"><img src="https://img.shields.io/badge/build-passing-green.svg" /></a>&nbsp;<a href="https://www.gnu.org/licenses/gpl-3.0"><img src="https://img.shields.io/badge/license-GPL v3-blue.svg" /></a>&nbsp;
</p>

Evo<sup>2</sup>Sim is a multi-scale evolutionary model which takes into account the (ultra-fast) dynamics of metabolic networks, the (fast) dynamics of gene regulatory networks, the (medium to slow) dynamics of resources in the ecosystem and the (slow) evolutionary dynamics of genes and genome structure. Evo<sup>2</sup>Sim is designed to perform <em>in silico</em> experimental evolution experiments, that mimick real experimental evolution protocols.

The software is equipped with the whole tool case of experimental setups, competition assays, phylogenetic analysis. Running simulations can be analyzed thanks to a detailed HTML graphical viewer, tracking and displaying on the fly every simulation events, statistics, phylogeny, and ecological interactions.

Evo<sup>2</sup>Sim was first developed during the EvoEvo project (http://www.evoevo.eu/), funded by the European Commission (FP7-ICT-610427, FET Proactive: Evolving Living Technologies).

## Table of contents
1. [Publications](#publications)
2. [Copyright](#copyright)
3. [License](#license)
4. [Community](#community)
5. [Installation instructions](#installation_instructions)
6. [Typical usage](#typical_usage)

## 1. Publications <a name="publications"></a>

&bull; Rocabert, C., Knibbe, C., Consuegra, J., Schneider, D., & Beslon, G. (2017, Sept.). Environmental seasonality drives digital populations towards stable cross-feeding. <em>Proceedings of the 14th European Conference on Artificial Life (ECAL)</em> (Villeurbanne, France). (https://hal.archives-ouvertes.fr/hal-01569093/)
<br />
&bull; Rocabert, C., Knibbe, C., Consuegra, J., Schneider, D., & Beslon, G. (2017). Beware batch culture: Seasonality and niche construction predicted to favor bacterial adaptive diversification. <em>PLoS computational biology</em>, 13(3), e1005459. (https://doi.org/10.1371/journal.pcbi.1005459)
<br />
&bull; Rocabert, C., Knibbe, C., & Beslon, G. (2015, Jul.). Towards an Integrated Evolutionary Model to Study Evolution of Evolution. <em>Proceedings of the EvoEvo Workshop, Satellite workshop of ECAL 2015</em> (York, UK). (https://hal.inria.fr/hal-01252796/)

## 2. Copyright <a name="copyright"></a>

Copyright &copy; 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon.
All rights reserved

## 3. License <a name="license"></a>

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

## 4. Community <a name="community"></a>

Evo<sup>2</sup>Sim was first developed by Charles Rocabert, Carole Knibbe and Guillaume Beslon, under the EvoEvo project (2013-2016). The list of contributors is displayed in [AUTHORS.md](AUTHORS.md). You shall find more details on EvoEvo community at http://www.evoevo.eu/community/.

## 5. Installation instructions <a name="installation_instructions"></a>

Installation instructions are also available in the [User Manual](doc/user_manual/user_manual.pdf).

Download the latest release of Evo<sup>2</sup>Sim, and save it to a directory of your choice. Open a terminal and use the <code>cd</code> command to navigate to this directory. Then follow the steps below to compile and build the executables.

### 5.1. Supported platforms
Evo<sup>2</sup>Sim software has been successfully tested on Ubuntu 12.04 LTS, Ubuntu 14.04 LTS, OSX 10.9.5 (Maverick) and OSX 10.10.1 (Yosemite).

### 5.2. Required dependencies
* A C++ compiler (GCC, LLVM, ...)
* CMake (command line version)
* zlib
* GSL
* CBLAS
* TBB
* R (packages ape and RColorBrewer are required)

### 5.3. Optional dependencies (for graphical outputs)
* X11 (or XQuartz on latest OSX versions)
* SFML 2
* matplotlib (this python library is needed for the script track_cell.py, see below)

### 5.4. HTML viewer dependencies
* Javascript must be activated in your favorite internet browser

Note, however, that Evo<sup>2</sup>Sim can be compiled without graphical outputs, and hence no needs for X and SFML libraries (see compilation instructions below for more information). This option is useful if you want to run Evo<sup>2</sup>Sim on a computer cluster, for example.

### 5.5. Software compilation

#### User mode
To compile Evo<sup>2</sup>Sim, run the following instructions on the command line:

    cd cmake/

and

    bash make.sh

To gain performances during large experimental protocols, or on computer cluster, you should compile the software without graphical outputs:

    bash make_no_graphics.sh

#### Debug mode
To compile the software in DEBUG mode, use <code>make_debug.sh</code> script instead of <code>make.sh</code>:

    bash make_debug.sh

When Evo<sup>2</sup>Sim is compiled in DEBUG mode, many tests are computed on the fly during a simulation (<em>e.g.</em> integrity tests on phylogenetic trees, or on the ODE solver . . . ). For this reason, this mode should only be used for test or development phases. Moreover, unitary and integrated tests must be ran in DEBUG mode (see below).

#### Executable files emplacement
Binary executable files are in <code>build/bin</code> folder.

## 6. Typical usage <a name="typical_usage"></a>

Evo<sup>2</sup>Sim includes three main executables (<code>evo2sim_create</code>, <code>evo2sim_bootstrap</code> and <code>evo2sim_run</code>), and a set of executables dedicated to post-treatments, data recovery or tests.

Everything in Evo<sup>2</sup>Sim relies on an ad-hoc file organization where all the data for a simulation is stored: populations in the <code>population</code> directory, environments in <code>environment</code>, phylogenetic and lineage trees in <code>tree</code> and so on. It is not recommended to manually modify these files since this may cause some inconsistency leading to undefined behavior. Besides, most of these files are compressed.

Open a terminal and use the <code>cd</code> command to navigate to Evo<sup>2</sup>Sim directory. A typical parameters file is provided in the folder <code>example</code> (an exhaustive description of the parameters is available in the [User Manual](doc/user_manual/user_manual.pdf)). Navigate to this folder using the <code>cd</code> command. Then follow the steps below for a first usage of the software.

### 6.1. Create a simulation
Create a fresh simulation from the parameters file (by default <code>parameters.txt</code>):

    ../build/bin/evo2sim_create

Several folders have been created. They mainly contain simulation backups (population, environment, trees, parameters, ...). Additional files and folders have also been created:
* <code>version.txt</code>: this file indicates the version of the software. This information is useful to ensure that the code version is compatible with the backup files (<em>e.g.</em>, in case of post-treatments).
* <code>track_cell.py</code>: when executed, this python script displays on the fly the internal protein and metabolic concentrations of the cell at position 0 × 0 on the grid. This script is useful to get an idea of internal cell’s dynamics (metabolic fluxes, regulation,).
* <code>viewer</code> folder: the viewer is central to the usage of Evo2Sim (see [User Manual](doc/user_manual/user_manual.pdf)). To access the viewer, open the html page <code>viewer/viewer.html</code> in an internet browser.

### 6.2. Generate viable initial conditions with a bootstrap
Alternatively to the <code>evo2sim_create</code> executable, use a bootstrap to find a simulation with good initial properties from the parameters file:

    ../build/bin/evo2sim_bootstrap

A fresh simulation with an updated parameters file will be automatically created if a suitable seed is found.

### 6.3. Run a simulation
In Evo<sup>2</sup>Sim, running a simulation necessitates to load it from backup files. Here, we will run a simulation from freshly created backups (see above):

    ../build/bin/evo2sim_run -b 0 -t 10000 -g


with <code>-b</code> the date of the backup, here 0 (fresh simulation), <code>-t</code> the simulation time, here 10,000 time-steps. Option <code>-g</code> activates the graphical output (does not work if the software has been compiled with the no-graphics option). At any moment during the simulation, you can take a closer look at the evolution of the system by opening <code>viewer/viewer.html</code> in an internet browser. You can track internal cell’s dynamics by executing the script <code>track_cell.py</code>.

Other main executables are described in the [User Manual](doc/user_manual/user_manual.pdf) (section “Main executables description”). You can also obtain help by running the executable with the <code>-h</code> option (<em>e.g.</em> <code>evo2sim_create -h</code>)



