<h1 align="center">Evo<sup>2</sup>Sim</h1>
<p align="center">
    <img src="logo/logo_evo2sim_small.png">
    <br/>
    <em>Evolution of Evolution Simulator</em>
    <br/><br/>
    <a href="https://github.com/charlesrocabert/Evo2Sim/releases/latest"><img src="https://img.shields.io/badge/version- 1.1.0-green.svg" /></a>&nbsp;<a href="https://github.com/charlesrocabert/Evo2Sim/releases/latest"><img src="https://img.shields.io/badge/build-passing-green.svg" /></a>&nbsp;<a href="https://www.gnu.org/licenses/gpl-3.0"><img src="https://img.shields.io/badge/license-GPL v3-blue.svg" /></a>&nbsp;
</p>

<p align="justify">
Evolution is the major source of complexity on Earth, at the origin of all the species we can observe, interact with or breed. On a smaller scale, evolution is at the heart of the adaptation process for many species, in particular micro-organisms (e.g. bacteria, viruses...). Microbial evolution results in the emergence of the species itself, and it also contributes to the organisms' adaptation to perturbations or environmental changes. These organisms are not only organised by evolution, they are also organised to evolve (EvoEvo, see www.evoevo.eu).
</p>

<p align="justify">
Evo<sup>2</sup>Sim is a digital evolution model which takes into account the (ultra-fast) dynamics of metabolic networks, the (fast) dynamics of gene regulatory networks, the (medium to slow) dynamics of resources in the ecosystem and the (slow) evolutionary dynamics of genes and genome structure. Evo<sup>2</sup>Sim is designed to perform <em>In Silico</em> Experimental Evolution (ISEE) experiments, that mimick real experimental evolution protocols.
</p>

<p align="justify">
The software is equipped with the whole tool case of experimental setups, competition assays, phylogenetic analysis. Simulations can be analyzed thanks to the HTML viewer, tracking and displaying on the fly every simulation events, from the phylogeny to ecological interactions.
</p>

<p align="justify">
Evo<sup>2</sup>Sim was first developed during the EvoEvo project (http://www.evoevo.eu/), funded by the European Commission (FP7-ICT-610427, FET Proactive: Evolving Living Technologies).
</p>

<p align="justify">
The richness of the genotype-to-phenotype mapping and of the environmental interactions implemented in Evo<sup>2</sup>Sim makes it a versatile model, with which many modeling questions can be tackled. With Evo<sup>2</sup>Sim, we studied the emergence of a stable polymorphism based on cross-feeding, and in which conditions a stable regulation network could evolve.

<strong>:arrow_right: If you want to try or use this model for research purpose, do not hesitate to contact <a href="mailto:charles[DOT]rocabert[AT]helsinki[DOT]fi">Charles Rocabert</a></strong>.
</p>

## Table of contents
- [Publications](#publications)
- [Copyright](#copyright)
- [License](#license)
- [Community](#community)
- [Download](#download)
- [Evo<sup>2</sup>Sim overview](#overview)
- [Installation instructions](#installation_instructions)
   - [Supported platforms](#supported_platforms)
   - [Required dependencies](#required_dependencies)
   - [Optional dependencies (for graphical outputs)](#optional_dependencies)
   - [HTML viewer dependencies](#viewer_dependencies)
   - [Software compilation](#compilation)
- [Typical usage](#typical_usage)
   - [Create a simulation](#create)
   - [Generate viable initial conditions with a bootstrap](#bootstrap)
   - [Run a simulation](#run)
- [Evo<sup>2</sup>Sim viewer](#viewer)
   - [Population viewer](#population_viewer)
   - [Best lineage viewer](#best_lineage_viewer)
   - [Best individual viewer](#best_individual_viewer)
   - [Environment viewer](#environment_viewer)
   - [Phylogeny viewer](#phylogeny_viewer)
- [Examples](#examples)
   - [Evolution of a stable polymorphism](#stable_polymorphism)
   - [Lactose-operon-like regulation](#lactose_operon)

## Publications <a name="publications"></a>

&bull; Rocabert, C., Knibbe, C., Consuegra, J., Schneider, D., & Beslon, G. (2017, Sept.). Environmental seasonality drives digital populations towards stable cross-feeding. <em>Proceedings of the 14th European Conference on Artificial Life (ECAL)</em> (Villeurbanne, France). (https://hal.archives-ouvertes.fr/hal-01569093/)
<br />
&bull; Rocabert, C., Knibbe, C., Consuegra, J., Schneider, D., & Beslon, G. (2017). Beware batch culture: Seasonality and niche construction predicted to favor bacterial adaptive diversification. <em>PLoS computational biology</em>, 13(3), e1005459. (https://doi.org/10.1371/journal.pcbi.1005459)
<br />
&bull; Rocabert, C., Knibbe, C., & Beslon, G. (2015, Jul.). Towards an Integrated Evolutionary Model to Study Evolution of Evolution. <em>Proceedings of the EvoEvo Workshop, Satellite workshop of ECAL 2015</em> (York, UK). (https://hal.inria.fr/hal-01252796/)

## Copyright <a name="copyright"></a>

Copyright &copy; 2014-2021 <a href="https://github.com/charlesrocabert">Charles Rocabert</a>, <a href="http://caroleknibbe.fr/">Carole Knibbe</a>, <a href="https://perso.liris.cnrs.fr/guillaume.beslon/G._Beslon_home_page/Welcome.html">Guillaume Beslon</a>.
All rights reserved.

## License <a name="license"></a>

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

## Community <a name="community"></a>

Evo<sup>2</sup>Sim was first developed by <a href="https://github.com/charlesrocabert">Charles Rocabert</a>, <a href="http://caroleknibbe.fr/">Carole Knibbe</a> and <a href="https://perso.liris.cnrs.fr/guillaume.beslon/G._Beslon_home_page/Welcome.html">Guillaume Beslon</a>, under the EvoEvo project (2013-2016).

The list of contributors is displayed in [AUTHORS.md](AUTHORS.md).

<table>
    <tr>
        <td><a href="https://www.inria.fr/"><img src="docs/readme_img/logo-inria.jpg" height="50px"></a></td>
        <td><a href="https://www.insa-lyon.fr/"><img src="docs/readme_img/logo-insa.png" height="50px"></a></td>
        <td><a href="https://liris.cnrs.fr/"><img src="docs/readme_img/logo-liris.png" height="50px"></a></td>
        <td><a href="https://www.univ-lyon1.fr/"><img src="docs/readme_img/logo-universite-de-lyon.png" height="50px"></a></td>
        <td><a href="https://ec.europa.eu/"><img src="docs/readme_img/logo-eu.jpg" height="50px"></a></td>
    </tr>
</table>

## Download <a name="download"></a>
Download the <a href="https://github.com/charlesrocabert/Evo2Sim/releases/latest">latest release</a>.

## Evo<sup>2</sup>Sim overview <a name="overview"></a>

<img src="https://github.com/charlesrocabert/Evo2Sim/blob/master/docs/user_manual/figures/general_algorithm.jpg">

### a. Description of the genotype-to-phenotype mapping.
Organisms own a coarse-grained genome made of units. This genome is a circular single-strand sequence, with a unique reading frame. Non coding (NC) units are not functional **(a.1)**. The arrangement of the units on the sequence defines functional regions, where a promoter (P, blue cross) controls the expression of enzyme coding units (E, red circles) or transcription factor coding units (TF, purple squares), thereby allowing for operons (here, one E and one TF). When coding units are expressed **(a.2)**, they contribute to the genetic regulatory network (for TFs) and the metabolic network (for Es). Depending on their attributes, transcription factors bind on binding sites. **(a.3)** If they bind on the enhancer sequence (binding sites flanking the promoter upstream), the promoter activity is up-regulated. If they bind on the operator sequence (binding sites flanking the promoter downstream), the promoter activity is down-regulated. **(a.4)** Metabolites can bind on a transcription factor as co-enzymes, and activate or inhibit it, depending on transcription factor attributes. Enzymes perform metabolic reactions in the cytoplasm **(a.5)**, or pump metabolites in or out **(a.6)**. The score of an organism is computed from its “essential metabolites” (usually the score is the sum of essential metabolite concentrations). Lethal toxicity thresholds are applied to each metabolic concentration and forbid organisms to accumulate resources.

### b. Description of the population and environment levels.
Organisms are placed on a 2D toroidal grid, and compete for resources and space. When an organism dies, it leaves its grid cell empty and organisms in the Moore neighborhood (if any) compete to divide in available space. The competition is based on scores, a minimal threshold being applied on scores to forbid worst organisms to divide. At division, daughters share cytoplasm content (enzymes and metabolites). At death, metabolites from the cytoplasm are released in the local environment, and diffuse on the grid **(b.1)**. On the largest scale, the population evolves on the environment by up-taking, transforming and releasing metabolites. Metabolites then diffuse and are degraded. This strong interaction between the population and the environment allows for the evolution of complex ecological situations, depending on environmental properties **(b.2)**.

## Installation instructions <a name="installation_instructions"></a>

Installation instructions are also available in the [User Manual](doc/user_manual/user_manual.pdf).

Download the latest release of Evo<sup>2</sup>Sim, and save it to a directory of your choice. Open a terminal and use the <code>cd</code> command to navigate to this directory. Then follow the steps below to compile and build the executables.

### Supported platforms <a name="supported_platforms"></a>
Evo<sup>2</sup>Sim software has been successfully tested on Ubuntu 12.04 LTS, Ubuntu 14.04 LTS, OSX 10.9.5 (Maverick) and OSX 10.10.1 (Yosemite).

### Required dependencies <a name="required_dependencies"></a>
* A C++ compiler (GCC, LLVM, ...)
* CMake (command line version)
* zlib
* GSL
* CBLAS
* TBB
* R (packages ape and RColorBrewer are required)

### Optional dependencies (for graphical outputs) <a name="optional_dependencies"></a>
* X11 (or XQuartz on latest OSX versions)
* SFML 2
* matplotlib (this python library is needed for the script track_cell.py, see below)

### HTML viewer dependencies <a name="viewer_dependencies"></a>
* Javascript must be activated in your favorite internet browser

Note, however, that Evo<sup>2</sup>Sim can be compiled without graphical outputs, and hence no needs for X and SFML libraries (see compilation instructions below for more information). This option is useful if you want to run Evo<sup>2</sup>Sim on a computer cluster, for example.

### Software compilation <a name="compilation"></a>

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

## Typical usage <a name="typical_usage"></a>

Evo<sup>2</sup>Sim includes three main executables (<code>evo2sim_create</code>, <code>evo2sim_bootstrap</code> and <code>evo2sim_run</code>), and a set of executables dedicated to post-treatments, data recovery or tests.

Everything in Evo<sup>2</sup>Sim relies on an ad-hoc file organization where all the data for a simulation is stored: populations in the <code>population</code> directory, environments in <code>environment</code>, phylogenetic and lineage trees in <code>tree</code> and so on. It is not recommended to manually modify these files since this may cause some inconsistency leading to undefined behavior. Besides, most of these files are compressed.

Open a terminal and use the <code>cd</code> command to navigate to Evo<sup>2</sup>Sim directory. A typical parameters file is provided in the folder <code>example</code> (an exhaustive description of the parameters is available in the [User Manual](doc/user_manual/user_manual.pdf)). Navigate to this folder using the <code>cd</code> command. Then follow the steps below for a first usage of the software.

### Create a simulation <a name="create"></a>
Create a fresh simulation from the parameters file (by default <code>parameters.txt</code>):

    ../build/bin/evo2sim_create

Several folders have been created. They mainly contain simulation backups (population, environment, trees, parameters, ...). Additional files and folders have also been created:
* <code>version.txt</code>: this file indicates the version of the software. This information is useful to ensure that the code version is compatible with the backup files (<em>e.g.</em>, in case of post-treatments).
* <code>track_cell.py</code>: when executed, this python script displays on the fly the internal protein and metabolic concentrations of the cell at position 0 × 0 on the grid. This script is useful to get an idea of internal cell’s dynamics (metabolic fluxes, regulation,).
* <code>viewer</code> folder: the viewer is central to the usage of Evo2Sim (see [User Manual](doc/user_manual/user_manual.pdf)). To access the viewer, open the html page <code>viewer/viewer.html</code> in an internet browser.

### Generate viable initial conditions with a bootstrap <a name="bootstrap"></a>
Alternatively to the <code>evo2sim_create</code> executable, use a bootstrap to find a simulation with good initial properties from the parameters file:

    ../build/bin/evo2sim_bootstrap

A fresh simulation with an updated parameters file will be automatically created if a suitable seed is found.

### Run a simulation <a name="run"></a>
In Evo<sup>2</sup>Sim, running a simulation necessitates to load it from backup files. Here, we will run a simulation from freshly created backups (see above):

    ../build/bin/evo2sim_run -b 0 -t 10000 -g


with <code>-b</code> the date of the backup, here 0 (fresh simulation), <code>-t</code> the simulation time, here 10,000 time-steps. Option <code>-g</code> activates the graphical output (does not work if the software has been compiled with the no-graphics option). At any moment during the simulation, you can take a closer look at the evolution of the system by opening <code>viewer/viewer.html</code> in an internet browser. You can track internal cell’s dynamics by executing the script <code>track_cell.py</code>.

Other main executables are described in the [User Manual](doc/user_manual/user_manual.pdf) (section “Main executables description”). You can also obtain help by running the executable with the <code>-h</code> option (<em>e.g.</em> <code>evo2sim_create -h</code>)

## Evo<sup>2</sup>Sim viewer <a name="viewer"></a>

Evo<sup>2</sup>Sim provides a HTML viewer displaying a very complete set of live statistics. Each new simulation owns a dedicated viewer, that is frequently updated on the fly (by default, every 500 simulation time-steps). This viewer has been developed using `Bootstrap`, `DyGraph`, `CytoscapeJS`, `ChartJS` and `JQuery`.

To access the viewer from a simulation folder, simply open the page `viewer/viewer.html` in your favorite internet browser (Javascript must be enabled).
You can see an example here: https://charlesrocabert.github.io/doc/evo2sim_simulation_example/viewer/viewer.html.

The different tabs are described below.

### Population viewer <a name="population_viewer"></a>
This page displays the evolution of main population statistics (population size, mean genome size, mean score, ...), as well as the evolution of the trophic network (graph of the ecological interactions).

### Best lineage viewer <a name="best_lineage_viewer"></a>
This page displays the evolution of the lineage of the last best aive individual. These informations are the most representative of evolutionary dynamics, since they contains all the mutations fixed since the beginning of the simulation.

### Best individual viewer <a name="best_individual_viewer"></a>
This page displays graphics about the last best alive individual (genome structure, genetic regulation network, metabolic network, internal metabolic state, ...).

### Environment viewer <a name="environment_viewer"></a>
This page displays the evolution of main environmental statistics, as well as its current state.

### Phylogeny viewer <a name="phylogeny_viewer"></a>
This page displays the current phylogenetic tree, as well as some phylogenetic properties through time (number of nodes, common ancestor age, ...).

## Examples <a name="examples"></a>

To test the following simulation examples, please download the attached packages. They contain simulation backups and the associated code version.

Then compile the software, create the simulation from the parameters files, or simply run it from backup files (see the User Manual provided in the package for detailed help).

You can track evolution on the fly thanks to the HTML viewer, and to the script `track_cell.py`, that displays the internal dynamics of a selected individual on the grid.

### Evolution of a stable polymorphism <a name="stable_polymorphism"></a>

In this example, a population is evolved in a periodic environment mimicking a batch culture setup. A stable polymorphism emerges, where one ecotype feeds on the exogenous food and releases by-products, while a second ecotype feeds on the by-product. Thanks to the seasonality of the environment, this interaction is negative frequency-dependent, and stable.

<a href="http://evoevo.liris.cnrs.fr/download/7_-_software/POLYMORPHISM_EXAMPLE.zip">Download the example</a>.

### Lactose-operon-like regulation <a name="lactose_operon"></a>

In this example, a population is initialized with a predefined genome, encoding for specific genetic regulation and metabolic networks. Due to strong energy trade-offs, the regulation of proteins expression is maintained for thousands of generations.

<a href="http://evoevo.liris.cnrs.fr/download/7_-_software/REGULATION_EXAMPLE.zip">Download the example</a>.

