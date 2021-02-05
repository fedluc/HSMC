# HSMC

hsmc performs Monte Carlo simulations of mono-disperse hard-sphere systems.

## Compiling

The code can be compiled with the [make file](src/Makefile) provided in the source directory. Please note that in order to correctly compile the program it is necessary 
that the [GNU scientific library](https://www.gnu.org/software/gsl/) is installed. 


## Running

In order to run the program it is necessary to create an input file that is passed through the command line via the argument `-i`. More information about the input file 
are given below. A full list of the accepted command line arguments can be retrieved via the command line argument `--help`.

## Input file

An example of input file can be retrieved via the command line argument `-e`. In the input file a combination of key-words and numerical values are used in order to specify 
the type of simulation that should be performed. Empty lines are skipped, lines starting with # are treated as comments. For the lines that are neither empty nor comments 
it is assumed that 

* Each line contains at most one key-word, if this is not the case then only the first key-word in the line is read. 
* All the input values are given in units in which the length are normalized to the hard-sphere diameter and the energy to k_BT, where k_B is the Boltzmann's constant and 
T is the temperature.

The key-words and the numerical inputs which are expected for each key-word are listed below:

* **Density**: The density is specified via the keyword `rho` followed by a numerical value for the density normalized.
                For NpT simulations this value of the density is only used to assign the total number of particles. Default: `rho 0.5`.

* **Isobaric calculations**: Isobaric calculations in the NpT ensemble can be performed if the command `npt` is specified. If `npt` is not specified then the code runs a
                              a standard isothermal simulation in the NVT ensemble. The command `npt` requires two numerical values: (1) the pressure and (2) the maximum 
                              volume deformation. In the case the optimizer is used (see section "Optimizer") the volume deformation specified with `npt` will only be used 
                              as an initial guess to the optimization procedure. Default: `npt 0 0.001`.                 
                   

* **Maximum particle displacement**: The maximum particle displacement in the course of standard Monte Carlo moves is specified with `dr_max`. In the
                                     case the optimizer is used (see section "Optimizer") the value speciefied with `opt` will only be used as an initial guess 
                                     for the optimization procedure. Default `dr_max 0.05`
                                    
* **Neighbor list**: For efficient calculations of the particle's interactions, the code uses a neighbor list based on cell lists that scales line O(N), with N being the 
                     number of particles. The cell list is speficied via the command `neigh_list` which requires two numerical values: (1) The minimum size of the cells
                     (sizes smaller than 1.0 are not allowed) and (2) The maximum number of particles which can be expected in each cell (a safe choice is typically 5 or 10).
                     Default: `neigh_list 1.0 10`
 
* **Optimizer**: The optimizer is used to optimize the particle and volume displacements in order to obtain a target acceptance ratio (which is typically 50%). The optimizer
                  is set with the command `opt` which requires five numerical values: (1) An activation flag which can be 0 or 1, (2) the number of sweeps (1 sweep = N moves)
                  used for the optimization, (3) the number of samples collected during optimization (<= the number of sweeps), (4) the target acceptance ratio for particle
                  moves and (5) the target acceptance ratio for volume moves (which must be specified also if `npt` is not used). Default: `opt 1 1000 10 0.5 0.5` 

* **Simulation box**: The simulation box is assumed to be a rectangular parallelepiped whose sides are specified by assign a certain amount of building blocks in each 
                       direction. The building blocks can be either simple cubic (SC) cells or face-centerd cubic (FCC) cells containing, respectively, 1 or 4 particles. 
                       The type of building block specifies also the starting configuration, i.e SC blocks will prodoce a starting configuration in which the particles are 
                       arranged on a cubic lattice while FCC blocks will prodoce a starting configuration in which the particles are arranged on a FCC lattice.
                       The number of building blocks in the x, y and z directions are specified with the key-words `cells_x`, `cells_y` and `cell_z`, while the type for 
                       the building blocks is specified as `type` followed by 1 for SC blocks or by 2 for FCC blocks. Default:  `cells_x 3`,  `cells_x 3`,  `cells_x 3`,
                       `type 1`
                          


             
