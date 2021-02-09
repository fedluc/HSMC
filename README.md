# HSMC

hsmc performs Monte Carlo simulations of mono-disperse hard-sphere systems. The default behavior of the code is to perform NVT simulations for which it is possible to monitor
a number of observables including pressure, chemical potential, radial distribution function and order parameters. The code can also perform simulations in the isobaric ensemble (NpT), cavity simulations as described by [Torrie and Patey](https://www.tandfonline.com/doi/abs/10.1080/00268977700102821) and simulations involving the cluster moves described by [Krauth](https://arxiv.org/abs/cond-mat/0311623). More information about the observables that can be monitored in the course of such simulations are provided in the section "Input file".

## Compiling

The code can be compiled with the [make file](src/Makefile) provided in the source directory. Please note that in order to correctly compile the program it is necessary 
that the [GNU scientific library](https://www.gnu.org/software/gsl/) and the [zlib library](https://zlib.net/) are installed. 


## Running

In order to run the program it is necessary to create an input file that is passed through the command line via the argument `-i`. More information about the input file 
are given below. A full list of the accepted command line arguments can be retrieved via the command line argument `--help`.

## Input file

An example of input file can be retrieved via the command line argument `-e`. In the input file a combination of key-words and numerical values are used in order to specify 
the type of simulation that should be performed. Empty lines are skipped, lines starting with # are treated as comments. For the lines that are neither empty nor comments 
it is assumed that 

* Each line contains at most one key-word, if this is not the case then only the first key-word in the line is read. 
* When one key-word is specified, then all the numerical valus associated to that key-word are also specified. i.e. it is not possible to use the default value only for some
  numerical values.
* All the input values are given in units in which the length are normalized to the hard-sphere diameter and the energy to k_BT, where k_B is the Boltzmann's constant and 
T is the temperature.

The key-words and the numerical inputs which are expected for each key-word are listed below:

* **Cavity simulations**: Cavity simulations as described by [Torrie and Patey](https://www.tandfonline.com/doi/abs/10.1080/00268977700102821) can be performed if the 
                          keyword `cavity` is used. This keyword requires five numerical values: (1) The probability of moving a cavity, (2) the maximum 
                          separation distance between the cavities, (3) the minimum separation distance between the cavities, (3) the saving interval (in number of sweeps)                             and (5) the resolution used to write the cavity interaction potential to file. The samples for the cavity separation are written to the file 
                          *cavity_distance.dat* while the interaction potential between the cavities is written to *cavity_psi.dat*. For the cavity simulations the                                     calculations of pressure, radial distribution function, chemical potential and order parameter are not available.
                          Default: `cavity 0 1.2 0.0 100 0.01`.

* **Chemical potential**: The chemical potential can be computed via the Widom insertion method via the keyword `widom` which requires two numerical values: (1) The number
                          of insertions that have to be performed each time the chemical potential is sampled and (2) the sampling interval. The output is written to the 
                          file *chem_pot.dat*. Default: `widom 100 0`.

* **Cluster moves**: Monte Carlo simulations with cluster moves as described by [Krauth](https://arxiv.org/abs/cond-mat/0311623) can be performed if the keyword `cluster` is 
                     used. This keyword requires three numerical values: (1) The activation flag (0 or 1), (2) the number of cluster moves which constitute a sweep and (3)
                     the number of standard (i.e. local, non cluster) moves which have to be performed at the beginning of the simulation in order to perturb the perfect 
                     crystal used in the initialization. For these simulations is it possible to collect samples for all the observables which can be monitored in the course
                     of standard NVT simulations. Default: `cluster 0 1 10000`
                     
* **Configuration**: The keyword `config_write` can be used to write the current configuration to file. `config_write` requires two numerical values: (1) the saving interval                      (in number of sweeps) and (2) the number of configurations to write per file. The configurations are written to files called *"config_%.dat.gz* where %                        is a numerical identifier for the file. Default `config_write 0 100`          

* **Density**: The density is specified via the keyword `rho` followed by a numerical value for the density normalized.
                For NpT simulations this value of the density is only used to assign the total number of particles. Default: `rho 0.5`.

* **Isobaric calculations**: Isobaric calculations in the NpT ensemble can be performed if the keyword `npt` is specified. If `npt` is not specified then the code runs a
                              a standard isothermal simulation in the NVT ensemble. The keyword `npt` requires two numerical values: (1) the pressure and (2) the maximum 
                              volume deformation. In the case the optimizer is used (see section "Optimizer") the volume deformation specified with `npt` will only be used 
                              as an initial guess to the optimization procedure. The density sampled in the course of NpT simulations is written on screen. In the course of
                              NpT simulations it is possible to monitor the pressure (via the thermodynamic route) and the order parameter.
                              Default: `npt 0 0.001`.                 
                   

* **Maximum particle displacement**: The maximum particle displacement in the course of standard Monte Carlo moves is specified with `dr_max`. In the
                                     case the optimizer is used (see section "Optimizer") the value speciefied with `opt` will only be used as an initial guess 
                                     for the optimization procedure. Default `dr_max 0.05`
                                    
* **Neighbor list**: For efficient calculations of the particle's interactions, the code uses a neighbor list based on cell lists that scales line O(N), with N being the 
                     number of particles. The cell list is speficied via the keyword `neigh_list` which requires two numerical values: (1) The minimum size of the cells
                     (sizes smaller than 1.0 are not allowed) and (2) The maximum number of particles which can be expected in each cell (a safe choice is typically 5 or 10).
                     Default: `neigh_list 1.0 10`
 
* **Optimizer**: The optimizer is used to optimize the particle and volume displacements in order to obtain a target acceptance ratio (which is typically 50%). The optimizer
                  is set with the keyword `opt` which requires five numerical values: (1) An activation flag which can be 0 or 1, (2) the number of sweeps (1 sweep = N moves)
                  used for the optimization, (3) the number of samples collected during optimization (<= the number of sweeps), (4) the target acceptance ratio for particle
                  moves and (5) the target acceptance ratio for volume moves (which must be specified also if `npt` is not used). Default: `opt 1 1000 10 0.5 0.5` 
                  
* **Output on screen**: The keyword `out` can be used to specify how often to print on screen the information regarding the current sweep. `out` requires one numerical value
                        which specifies how often to print on screen. Default: `out 0`

* **Order parameter**: The order parameter (as defined by [Lechner and Dellago](https://aip.scitation.org/doi/full/10.1063/1.2977970?casa_token=Pq7x6pG7HZAAAAAA%3ASDGRrjz3OL1_tOC1qLBvvrGSDBdJEBMMDUxZcJSoyTAOEBNouzPmfF23Z25dh7R4D91fsr_0dEJw)) can be computed with the keyword `ql` which requires three 
                        numerical values: (1) The order of the order parameter, (2) The cutoff distance for the neighbor interaction (typically 1.5 is a good choice) and (3)
                        the saving interval (in number of sweeps). This keyword compute an average order parameter for the whole system by averaging over the order parameters
                        of all the particles inside the system. The output is written to the file *oder_param.dat*. Default: `ql 6 1.5 0`.

* **Pressure (thermodynamic)**: The pressure calculated via the thermodynamic definition (i.e. via virtual volume compressions) can be computed with the keyword                                                `press_thermo` which requires three numerical values: (1) The resolution for volume perturbations, (2) maximimum virtual compression
                                  and (3) the sampling interval (in number of sweeps). The output is written to the files *press_thermo.dat* and *density.dat* (the latter
                                  only for isobaric calcualtions). Default: `press_thermo 0.0001 0.002 0`.

* **Pressure (virial)**: The pressure calculated via the virial equation of state (i.e via extrapolating the radial distribution function (RDF) to contact) can be computed                            with the keyword `press_virial` which requires two numerical values: (1) the resolution for the RDF calculation and (2) the sampling intervals (in                             number of sweeps). This calcuation is available only for NVT simulations. The output is written to the file *press_virial.dat*. 
                          Default: `press_virial 0.01 0`.

* **Radial distribution function**: The radial distribution function (RDF) can be computed with the keyword `rdf` which requires four numerical values: (1) The resolution
                                    used in the RDF calculation, (2) the cutoff for the RDF calculation (should not be larger than half the smallest simulation box size), 
                                    (3) the saving interval (in number of sweeps) and (4) the number of samples per output file. All the RDF samples are written to a series                                       of files called *rdf_%.dat.gz* (the total number of files depends on the number of samples that are collected and on the parameter 4 of                                       the keyword). Each sample contains the information extracted from one configuration, so it is advisable to average the RDF samples before 
                                    using them for any sort of calculation. This calcuation is available only for NVT simulations. 
                                    WARNING: The calculation of the RDF can have a significant impact on performance since it scales like O(N^2). Default: `rdf 0.01 0 10 100`

* **Random number generator**: The random number generator used by hsmc is the [Marsenne Twister implemented in the GNU Scientific Library](https://www.gnu.org/software/gsl/doc/html/rng.html). The see used to initialize the random number generator can be specified with the keyword `seed`. Default: `seed 0`. 

* **Restart files**: The keywords `restart_write` and `restart_read` can be used to specify how to write and read restart files. `restart_write` requires one  numerical value                      which is  the interval (in number of sweeps) after which a restart file is written. The restart files are in binary format and have the name                                  *restart_%.bin*, where % is the sweep at which the file is written.  `restart_read` can be used to read a restart file and start a simulation from where                      another previous simulations had left off. If the optmizer is invoked it might modify the maximum displacements read from the restart file and, hence, it                      might not allow a perfect restart. `restart_read` requires two numerical values: (1) The activation flag (0 or 1) and (2) the name of the file used to                        read the restart from (note: it should be a file previously produced with `restart_write`. Default:  `restart_write 0` (no restart file is written),                          `restart_read 0 restart_000000.bin` (no restart is performed)


* **Simulation box**: The simulation box is assumed to be a rectangular parallelepiped whose sides are specified by assign a certain amount of building blocks in each 
                      direction. The building blocks can be either simple cubic (SC) cells or face-centerd cubic (FCC) cells containing, respectively, 1 or 4 particles. 
                      The type of building block specifies also the starting configuration, i.e SC blocks will prodoce a starting configuration in which the particles are 
                      arranged on a cubic lattice while FCC blocks will prodoce a starting configuration in which the particles are arranged on a FCC lattice.
                      The number of building blocks in the x, y and z directions are specified with the key-words `cells_x`, `cells_y` and `cell_z`, while the type for 
                      the building blocks is specified as `type` followed by 1 for SC blocks or by 2 for FCC blocks. Default:  `cells_x 3`,  `cells_x 3`,  `cells_x 3`,
                       `type 1`
                          

* **Sweeps**: The number of sweeps to perform can be specified with the command `sweep_eq` and `sweep_stat`. Both these commands require one numerical value which specifies                 the number of sweeps to perform. `sweep_eq` is used to specify the number of sweeps for equilibration during which no quantity is sampled (except for the                     density in the course of isobaric simulations) and `sweep_stat` is used to specify the number of steps used to collect the statistics for the desired quantities               Default: `sweep_eq 1000`, `sweep_stat 1000`
             
