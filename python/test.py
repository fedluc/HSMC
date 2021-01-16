#!/home/x_fedlu/myPython/bin/python2.7 -u

import hsmc_density as density
import hsmc_op as op
import hsmc_pressure as press
import hsmc_chem_pot as mu
import hsmc_cavity as cavity
import hsmc_rdf as rdf

# NVT simulations
# data_dir = '../../HSMC_tests'
# op.ave(data_dir)
# op.std(data_dir)
# op.hist(data_dir,hist_bins=10)
# op.acf(data_dir,samples_block=10)
# press.pressure_virial(data_dir,samples_block=128)
# press.pressure_thermo(data_dir,samples_block=128)

# # NpT simulations
# data_dir = '../../HSMC_tests'
# op.ave(data_dir)
# op.std(data_dir)
# op.hist(data_dir)
# op.acf(data_dir)
# density.ave(data_dir)
# density.std(data_dir)
# density.hist(data_dir)
# density.acf(data_dir)
# press.pressure_thermo(data_dir,samples_block=128)

# NVT simulations with Widom insertion
# data_dir = '../../HSMC_tests'
# mu.ave(data_dir,plot=False)
# mu.hist(data_dir)

# cavity NVT simulations
# data_dir = '../../HSMC_tests'
# cavity.rdf(data_dir,normalize=True)
# cavity.lny(data_dir)
# cavity.psi(data_dir)

# Radial distribution function
data_dir = '../../HSMC_tests'
rdf.rdf(data_dir)
