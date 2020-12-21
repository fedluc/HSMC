#!/home/x_fedlu/myPython/bin/python2.7 -u

import hsmc_density as density
import hsmc_op as op
import hsmc_pressure as press

# NVT simulations
# data_dir = '../../HSMC_tests'
# op.ave(data_dir)
# op.std(data_dir)
# op.hist(data_dir,hist_bins=10)
# op.acf(data_dir,samples_block=10)
# press.pressure_virial(data_dir,samples_block=128)
# press.pressure_thermo(data_dir,samples_block=128)

# NpT simulations
data_dir = '../../HSMC_tests'
op.ave(data_dir)
op.std(data_dir)
op.hist(data_dir)
op.acf(data_dir)
density.ave(data_dir)
density.std(data_dir)
density.hist(data_dir)
density.acf(data_dir)
press.pressure_thermo(data_dir,samples_block=128)
