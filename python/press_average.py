#!/home/x_fedlu/myPython/bin/python2.7 -u

import hsmc_pressure as hsp

hsp.pressure_virial(data_dir="../../HSMC_tests/nvt_simulation_1",samples_block=2000)
hsp.pressure_thermo(data_dir="../../HSMC_tests/nvt_simulation_1",samples_block=2000)
hsp.pressure_virial(data_dir="../../HSMC_tests/nvt_simulation_2",samples_block=2000)
hsp.pressure_thermo(data_dir="../../HSMC_tests/nvt_simulation_2",samples_block=2000)
# hsp.pressure_thermo(data_dir="../tests",samples_block=2000,npt=True)
