#!/home/x_fedlu/myPython/bin/python2.7 -u

import hsmc_pressure as hsp

hsp.pressure_virial(data_dir="../tests/nvt_simulation",samples_block=20000)
hsp.pressure_thermo(data_dir="../tests/nvt_simulation",samples_block=20000)
