#!/home/x_fedlu/myPython/bin/python2.7 -u

import hsmc_pressure as hsp

hsp.pressure_virial(data_dir="../tests",samples_block=2000)
hsp.pressure_thermo(data_dir="../tests",samples_block=2000)

# hsp.pressure_virial(data_dir="../tests",samples_block=10)
# hsp.pressure_virial(data_dir="../tests",samples_block=100)
# hsp.pressure_virial(data_dir="../tests",samples_block=1000)
# hsp.pressure_virial(data_dir="../tests",samples_block=5000)
