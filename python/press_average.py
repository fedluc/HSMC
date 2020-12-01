#!/home/x_fedlu/myPython/bin/python2.7 -u

import hsmc_pressure as hsp

hsp.pressure_virial(data_dir="../tests",samples_block=1)
hsp.pressure_virial(data_dir="../tests",samples_block=2)
hsp.pressure_virial(data_dir="../tests",samples_block=4)
hsp.pressure_virial(data_dir="../tests",samples_block=50)
hsp.pressure_virial(data_dir="../tests",samples_block=100)
hsp.pressure_virial(data_dir="../tests",samples_block=200)
hsp.pressure_virial(data_dir="../tests",samples_block=500)
hsp.pressure_virial(data_dir="../tests",samples_block=1000)
hsp.pressure_virial(data_dir="../tests",samples_block=2000)
hsp.pressure_virial(data_dir="../tests",samples_block=5000)
hsp.pressure_virial(data_dir="../tests",samples_block=10000)

hsp.pressure_thermo(data_dir="../tests",samples_block=50)
hsp.pressure_thermo(data_dir="../tests",samples_block=100)
hsp.pressure_thermo(data_dir="../tests",samples_block=200)
hsp.pressure_thermo(data_dir="../tests",samples_block=500)
hsp.pressure_thermo(data_dir="../tests",samples_block=1000)
hsp.pressure_thermo(data_dir="../tests",samples_block=2000)
hsp.pressure_thermo(data_dir="../tests",samples_block=5000)
hsp.pressure_thermo(data_dir="../tests",samples_block=10000)



hsp.pressure_virial(data_dir="../tests",samples_block=40000)
hsp.pressure_thermo(data_dir="../tests",samples_block=40000)
