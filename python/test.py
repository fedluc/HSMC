#!/home/x_fedlu/myPython/bin/python2.7 -u

import hsmc_pressure as hsp
import hsmc_acf as acf
import matplotlib.pyplot as plt
#hsp.pressure_virial(data_dir="../tests",samples_block=40000)
#hsp.pressure_thermo(data_dir="../../HSMC_tests",samples_block=2000,npt=True) 
corr1 = acf.acf_ave("../../HSMC_tests","density.dat",samples_block=500)
corr2 = acf.acf_ave("../../HSMC_tests","density.dat",samples_block=1000)
corr3 = acf.acf_ave("../../HSMC_tests","density.dat",samples_block=2500)
corr4 = acf.acf_ave("../../HSMC_tests","density.dat",samples_block=-1)
plt.plot(corr1)
plt.plot(corr2)
plt.plot(corr3)
plt.plot(corr4)
plt.ylabel('C(t)')
axes = plt.gca()
#axes.set_xlim([0,10000])
plt.show()
