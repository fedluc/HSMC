#!/home/x_fedlu/myPython/bin/python2.7 -u

import hsmc_pressure as pressure
import hsmc_density as density
import hsmc_op as op
#density.average('../../HSMC_tests')
#density.var_blocking('../../HSMC_tests')
op.ave('../../HSMC_tests')
op.std('../../HSMC_tests')
op.std('../../HSMC_tests',jackknife=True)
