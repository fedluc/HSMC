import os, sys
import numpy as np
import matplotlib.pyplot as plt

def estimate_var(data, plot=True, output=True):

    # Total length of the data
    nn_tot = len(data)

    # Use only a fraction of the data if nn_tot is not a power of 2
    maxp2 = nn_tot.bit_length() - 1
    nn_max = 2**(maxp2)
    if (nn_max < nn_tot):
        print("WARNING: Used %d samples out of %d for variance estimation" % (nn_max,nn_tot))
        nn_tot = nn_max
        data = data[:nn_tot]

    # Perform blocking transformations defined by Flyvberg
    sigma = np.zeros((maxp2-1,2))
    ii = 0
    while ii < maxp2-1:

        # Estimate standard deviation
        sigma[ii,0] = np.std(data)/np.sqrt(nn_tot-1)
        
        # Estimate uncertainty of standard deviation
        sigma[ii,1] = sigma[ii,0]/np.sqrt(2*(nn_tot-1))
        
        # Decimation
        data = np.mean(data.reshape(-1, 2), axis=1)
        nn_tot = len(data)
    
        # Update counter
        ii+=1

    # Print results of blocking on screen
    output_var(sigma)
        
    # Plot result of blocking
    plot_var(sigma)

def output_var(sigma):

    print('Transformation, variance')
    nn_trans = len(sigma[:,0])
    for ii in range(nn_trans):
        print('%d %.8f' % (ii, sigma[ii,0]))

def plot_var(sigma):

    nn_trans = len(sigma[:,0])
    yerr = sigma[:,1]
    plt.errorbar(np.arange(nn_trans),sigma[:,0],yerr,
                 color='b',marker='o',linestyle='none')
    plt.ylabel('Standard deviation')
    plt.xlabel('Number of applied blocking transformations')
    plt.show()
