import os, sys
import numpy as np
import matplotlib.pyplot as plt
import glob

def blocking_std(data, plt_flag=True, print_flag=True, jk_flag=False):

    # Total length of the data
    nn = len(data)

    # Use only a fraction of the data if nn is not a power of 2
    maxp2 = nn.bit_length() - 1
    nn_max = 2**(maxp2)
    if (nn_max < nn):
        print("WARNING: Used %d samples out of %d for variance estimation" % (nn_max,nn))
        nn = nn_max
        data = data[:nn]

    # Perform blocking transformations defined by Flyvberg
    sigma = np.zeros((maxp2-1,2))
    ii = 0
    while ii < maxp2-1:

        # Estimate standard deviation (with jackknife if necessary)
        if jk_flag: 
            sigma[ii,0] = jackknife(data)[1]
        else:
            sigma[ii,0] = np.std(data)/np.sqrt(nn-1)
            
        
        # Estimate uncertainty of standard deviation
        sigma[ii,1] = sigma[ii,0]/np.sqrt(2*(nn-1))
        
        # Decimation
        data = np.mean(data.reshape(-1, 2), axis=1)
        nn = len(data)
    
        # Update counter
        ii+=1

    # Print results of blocking on screen
    if print_flag:
        blocking_print(sigma)
        
    # Plot result of blocking
    if plt_flag:
        blocking_plot(sigma)

def blocking_print(sigma):

    print('Transformation, variance')
    nn_trans = len(sigma[:,0])
    for ii in range(nn_trans):
        print('%d %.8f' % (ii, sigma[ii,0]))

def blocking_plot(sigma):

    nn_trans = len(sigma[:,0])
    yerr = sigma[:,1]
    plt.errorbar(np.arange(nn_trans),sigma[:,0],yerr,
                 color='b',marker='o',linestyle='none')
    plt.ylabel('Standard deviation')
    plt.xlabel('Number of applied blocking transformations')
    plt.show()


def acf(data, norm_flag=True):

    # Maximum lag equal to length of data
    nn = len(data)

    # Initialize correlation
    corr = np.zeros(nn)
    
    # Mean of the input data
    mu = np.mean(data)

    # Compute autocorrelation
    for ii in range(nn):
        for jj in range(nn-ii):
            corr[ii] += (data[jj]-mu)*(data[ii+jj]-mu)/nn
            
    # Normalize
    if norm_flag:
        corr /= corr[0]
  
    # Output
    return corr


def acf_fft(data, norm_flag=True):

    # Compute autocorrelatiion via FFT 
    nn = next_pow_two(len(data))
    data_ft = np.fft.fft(data - np.mean(data), n=2*nn)
    corr = np.fft.ifft(data_ft * np.conjugate(data_ft))[:len(data)].real
    corr /= len(data)

    # Normalize
    if norm_flag:
        corr /= corr[0]

    return corr

def next_pow_two(n):
    i = 1
    while i < n:
        i = i << 1
    return i

def jackknife(data):

    nn = len(data)
    xi = np.zeros(nn)
    
    # Mean by leaving out one sample
    for ii in range(nn):
        xi[ii] = np.mean(np.delete(data,ii))

    # jackknife estimate for the mean
    xi_ave = np.mean(xi)

    # jackknife estimate for the variance
    xi_std = np.std(xi) * np.sqrt(nn-1)

    # Output
    return np.array([xi_ave, xi_std])







