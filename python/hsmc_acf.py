import os, sys
import numpy as np
import glob

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











