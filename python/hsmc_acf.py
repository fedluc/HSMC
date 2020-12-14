import os, sys
import numpy as np
import glob

def acf_ave(data_dir,file_id,mode='fft',samples_block=-1,file_comments='#'):

    # Get names of files in data directory
    file_names = glob.glob(os.path.join(data_dir,file_id))
    out_dir = data_dir
    if len(file_names) == 0:
        sys.exit('hsmc_acf.pressv: No data file was found')

    # Read data
    data = read_hsmc_output(file_names, file_comments, samples_block)
    n_blocks = data.shape[0]

    # Compute average autocorrelation function
    corr = np.zeros(data.shape[1])
    for ii in range(n_blocks):
        if mode == 'fft':
            corr += acf_fft(data[ii,:])/n_blocks
        elif mode == 'dsum':
            corr += acf(data[ii,:])/n_blocks
        else:
            sys.exit('Unkown mode option')
            
    # Output
    return corr

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

def read_hsmc_output(file_names,file_comments,samples_block=-1):

    # Read files listed in file_names
    init_file_flag = True
    for name in file_names:

        # Open file, read all lines and close
        data_file = np.loadtxt(name, comments=file_comments)

        # Update output
        if init_file_flag:
            data = data_file
            init_file_flag = False
        else:
            data = np.append(data,data_file,axis=0)


    # Re-shape the data in blocks containing samples_blocks elements
    data_samples = len(data)
    if samples_block > data_samples:
        print("WARNING: Not enough data to fill one block, samples_block limited to %d" % data_samples)
        samples_block = data_samples
    elif samples_block < 0:
        samples_block = data_samples
    n_samples = data_samples - (data_samples % samples_block)
    data = data[:n_samples]
    data = data.reshape(-1, samples_block)

    # Output
    return data











