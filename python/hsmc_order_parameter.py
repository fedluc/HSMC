import os, sys
import numpy as np
import glob
import hsmc_acf
import hsmc_blocking
import matplotlib.pyplot as plt

def acf(data_dir,file_id='order_param.dat',mode='fft',
        samples_block=-1,file_comments='#',plot=True):

    # Get names of files in data directory
    file_names = glob.glob(os.path.join(data_dir,file_id))
    out_dir = data_dir
    if len(file_names) == 0:
        sys.exit('hsmc_op.acf: No data file was found')

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
            sys.exit('hsmc_op.acf: Unkown mode option')

    # Plot 
    if plot:
        plt.plot(corr,'b')
        plt.ylabel('Autocorrelation function')
        plt.xlabel('Lag (sweeps)')
        plt.show()
        
    # Output
    return corr


def average(data_dir,file_id='order_param.dat',file_comments='#',
        plot=True):

    # Get names of files in data directory
    file_names = glob.glob(os.path.join(data_dir,file_id))
    out_dir = data_dir
    if len(file_names) == 0:
        sys.exit('hsmc_op.ave: No data file was found')

    # Read data
    data = read_hsmc_output(file_names, file_comments,squeeze=True)

    # Average density
    dens_ave = np.average(data)

    # Plot 
    if plot:
        dens_ave_plot = np.full(len(data), dens_ave)    
        plt.plot(data,'b')
        plt.plot(dens_ave_plot,'r--')
        plt.ylabel('Order parameter')
        plt.xlabel('Sweep number')
        plt.show()
        
    # Output
    return dens_ave

def var_blocking(data_dir,file_id='order_param.dat',file_comments='#',
                plot = True, output=True):

    # Get names of files in data directory
    file_names = glob.glob(os.path.join(data_dir,file_id))
    out_dir = data_dir
    if len(file_names) == 0:
        sys.exit('hsmc_op.var_blocking: No data file was found')

    # Read data
    data = read_hsmc_output(file_names, file_comments, squeeze=True)

    # Blocking to estimate the variance
    sigma = hsmc_blocking.estimate_var(data,plot,output)

def read_hsmc_output(file_names,file_comments,samples_block=-1,squeeze=False):

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

    # Remove one dimension for one-dimensional output
    if squeeze: 
        data = np.squeeze(data)

    # Output
    return data











