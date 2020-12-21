import os, sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import hsmc_stat

# ---------------------------------------------------------
# Note: these functions are only applicable to the density
# output produced when hsmc runs in the NpT ensemble
# ---------------------------------------------------------

# ------ Average density -------

def ave(data_dir,file_id='density.dat',file_comments='#',
        output=True, plot=True):

    # Read data
    data = read_hsmc_output(data_dir, file_id, file_comments,squeeze=True)

    # Average density
    dens_ave = np.average(data)
    
    # Print result on screen
    if output:
        print("Average density: %.8f" % dens_ave)

    # Plot 
    if plot:
        dens_ave_plot = np.full(len(data), dens_ave)    
        plt.plot(data,'b')
        plt.plot(dens_ave_plot,'r--')
        plt.ylabel('Density [1/sigma^3]')
        plt.xlabel('Sweep number')
        plt.show()
        
    # Output
    return dens_ave


# ------ Histograms ------

def hist(data_dir,file_id='density.dat',file_comments='#',
        hist_bins=100,output=True, plot=True):

    # Read data
    data = read_hsmc_output(data_dir, file_id, file_comments,squeeze=True)
    
    # Density histograms
    dens_hist = np.histogram(data, bins=hist_bins)

    # Plot
    plt.hist(data, bins=hist_bins)
    plt.xlabel('Density [1/sigma^3]')
    plt.show()
    
    # Output
    return dens_hist

# ------ Standard deviation of the density (via blocking)  -------

def std(data_dir,file_id='density.dat',file_comments='#',
        output=True, plot=True, jackknife=False):

    # Read data
    data = read_hsmc_output(data_dir, file_id, file_comments, squeeze=True)

    # Blocking to estimate the variance
    sigma = hsmc_stat.blocking_std(data,plt_flag=plot,
                                   print_flag=output,jk_flag=jackknife)

# ------ Density autocorrelation function  -------

def acf(data_dir,file_id='density.dat',mode='fft',
        samples_block=-1,file_comments='#',plot=True):

    # Read data
    data = read_hsmc_output(data_dir, file_id, file_comments, samples_block)
    n_blocks = data.shape[0]

    # Compute average autocorrelation function
    corr = np.zeros(data.shape[1])
    for ii in range(n_blocks):
        if mode == 'fft':
            corr += hsmc_stat.acf_fft(data[ii,:])/n_blocks
        elif mode == 'dsum':
            corr += hsmc_stat.acf(data[ii,:])/n_blocks
        else:
            sys.exit('hsmc_density.acf: Unkown mode option')

    # Plot 
    if plot:
        plt.plot(corr,'b')
        plt.ylabel('Autocorrelation function')
        plt.xlabel('Lag (sweeps)')
        plt.show()
        
    # Output
    return corr


# ------ Read density output of hsmc ------

def read_hsmc_output(data_dir,file_id,file_comments,samples_block=-1,squeeze=False):

    # Get names of files in data directory
    file_names = glob.glob(os.path.join(data_dir,file_id))
    out_dir = data_dir
    if len(file_names) == 0:
        sys.exit('hsmc_density.read_hsmc_output: No data file was found')

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











