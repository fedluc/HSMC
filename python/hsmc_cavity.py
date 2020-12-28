import os, sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import hsmc_stat

# ---------------------------------------------------------
# Note: these functions analyze the cavity output produced
# when hsmc is run with the cavity option
# ---------------------------------------------------------

# ------ Histograms ------

def hist(data_dir,file_id='cavity_distance.dat',file_comments='#',
        hist_bins=50,output=True, plot=True):

    # Read data
    data = read_hsmc_output(data_dir, file_id, file_comments)
    
    # Chemical potential histograms
    cavity_hist = np.histogram(data, bins=hist_bins)

    # Plot
    plt.hist(data, bins=hist_bins)
    plt.xlabel('Cavity distance [sigma]')
    plt.show()
    
    # Output
    return cavity_hist


# ------ Read density output of hsmc ------

def read_hsmc_output(data_dir,file_id,file_comments):

    # Get names of files in data directory
    file_names = glob.glob(os.path.join(data_dir,file_id))
    out_dir = data_dir
    if len(file_names) == 0:
        sys.exit('hsmc_cavity.read_hsmc_output: No data file was found')

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

    # Output
    return data











