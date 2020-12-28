import os, sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import hsmc_stat

# ---------------------------------------------------------
# Note: these functions analyze the chemical potential
# output produced when hsmc is run with the widom option
# ---------------------------------------------------------

# ------ Average chemical potential  -------

def ave(data_dir,file_id='chem_pot.dat',file_comments='#',
        output=True, plot=True):

    # Read data
    data = read_hsmc_output(data_dir, file_id, file_comments)

    # Average chemical potential
    mu_ave = np.average(data[:,0])

    # Average percentage  of accepted moves 
    moves_ave = np.average(data[:,1])*100

    # Print result on screen
    if output:
        print("Average chemical potential: %.8f" % mu_ave)
        print("Percentage of accepted moves: %.8f" % moves_ave)

    # Plot 
    if plot:
        dens_ave_plot = np.full(len(data[:,0]), mu_ave)    
        plt.plot(data[:,0],'b')
        plt.plot(dens_ave_plot,'r--')
        plt.ylabel('Chemical potential [k_BT]')
        plt.xlabel('Sweep number')
        plt.show()
        
    # Output
    return mu_ave


# ------ Standard deviation (via Jackknife)  -------

def std(data_dir,file_id='chem_pot.dat',file_comments='#',
        output=True):

    # Read data
    data = read_hsmc_output(data_dir, file_id, file_comments)

    # Standard deviation of the chemical potential (via jackknife)
    mu_std = hsmc_stat.jackknife(data[:,0])[1]

    # Print result on screen
    if output:
        print("Standard deviation of chemical potential: %.8f" % mu_std)
        
    # Output
    return mu_std


# ------ Histograms ------

def hist(data_dir,file_id='chem_pot.dat',file_comments='#',
        hist_bins=50,output=True, plot=True):

    # Read data
    data = read_hsmc_output(data_dir, file_id, file_comments)
    
    # Chemical potential histograms
    mu_hist = np.histogram(data[:,0], bins=hist_bins)

    # Plot
    plt.hist(data[:,0], bins=hist_bins)
    plt.xlabel('Chemical potential [k_BT]')
    plt.show()
    
    # Output
    return mu_hist


# ------ Read density output of hsmc ------

def read_hsmc_output(data_dir,file_id,file_comments):

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

    # Remove all lines where no move was accepted
    mask = data[:,1] > 0
    data = data[mask,:]

    # Output
    return data











