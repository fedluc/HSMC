import os, sys
import numpy as np
import math
import matplotlib.pyplot as plt
import glob
import gzip
import hsmc_stat

# ------ Compute average radial distribution function ------

def rdf(data_dir):

    # Lines for the header of each sample (assumed fixed)
    lines_header = 7

    # Get names of files in data directory
    file_names = glob.glob(os.path.join(data_dir,'rdf_*.dat.gz'))
    out_dir = data_dir
    n_files = len(file_names)
    if n_files == 0:
        sys.exit('hsmc_rdf.rdf: No rdf file was found')
    
    # Compute average rdf from HSMC output
    rdf = read_rdf_output(file_names,lines_header)

    # Plot
    plt.plot(rdf[:,0],rdf[:,1],'b')
    plt.ylabel('g(x)')
    plt.xlabel('x = r/sigma')
    #plt.show()

# ------ Read  HSMC output for the radial distribution function -----

def read_rdf_output(file_names,lines_header):

    # Global variable with lines read from file
    global lines_data
    
    # Number of samples 
    n_samples = 0.0

    # Read files listed in file_names
    for name in file_names:

        print(name)

        # Open file, read all lines and close
        with gzip.open(name) as ff:
            lines_data = ff.readlines()
        lines_file = len(lines_data)

        # Extract number of bins, volume and number of particles
        [n_bins, vol, n_part]  = [float(kk) for kk in lines_data[3].split()]
        n_bins = int(n_bins)

        # Initialize rdf
        rdf = np.zeros((n_bins,2))

        # Number of lines per sample and per block
        lines_sample = n_bins + lines_header
       
        # Read file
        lines_read = 0        
        while lines_read <= lines_file - lines_sample:
            
            # Read one sample
            rdf_tmp = read_rdf_sample(lines_read+lines_header, n_bins)

            # Update sample counter
            n_samples += 1.0

            # Update running average
            rdf[:,1] += rdf_tmp[:,1]

            # Update the number of read lines
            lines_read += lines_sample
            
            
    # Normalize the rdf
    rdf_file[:,1] /= n_samples

    # Define the set of interparticle distances
    rdf[:,0] = rdf_tmp[:,0]
            
    # Output
    return rdf


def read_rdf_sample(line_start, n_bins):

    rdf = np.zeros((n_bins,2))

    for ii in range(n_bins):
        rdf_tmp = np.array([float(kk) for kk in lines_data[ii+line_start].split()])
        rdf[ii,0] = rdf_tmp[0]
        rdf[ii,1] = rdf_tmp[1] 
    
        
    return rdf
