import os, sys
import numpy as np
import math
from scipy import optimize as opt
from scipy import special as spec
import glob
import gzip as gz

# ------ Compute pressure from virial route ------
def pressure_virial(data_dir=os.getcwd(),samples_block=1):

    # Lines for the header of each sample (assumed fixed)
    lines_header = 7

    # Get names of files in data directory
    file_names = glob.glob(os.path.join(data_dir,'press_virial.dat'))
    out_dir = data_dir
    n_files = len(file_names)
    if n_files == 0:
        sys.exit('hsmc_pressure.pressv: No pressure file was found')
    
    # Samples for the radial distribution function at contact
    [coeff, vol, n_part] = read_output(file_names,samples_block,lines_header,"virial")
    rdf_contact = coeff[:,0] + coeff[:,1]

    # Average and standard deviation of the pressure
    press = 1 + (2.*math.pi/3.) * n_part/vol * rdf_contact
    print(len(coeff),np.average(press),np.var(press)/np.sqrt(len(coeff)))
    

# ------ Compute pressure from virial route ------
def pressure_thermo(data_dir=os.getcwd(),samples_block=1):

    # Lines for the header of each sample (assumed fixed)
    lines_header = 7

    # Get names of files in data directory
    file_names = glob.glob(os.path.join(data_dir,'press_thermo.dat'))
    out_dir = data_dir
    n_files = len(file_names)
    if n_files == 0:
        sys.exit('hsmc_pressure.pressv: No pressure file was found')
    
    # Samples for the excess pressure
    [coeff, vol, n_part] = read_output(file_names,samples_block,lines_header,"thermo")
    p_ex = coeff[:,0]/n_part

    # Average and standard deviation of the pressure
    press = 1 + p_ex
    print(len(coeff),np.average(press),np.var(press)/np.sqrt(len(coeff)))


# ------ Read  HSMC output -----
def read_output(file_names,samples_block,lines_header,press_flag):

    # Global variables declaration
    global lines

    initFlag = True
    for name in file_names:
            
        # Open file, read all lines and close
        ff = open(name)
        lines = ff.readlines()
        ff.close()
        lines_file = len(lines)
                
        # Extract number of bins, volume and number of particles
        [n_bins, vol, n_part]  = [float(kk) for kk in lines[3].split()]
        n_bins = int(n_bins)

        # Initialize histograms and pressure
        hist = np.zeros((n_bins,2))

        # Number of lines per sample and per block
        lines_sample = n_bins + lines_header
        lines_block = samples_block*lines_sample

        # Number of usable samples in file
        n_samples_file = lines_file//lines_block
        coeff = np.zeros((n_samples_file,2))
        
        # Read file
        n_samples = 0
        n_blocks = 0;
        lines_read = 0;        
        while lines_read <= lines_file - lines_block:

            # Collect samples for average
            while n_samples < samples_block:
            
                # Read one sample
                histTmp = read_sample(lines_read+lines_header, n_bins)
                # Create vector with abscissa for fit 
                if initFlag:
                    hist[:,0] = histTmp[:,0]
                    initFlag = False
                # Update running average
                hist[:,1] += histTmp[:,1]/samples_block

                # Update number of samples
                n_samples += 1
            
                
            # Prepare function to fit
            if (press_flag == "thermo"):

                # Remove indeterminate forms
                del_idx = np.argwhere(hist[:,1]==0)
                hist = np.delete(hist,del_idx,axis=0)
                # Create function to fit
                hist[:,1] = -np.log(hist[:,1])/(hist[:,0])

            # Fit data
            coeff[n_blocks,:] = data_fit(hist)

            # Update number of blocks 
            n_blocks+=1

            # Update the number of read lines
            lines_read += lines_block
                     
            # Reset histogram
            hist = np.zeros((n_bins,2))
            initFlag = True

            # Reset number of samples
            n_samples = 0

            
        # Remove nan values 
        del_idx = np.argwhere(np.isnan(coeff))
        coeff = np.delete(coeff,del_idx,axis=0)

        # Output
        return [coeff, vol, n_part]



def read_sample(line_start, n_bins):

    hist = np.zeros((n_bins,2))

    for jj in range(n_bins):
        histTmp = np.array([float(kk) for kk in lines[jj+line_start].split()])
        hist[jj,0] = histTmp[0]
        hist[jj,1] = histTmp[1] 
    
    return hist

def linear_fit(x, b0, b1):
    return b0 + b1*x

def data_fit(hist):

    if (np.shape(hist)[0] >= 2):
        [fit_par, fit_par_cov] = opt.curve_fit(linear_fit,
                                               hist[:,0], hist[:,1],
                                               p0=[10.0,1.0])
    else:

        fit_par = [np.nan, np.nan]

    
    return fit_par
